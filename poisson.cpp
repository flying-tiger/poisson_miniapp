#include <cassert>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include "mdspan.hpp"
#include "stopwatch.hpp"

namespace stdex = std::experimental;
using span2d = stdex::mdspan<double, stdex::extents<
    stdex::dynamic_extent,
    stdex::dynamic_extent
>>;

struct Metrics {
    double dx, dy;
    double ax, ay, dxy;
    double inv_dxy;
    Metrics(int nx, int ny) {
        dx  = 1.0 / (nx-1);
        dy  = 1.0 / (ny-1);
        ax  = 1.0 / (dx*dx);
        ay  = 1.0 / (dy*dy);
        dxy = 0.5 * (ax + ay);
        inv_dxy = 1.0/dxy;
    }
};

double boundary_value(double x, double y) {
    double a2 = 0.1;
    double r2 = x*x + y*y;
    return a2 / (a2 + r2);
}

auto simple_update(const span2d f0, const span2d f1, const Metrics& m) {
    const int ni = f0.extent(0);
    const int nj = f0.extent(1);

    auto norm_df = 0.0 * f0(1,1);
    for (int i = 1; i < ni-1; ++i) {
    for (int j = 1; j < nj-1; ++j) {
        f1(i,j) = 0.25*m.inv_dxy*(
            m.ax*(f0(i-1,j) + f0(i+1,j)) +
            m.ay*(f0(i,j-1) + f0(i,j+1))
        );
        auto df = f1(i,j) - f0(i,j);
        norm_df += df*df;
    }}

    return std::sqrt(norm_df)/(ni*nj);
}

auto blocked_update(const span2d f0, const span2d f1, const Metrics& m) {
    const int ni = f0.extent(0);
    const int nj = f0.extent(1);

    constexpr int nb = 8;
    assert((ni-2) % nb == 0);
    assert((nj-2) % nb == 0);
    const int nbi = (ni-2) / nb;
    const int nbj = (nj-2) / nb;

    auto norm_df = 0.0 * f0(1,1);
    for (int bi = 0; bi < nbi; ++bi) {
    for (int bj = 0; bj < nbj; ++bj) {
        const int imin = bi*nb+1, imax = imin+nb;
        const int jmin = bj*nb+1, jmax = jmin+nb;
        for (int i = imin; i < imax; ++i) {
        for (int j = jmin; j < jmax; ++j) {
            f1(i,j) = 0.25*m.inv_dxy*(
                m.ax*(f0(i-1,j) + f0(i+1,j)) +
                m.ay*(f0(i,j-1) + f0(i,j+1))
            );
            auto df = f1(i,j) - f0(i,j);
            norm_df += df*df;
        }}
    }}

    return std::sqrt(norm_df)/(ni*nj);
}

void print(span2d f) {
    for (int j = 0; j < f.extent(1); ++j) {
        for (int i = 0; i < f.extent(0); ++i) {
            printf("%6.4f ", f(i,j));
        }
        printf("\n");
    }
}

int main(int argc, char* argv[]) {

    // Arugment processing
    const int nx = (argc > 1)? std::stoi(argv[1])+2 : 18;  // +2 for boundary nodes
    const int ny = (argc > 2)? std::stoi(argv[2])+2 : nx;
    const int nt = (argc > 3)? std::stoi(argv[3])   : 1000;
    const int np = (argc > 4)? std::stoi(argv[4])   : 100;

    // Allocate data arrays
    auto storage0 = std::vector<double>(nx*ny, 0.0);
    auto storage1 = std::vector<double>(nx*ny, 0.0);
    auto f0 = span2d(storage0.data(), nx, ny);
    auto f1 = span2d(storage1.data(), nx, ny);

    // Apply boundary conditions
    auto metrics = Metrics(nx,ny);
    for (int i = 0; i < nx; ++i) {
        f0(i,0)    = f1(i,0)    = boundary_value(i*metrics.dx, 0.0);
        f0(i,ny-1) = f1(i,ny-1) = boundary_value(i*metrics.dx, 1.0);
    }
    for (int j = 0; j < ny; ++j) {
        f0(0,j)    = f1(0,j)    = boundary_value(0.0, j*metrics.dy);
        f0(nx-1,j) = f1(nx-1,j) = boundary_value(1.0, j*metrics.dy);
    }

    // Iterate
    Stopwatch<> sw;
    sw.start();
    for (int n = 1; n <= nt; ++n) {
        auto norm_df = simple_update(f0, f1, metrics);
        //auto norm_df = blocked_update(f0, f1, metrics);
        //auto norm_df = stdpar_update(f0, f1, metrics);
        //auto norm_df = kokkos_update(f0, f1, metrics);
        if (n % np == 0) printf("%6i %12.6e\n", n, norm_df);
        std::swap(f1,f0);
    }
    sw.stop();

    // Final solution
    int imid = nx/2;
    int jmid = ny/2;
    auto dt = sw.get_elapsed();
    printf("\nFinal Value at (%i,%i): %6.4f\n", imid, jmid, f0(imid,jmid));
    printf("Elapsed Time:             %.2f sec\n", dt);
    printf("Throughput:               %.2f MUPS\n\n", 1.0e-6*nt*(nx-2)*(ny-2)/dt);
    return 0;

}
