#include <array>
#include <cassert>
#include <cstdio>
#include <cmath>
#include <string>
#include <map>
#include <vector>
#include "mdspan.hpp"
#include "stopwatch.hpp"

using namespace std::string_literals;
namespace stdex = std::experimental;
using Span3d = stdex::mdspan<double, stdex::extents<
    stdex::dynamic_extent,
    stdex::dynamic_extent,
    stdex::dynamic_extent
>>;
using CSpan3d = stdex::mdspan<const double, stdex::extents<
    stdex::dynamic_extent,
    stdex::dynamic_extent,
    stdex::dynamic_extent
>>;
using Double2 = std::array<double,2>;
using Int3 = std::array<int,3>;
using std::pair;

struct Metrics {
    double dx, dy, dz;
    double ax, ay, az, dxyz;
    double inv_dxyz;
    Metrics(int nx, int ny, int nz) {
        dx = 1.0 / (nx-1);
        dy = 1.0 / (ny-1);
        dz = 1.0 / (nz-1);
        ax = 1.0 / (dx*dx);
        ay = 1.0 / (dy*dy);
        az = 1.0 / (dz*dz);
        dxyz = 2.0 * (ax + ay + az);
        inv_dxyz = 1.0/dxyz;
    }
};

inline Double2 solution(int i, int j, int k, const Metrics& m) {
    auto x  = i*m.dx;
    auto y  = j*m.dy;
    auto z  = k*m.dz;
    auto r2 = x*x + y*y + z*z;
    auto u  = std::exp(-r2);
    auto f  = (4.0*r2 - 6.0)*u;
    return {u, f};
}

inline Double2 stencil(int i, int j, int k, const CSpan3d u0, const CSpan3d f, const Metrics& m) {
    double u = -f(i,j,k);
    u += m.ax*u0(i-1,j,k);
    u += m.ax*u0(i+1,j,k);
    u += m.ay*u0(i,j-1,k);
    u += m.ay*u0(i,j+1,k);
    u += m.az*u0(i,j,k-1);
    u += m.az*u0(i,j,k+1);
    u *= m.inv_dxyz;
    double du = u - u0(i,j,k);
    return {u, du};
}

double simple_update(const CSpan3d u0, const Span3d u1, const CSpan3d f, const Metrics& m) {
    double norm_du = 0.0;
    for (int i = 1; i < u0.extent(0)-1; ++i) {
    for (int j = 1; j < u0.extent(1)-1; ++j) {
    for (int k = 1; k < u0.extent(2)-1; ++k) {
        auto [u, du] = stencil(i, j, k, u0, f, m);
        u1(i,j,k) = u;
        norm_du += du*du;
    }}}
    return std::sqrt(norm_du);
}

double simple_update_omp(const CSpan3d u0, const Span3d u1, const CSpan3d f, const Metrics& m) {
    double norm_du = 0.0;
    #pragma omp parallel for collapse(2) reduction(+:norm_du)
    for (int i = 1; i < u0.extent(0)-1; ++i) {
    for (int j = 1; j < u0.extent(1)-1; ++j) {
    for (int k = 1; k < u0.extent(2)-1; ++k) {
        auto [u, du] = stencil(i, j, k, u0, f, m);
        u1(i,j,k) = u;
        norm_du += du*du;
    }}}
    return std::sqrt(norm_du);
}

double blocked_update(const CSpan3d u0, const Span3d u1, const CSpan3d f, const Metrics& m) {
    constexpr int block_size = 8;
    const int ni = u0.extent(0);
    const int nj = u0.extent(1);
    const int nk = u0.extent(2);
    assert((ni-2) % block_size == 0);
    assert((nj-2) % block_size == 0);
    double norm_du = 0.0;
    for (int ib = 1; ib < ni-1; ib+=block_size) {
    for (int jb = 1; jb < nj-1; jb+=block_size) {
        for (int i = ib; i < ib+block_size; ++i) {
        for (int j = jb; j < jb+block_size; ++j) {
        for (int k = 1;  k < nk-1;          ++k) {
            auto [u, du] = stencil(i, j, k, u0, f, m);
            u1(i,j,k) = u;
            norm_du += du*du;
        }}}
    }}
    return std::sqrt(norm_du);
}

double blocked_update_omp(const CSpan3d u0, const Span3d u1, const CSpan3d f, const Metrics& m) {
    constexpr int block_size = 8;
    const int ni = u0.extent(0);
    const int nj = u0.extent(1);
    const int nk = u0.extent(2);
    assert((ni-2) % block_size == 0);
    assert((nj-2) % block_size == 0);
    double norm_du = 0.0;
    #pragma omp parallel for collapse(2) reduction(+:norm_du)
    for (int ib = 1; ib < ni-1; ib+=block_size) {
    for (int jb = 1; jb < nj-1; jb+=block_size) {
        for (int i = ib; i < ib+block_size; ++i) {
        for (int j = jb; j < jb+block_size; ++j) {
        for (int k = 1;  k < nk-1;          ++k) {
            auto [u, du] = stencil(i, j, k, u0, f, m);
            u1(i,j,k) = u;
            norm_du += du*du;
        }}}
    }}
    return std::sqrt(norm_du);
}

int main(int argc, char* argv[]) {

    // Argument processing
    const std::string kname = (argc > 1)? argv[1] : "simple";
    const int nx = (argc > 2)? std::stoi(argv[2])+2 : 18;  // +2 for boundary nodes
    const int ny = (argc > 3)? std::stoi(argv[3])+2 : nx;
    const int nz = (argc > 4)? std::stoi(argv[4])+2 : nx;
    const int nt = (argc > 5)? std::stoi(argv[5])   : 1000;
    const int np = (argc > 6)? std::stoi(argv[6])   : 100;

    // Select update kernel
    auto kernels = std::map{
        pair{"simple"s,     simple_update},
        pair{"simple_omp"s, simple_update_omp},
        pair{"block"s,      blocked_update},
        pair{"block_omp"s,  blocked_update_omp}
    };
    assert(kernels[kname]);
    auto update_kernel = kernels[kname];

    // Allocate data arrays
    auto storage_u0 = std::vector<double>(nx*ny*nz, 0.0);
    auto storage_u1 = std::vector<double>(nx*ny*nz, 0.0);
    auto storage_f  = std::vector<double>(nx*ny*nz, 0.0);
    auto u0 = Span3d(storage_u0.data(), nx, ny, nz);
    auto u1 = Span3d(storage_u1.data(), nx, ny, nz);
    auto f  = Span3d(storage_f.data(),  nx, ny, nz);

    // Apply boundary conditions, source term
    auto metrics = Metrics(nx,ny,nz);
    for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
    for (int k = 0; k < nz; ++k) {
        const auto [value, source] = solution(i, j, k, metrics);
        bool iboundary = (i == 0) || (i == nx-1);
        bool jboundary = (j == 0) || (j == ny-1);
        bool kboundary = (k == 0) || (k == nz-1);
        if (iboundary || jboundary || kboundary) {
            u0(i,j,k) = value;
            u1(i,j,k) = value;
        }
        f(i,j,k) = source;
    }}}

    // Iterate
    Stopwatch<> sw;
    sw.start();
    for (int n = 1; n <= nt; ++n) {
        auto norm_du = update_kernel(u0, u1, f, metrics);
        if (n % np == 0) printf("%6i %12.6e\n", n, norm_du/(nx*ny*nz));
        std::swap(u1,u0);
    }
    sw.stop();

    // Extract check value
    const auto u_check = u0(nx/2, ny/2, nz/2);
    const auto [u_exact, f_exact] = solution(nx/2, ny/2, nz/2, metrics);

    // Result summary
    auto dt = sw.get_elapsed();
    printf("\nKernel:            %s\n", kname.c_str());
    printf("Check Value:       %8.6f\n", u_check);
    printf("Exact Value:       %8.6f\n", u_exact);
    printf("Check Value Error: %8.6f\n", std::abs(u_check - u_exact));
    printf("Elapsed Time:      %.2f sec\n", dt);
    printf("Throughput:        %.2f MUPS\n\n", (nt/1.0e6)*(nx-2)*(ny-2)*(nz-2)/dt);
    return 0;

}
