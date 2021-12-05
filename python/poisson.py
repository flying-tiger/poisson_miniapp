#!/usr/bin/env python3
import sys
import time
import numpy as np
import numba as nb

class Metrics:
    def __init__(self, nx, ny, nz):
        self.dx = 1.0 / (nx-1)
        self.dy = 1.0 / (ny-1)
        self.dz = 1.0 / (nz-1)
        self.ax = 1.0 / (self.dx*self.dx)
        self.ay = 1.0 / (self.dy*self.dy)
        self.az = 1.0 / (self.dz*self.dz)
        self.dxyz = 2.0 * (self.ax + self.ay + self.az)
        self.inv_dxyz = 1.0/self.dxyz

def simple_update(u0, u1, f, ax, ay, az, inv_dxyz):
    norm_du = 0.0
    for i in nb.prange(1, u0.shape[0]-1): # prange => parallel
      for j in range(1, u0.shape[1]-1):
        for k in range(1, u0.shape[2]-1):
          u1[i,j,k] = -f[i,j,k]
          u1[i,j,k] += ax*u0[i-1, j, k]
          u1[i,j,k] += ax*u0[i+1, j, k]
          u1[i,j,k] += ay*u0[i, j-1, k]
          u1[i,j,k] += ay*u0[i, j+1, k]
          u1[i,j,k] += az*u0[i, j, k-1]
          u1[i,j,k] += az*u0[i, j, k+1]
          u1[i,j,k] *= inv_dxyz
          du = u1[i,j,k] - u0[i,j,k]
          norm_du += du*du
    return np.sqrt(norm_du)

def vectorized_update(u0, u1, f, ax, ay, az, inv_dxyz):
    ni, nj, nk = u0.shape
    i = slice(1,ni-1)
    j = slice(1,nj-1)
    k = slice(1,nk-1)
    u1[i,j,k] = -f[i,j,k]
    u1[i,j,k] += ax*u0[:ni-2,j,k]
    u1[i,j,k] += ax*u0[2:,j,k]
    u1[i,j,k] += ay*u0[i,:nj-2,k]
    u1[i,j,k] += ay*u0[i,2:,k]
    u1[i,j,k] += az*u0[i,j,:nk-2]
    u1[i,j,k] += az*u0[i,j,2:]
    u1[i,j,k] *= inv_dxyz
    # u1=u0 on boundary so don't need to exclude from norm()
    # Need to ravel so numba will use 1D c-impl of norm()
    return np.linalg.norm(np.ravel(u1 - u0))

@nb.stencil
def poisson_stencil(u, f, ax, ay, az, inv_dxyz):
    return inv_dxyz * (
        ax*u[-1, 0, 0] +
        ax*u[ 1, 0, 0] +
        ay*u[ 0,-1, 0] +
        ay*u[ 0, 1, 0] +
        az*u[ 0, 0,-1] +
        az*u[ 0, 0, 1] -
        f[0, 0, 0]
    )

def stencil_update(u0, u1, f, ax, ay, az, inv_dxyz):
    poisson_stencil(u0, f, ax, ay, az, inv_dxyz, out=u1)
    # u1=u0 on boundary so don't need to exclude from norm()
    # Need to ravel so numba can find 1D c-impl of norm()
    return np.linalg.norm(np.ravel(u1-u0))

def main(argv):

    # Argument Processing
    kn  = argv.pop(0)        if argv else 'vector'
    nx  = int(argv.pop(0))+2 if argv else 18 # +2 for boundary nodes
    ny  = int(argv.pop(0))+2 if argv else nx
    nz  = int(argv.pop(0))+2 if argv else nx
    nt  = int(argv.pop(0))   if argv else 1000
    npr = int(argv.pop(0))   if argv else 100

    # Allocate and initialize solution, source term
    x = np.linspace(0.0, 1.0, nx).reshape(nx,1,1)
    y = np.linspace(0.0, 1.0, ny).reshape(1,ny,1)
    z = np.linspace(0.0, 1.0, nz).reshape(1,1,nz)
    r2 = x*x + y*y + z*z
    u0 = np.exp(-r2)
    f  = (4.0*r2 - 6.0)*u0
    u0[1:nx-1, 1:ny-1, 1:nz-1] = 0.0
    u1 = np.copy(u0)

    # Iterate
    update_kernel = {
        'simple':      simple_update,
        'vector':      vectorized_update,
        'stencil':     stencil_update,
        'simple_jit':  nb.njit(simple_update),
        'vector_jit':  nb.njit(vectorized_update),
        'stencil_jit': nb.njit(stencil_update),
        'simple_par':  nb.njit(simple_update, parallel=True),
        'vector_par':  nb.njit(vectorized_update, parallel=True),
        'stencil_par': nb.njit(stencil_update, parallel=True),
    }
    metrics = Metrics(nx, ny, nz)
    mvalues = (metrics.ax, metrics.ay, metrics.az, metrics.inv_dxyz)
    start = time.perf_counter()
    for n in range(1,nt+1):
        # Pass metrics as scalars so we can jit with Numba
        norm_du = update_kernel[kn](u0, u1, f, *mvalues)
        if n % npr == 0:
            print(f"{n:6d} {norm_du/(nx*ny*nz):12.6e}")
        u0, u1 = u1, u0
    end = time.perf_counter()

    # Extract check value
    u_check = u0[nx//2, ny//2, nz//2]
    u_exact = np.exp(-r2[nx//2, ny//2, nz//2])
    u_error = np.abs(u_check - u_exact)

    # Result summary
    dt = end - start
    mups = (nt/1.0e6)*(nx-2)*(ny-2)*(nz-2)/dt
    print(f"\nCheck Value:       {u_check:8.6f}")
    print(f"Exact Value:       {u_exact:8.6f}")
    print(f"Check Value Error: {u_error:8.6f}")
    print(f"Elasped Time:      {dt:.2f} sec")
    print(f"Throughput:        {mups:.2f} MUPS\n\n")

if __name__ == '__main__':
    main(sys.argv[1:])
