program poisson
    use utils
    implicit none

    real(8), allocatable, target :: data0(:,:), data1(:,:)
    real(8), pointer :: f0(:,:), f1(:,:)
    real(8) :: t0, tf, dt, df, norm_df
    type(metrics_t) :: metrics
    type(timer_t) :: timer
    integer :: narg, nx, ny, nt, np
    integer :: i, j, n

    ! Argument processing
    narg = command_argument_count()
    nx = 18;   if (narg >= 1) nx = get_arg(1)+2
    ny = nx;   if (narg >= 2) ny = get_arg(2)+2
    nt = 1000; if (narg >= 3) nt = get_arg(3)
    np = 100;  if (narg >= 4) np = get_arg(4)

    ! Allocate data arrays
    allocate(data0(nx,ny), data1(nx,ny))
    f0 => data0
    f1 => data1

    ! Apply boundary conditions
    metrics = metrics_t(nx,ny)
    do i = 1,nx
        f0(i,1)  = boundary_value((i-1)*metrics%dx, 0.0d0)
        f0(i,ny) = boundary_value((i-1)*metrics%dx, 1.0d0)
        f1(i,1)  = f0(i,1)
        f1(i,ny) = f0(i,ny)
    end do
    do j = 1,ny
        f0(1,j)  = boundary_value(0.0d0, (j-1)*metrics%dy)
        f0(nx,j) = boundary_value(1.0d0, (j-1)*metrics%dy)
        f1(1,j)  = f0(1,j)
        f1(nx,j) = f0(nx,j)
    end do

    ! Iterate
    call timer%tick()
    do n = 1,nt
        norm_df = 0.0d0
        do j = 2,ny-1
        do i = 2,nx-1
            f1(i,j) = 0.25d0*metrics%inv_dxy*(         &
                metrics%ax*(f0(i-1,j) + f0(i+1,j)) +   &
                metrics%ay*(f0(i,j-1) + f0(i,j+1))     &
            )
            df = f1(i,j) - f0(i,j);
            norm_df = norm_df + df*df;
        end do
        end do
        norm_df = sqrt(norm_df)/(nx*ny)
        if (mod(n,np) == 0) write(*,'(i6,1p,e14.6)') n, norm_df
        call swap(f1, f0)
    end do
    call timer%tock()

    ! Result summary
    i  = nx/2 + 1
    j  = ny/2 + 1
    dt = timer%last_elapsed
    write(*,*)
    write(*,'(a,f6.4)')   'Midpoint Value:  ', f0(i,j)
    write(*,'(a,f0.2,a)') 'Elapsed Time:    ', dt, ' sec'
    write(*,'(a,f0.2,a)') 'Througput:       ', (nt/1.0d6)*(nx-2)*(ny-2)/dt, ' MUPS'
    write(*,*)

contains

    function get_arg(pos) result(val)
        integer, intent(in) :: pos
        integer :: val
        character(256) :: buffer
        call get_command_argument(pos, buffer)
        read(buffer,*) val
    end function

    function boundary_value(x,y) result(bv)
        real(8), intent(in) :: x,y
        real(8), parameter  :: a2 = 0.1d0
        real(8) :: bv, r2
        r2 = x*x + y*y
        bv = a2 / (a2 + r2)
    end function

    function simple_update(f0, f1, metrics) result(norm_df)
        real(8), intent(in) :: f0(:,:), f1(:,:)
        type(metrics_t), intent(in) :: metrics
        real(8) :: norm_df
        norm_df = 0.0d0
    end function

    subroutine swap(p1, p2)
        real(8), pointer, intent(inout) :: p1(:,:), p2(:,:)
        real(8), pointer :: ptmp(:,:)
        ptmp => p1
        p1 => p2
        p2 => ptmp
    end subroutine

    subroutine print(f)
        real(8), intent(in) :: f(:,:)
        integer :: i, j
        do i = 1, size(f,1)
            do j = 1, size(f,2)
                write(*,'(f6.4,1x)', advance='no') f(i,j)
            end do
            write(*,*)
        end do
    end subroutine

end program