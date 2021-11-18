module utils
    implicit none

    type :: metrics_t
        real(8) :: dx, dy
        real(8) :: ax, ay, dxy
        real(8) :: inv_dxy
    end type
    interface metrics_t
        module procedure :: metrics_init
    end interface

    type :: timer_t
        integer :: count = 0
        integer :: last_tick = 0.0d0
        integer :: last_tock = 0.0d0
        double precision :: clock_rate
        double precision :: last_elapsed = 0.0d0
        double precision :: average_elapsed = 0.0d0
    contains
        procedure :: tick => timer_tick
        procedure :: tock => timer_tock
    end type

contains

    function metrics_init(nx,ny) result(m)
        integer, intent(in) :: nx, ny
        type(metrics_t) :: m
        m%dx  = 1.0d0 / (nx-1);
        m%dy  = 1.0d0 / (ny-1);
        m%ax  = 1.0d0 / (m%dx*m%dx);
        m%ay  = 1.0d0 / (m%dy*m%dy);
        m%dxy = 0.5d0 * (m%ax + m%ay);
        m%inv_dxy = 1.0d0/m%dxy;
    end function

    subroutine timer_tick(self)
        class(timer_t), intent(inout) :: self
        integer :: current_tick
        call system_clock(current_tick)
        if (self%last_tick > self%last_tock) then
            write(*,*) "Warning: No tock() since last tick(). Reseting timer."
        endif
        self%last_tick = current_tick
        return
    end subroutine

    subroutine timer_tock(self)
        class(timer_t), intent(inout) :: self
        double precision :: k
        integer :: current_tick
        call system_clock(current_tick, self%clock_rate)
        if (self%last_tock > self%last_tick) then
            write(*,*) "Warning: No tick() since last tock(). Ignoring tock()."
            return
        end if

        ! Save latest data
        self%count        = self%count+1
        self%last_tock    = current_tick
        self%last_elapsed = (self%last_tock - self%last_tick)/self%clock_rate

        ! Compute the running mean
        k = 1.0d0/self%count
        self%average_elapsed = (1.0d0-k)*self%average_elapsed + k*self%last_elapsed

        return
    end subroutine

end module