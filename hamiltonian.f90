module hamiltonian
    use global
    implicit none
    ! note: hamiltonian > laser field: freq, amp_max, amp_min
    real(dp) :: freq, amp_max, amp_min
    real(dp), private :: damp
    ! note: hamiltonian > potential: poten_type
    integer :: poten_type
    ! note: hamiltonian > coodinate r: r_inner, dr_pdrho, weight_grid
    real(dp) :: r_inner, dr_pdrho
    real(dp), allocatable :: weight_grid(:)
    real(dp), allocatable, private :: rho_grid(:), dshape_grid(:, :)

contains










! ==============================================================================
! COORDINATES & OPERATORS: amp_grid, r_grid, delta_grid, poten_r
! ==============================================================================
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    function amp_grid(i) ! :::::::::::::::::::::::::::::::::::::::::::::::::::::
        integer, intent(in) :: i
        real(dp) :: amp_grid

        amp_grid = amp_min +dble(i)*damp
    end function amp_grid ! ::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    function r_grid(i) ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::
        integer, intent(in) :: i
        real(dp) :: r_grid

        r_grid = dr_pdrho*(rho_grid(i) -rho_grid(0))
    end function r_grid ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    function delta_grid(i, j) ! ::::::::::::::::::::::::::::::::::::::::::::::::
        integer :: k
        integer, intent(in) :: i, j
        integer, pointer :: N
        real(dp) :: delta_grid
        real(qp) :: sum

        N => r_N
        sum = 0.d0
        do k = 0, N
            sum = sum -dshape_grid(i, k)*dshape_grid(j, k)
        end do
        delta_grid = sum
    end function delta_grid ! ::::::::::::::::::::::::::::::::::::::::::::::::::
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    function poten_r(r) ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::
        real(dp), intent(in) :: r
        real(dp) :: poten_r

        poten_r = -2.5d0
    end function poten_r ! :::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! END COORDINATES & OPERATORS --------------------------------------------------










! ==============================================================================
! PROCESS: process_coord
! ==============================================================================
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    subroutine process_coord ! :::::::::::::::::::::::::::::::::::::::::::::::::
        ! It builds a coordinate system of amplitude and r. we use
        ! Lobatto Quadrature for coordinate r. If you want to check r, visit
        ! this page: http://mathworld.wolfram.com/LobattoQuadrature.html.
        use fgsl, only: fgsl_sf_legendre_Pl
        use mylinear, only:diag_sym_band
        integer :: i, j, k
        integer, pointer :: N
        real(dp) :: tmp
        real(dp), allocatable :: X(:, :)
        real(qp) :: sum

        N => r_N
        if(allocated(rho_grid))    deallocate(rho_grid)
        if(allocated(weight_grid)) deallocate(weight_grid)
        if(allocated(dshape_grid)) deallocate(dshape_grid)
        if(allocated(X))           deallocate(X)
        allocate(rho_grid(0:N))
        allocate(weight_grid(0:N))
        allocate(dshape_grid(0:N, 0:N))
        allocate(X(1:2, 1:N -1))

        damp = (amp_max -amp_min)/amp_N

        X = 0.d0
        do i = 1, N -2
            tmp = (dble(i)*dble(i +2))/(dble(2*i +1)*dble(2*i +3))
            tmp = tmp**0.5d0
            X(1, i +1) = tmp
        end do
        call diag_sym_band(X)

        tmp = 2.d0/(dble(N +1)*dble(N))
        rho_grid   (0) = -1.d0
        weight_grid(0) = tmp
        do i = 1, N -1
            rho_grid   (i) = X(1, i)
            weight_grid(i) = tmp/(fgsl_sf_legendre_Pl(N, rho_grid(i)))**2.d0
        end do
        rho_grid   (N) = 1.d0
        weight_grid(N) = tmp
        dr_pdrho = r_inner/(rho_grid(N) -rho_grid(0))
        if(allocated(X)) deallocate(X)

        do i = 0, N
            do j = 0, N
                if(i == j) then
                    sum = 0.d0
                    do k = 0, N
                        if(.not. k == i) then
                            sum = sum +1.d0/(rho_grid(i) -rho_grid(k))
                        end if
                    end do
                else if(i /= j) then
                    sum = 1.d0/(rho_grid(i) -rho_grid(j))
                    do k = 0, N
                        if(.not. (k == i .or. k == j)) then
                            sum = sum*(rho_grid(j) -rho_grid(k))/(rho_grid(i) -rho_grid(k))
                        end if
                    end do
                end if
                sum = sum*(weight_grid(j)/weight_grid(i))**0.5d0
                dshape_grid(i, j) = sum
            end do
        end do

!         ! These codes are for test. If r_N is 5, the results should be...
!         !     +/-1.000000    +/-0.765055    +/-1.000000
!         !     +/-0.666667    +/-0.378475    +/-0.554858
!         print *, rho_grid(:)
!         print *, weight_grid(:)
    end subroutine process_coord ! :::::::::::::::::::::::::::::::::::::::::::::
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! END PROCESS ------------------------------------------------------------------

end module hamiltonian




















































