module inner
    use global
    use hamiltonian, only: &
        freq, r_inner, dr_pdrho, weight_grid, &
        amp_grid, r_grid, delta_grid, poten_r
    implicit none
    real(dp), allocatable :: basis_e(:)
    real(dp), allocatable, private :: H1(:, :)
    real(dp), allocatable, target, private :: basis_p(:)
    real(dp), pointer :: basis_f(:, :, :)
contains










! ==============================================================================
! HAMILTONIAN: term_kinet, term_poten, term_floq, term_dipole
! ==============================================================================
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    function term_kinet(i, j) ! ::::::::::::::::::::::::::::::::::::::::::::::::
        integer, intent(in) :: i, j
        real(dp) :: term_kinet

        term_kinet = -Delta_grid(i, j)/dr_pdrho**2.d0/2.d0
    end function term_kinet ! ::::::::::::::::::::::::::::::::::::::::::::::::::
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    function term_poten(i) ! :::::::::::::::::::::::::::::::::::::::::::::::::::
        integer, intent(in) :: i
        real(dp) :: term_poten

        term_poten = poten_r(r_grid(i))
    end function term_poten ! ::::::::::::::::::::::::::::::::::::::::::::::::::
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    function term_floq(n)
        integer, intent(in) :: n
        real(dp) :: term_floq

        term_floq = -dble(n)*freq
    end function term_floq ! :::::::::::::::::::::::::::::::::::::::::::::::::::
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    function term_dipole(i, m) ! :::::::::::::::::::::::::::::::::::::::::::::::
        integer, intent(in) :: i, m
        real(dp) :: term_dipole

        term_dipole = 0.5d0*amp_grid(m)*r_grid(i)
    end function term_dipole ! :::::::::::::::::::::::::::::::::::::::::::::::::
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! END HAMILTONIAN --------------------------------------------------------------










! ==============================================================================
! SUB-PROCESS: sub_build_H1, sub_descale_H1, sub_build_HF, sub_descale_HF
! ==============================================================================
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    subroutine sub_build_H1 ! ::::::::::::::::::::::::::::::::::::::::::::::::::
        use mylinear, only: diag_sym
        integer :: i, j
        integer, pointer :: N

        N => r_N
        H1(:, :) = 0.d0
        do j = 1, N
            do i = 1, j
                H1(i, j) = term_kinet(i, j)
            end do
            H1(j, j) = H1(j, j) +term_poten(j)
        end do
        basis_f(:, 0, :) = H1(:, :)
        call diag_sym(basis_f(:, 0, :), basis_e(:))
    end subroutine sub_build_H1 ! ::::::::::::::::::::::::::::::::::::::::::::::
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    subroutine sub_descale_H1 ! ::::::::::::::::::::::::::::::::::::::::::::::::
        integer :: j
        integer, pointer :: N

        N => r_N
        do j = 1, N
            basis_f(1:N, 0, j) = basis_f(1:N, 0, j)/(weight_grid(1:N)*dr_pdrho)**0.5d0
        end do
    end subroutine sub_descale_H1 ! ::::::::::::::::::::::::::::::::::::::::::::
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    subroutine sub_build_HF(m) ! :::::::::::::::::::::::::::::::::::::::::::::::
        use mylinear, only: diag_sym
        integer :: j, k
        integer, intent(in) :: m
        integer, pointer :: N, F
        real(dp) :: tmp
        real(dp), pointer :: HF(:, :, :, :), HC(:, :)

        N => r_N
        F => floq_N
        nullify(HF)
        nullify(HC)
        allocate(HF(1:N, -F:F, 1:N, -F:F))
        allocate(HC(1:N*(2*F +1), 1:N*(2*F +1)))
        HF(1:N, -F:F, 1:N, -F:F) => basis_p(1:N*(2*F +1)*N*(2*F +1))
        HC(1:N*(2*F +1), 1:N*(2*F +1)) => basis_p(1:N*(2*F +1)*N*(2*F +1))

        HF(:, :, :, :) = 0.d0
        do k = -F, F
            tmp = term_floq(k)
            HF(:, k, :, k) = H1(:, :)
            do j = 1, N
                HF(j, k, j, k) = HF(j, k, j, k) +tmp
            end do
            if(.not. k == -F) then
                do j = 1, N
                    HF(j, k -1, j, k) = term_dipole(j, m)
                end do
            end if
        end do
        call diag_sym(HC(:, :), basis_e(:))
        nullify(HF)
        nullify(HC)
    end subroutine sub_build_HF ! ::::::::::::::::::::::::::::::::::::::::::::::
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    subroutine sub_descale_HF ! ::::::::::::::::::::::::::::::::::::::::::::::::
        integer :: i, j
        integer, pointer :: N, F

        N => r_N
        F => floq_N
        do j = 1, N*(2*F +1)
            do i = -F, F
                basis_f(1:N, i, j) = basis_f(1:N, i, j)/(weight_grid(1:N)*dr_pdrho)**0.5d0
            end do
        end do
    end subroutine sub_descale_HF ! ::::::::::::::::::::::::::::::::::::::::::::
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! END SUB-PROCESS --------------------------------------------------------------



! ==============================================================================
! PROCESS: process_basis_H1, process_basis_HF
! ==============================================================================
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    subroutine process_basis_H1 ! ::::::::::::::::::::::::::::::::::::::::::::::
        integer :: j
        integer, pointer :: N

        N => r_N
        if(allocated(H1)) deallocate(H1)
        if(allocated(basis_e)) deallocate(basis_e)
        if(allocated(basis_p)) deallocate(basis_p)
        nullify(basis_f)
        allocate(H1(1:N, 1:N))
        allocate(basis_e(1:N))
        allocate(basis_p(1:N*N))
        allocate(basis_f(1:N, 0, 1:N))
        basis_f(1:N, 0:0, 1:N) => basis_p(1:N*N)

        call sub_build_H1
        call sub_descale_H1

!         ! These codes are for test. The results should be...
!         !       -1.2663       8.6033      28.342
!         !     +/-1.4142    +/-1.4142    +/-1.4142
!         print *, (basis_e(j), j = 1, 3)
!         print *, (basis_f(N, 0, j), j = 1, 3)
    end subroutine process_basis_H1 ! ::::::::::::::::::::::::::::::::::::::::::
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    subroutine process_basis_HF(m) ! :::::::::::::::::::::::::::::::::::::::::::
        integer :: i, j, k
        integer, intent(in) :: m
        integer, pointer :: N, F

        N => r_N
        F => floq_N
        if(allocated(basis_e)) deallocate(basis_e)
        if(allocated(basis_p)) deallocate(basis_p)
        nullify(basis_f)
        allocate(basis_e(1:N*(2*F +1)))
        allocate(basis_p(1:N*(2*F +1)*N*(2*F +1)))
        allocate(basis_f(1:N, -F:F, 1:N*(2*F +1)))
        basis_f(1:N, -F:F, 1:N*(2*F +1)) => basis_p(1:N*(2*F +1)*N*(2*F +1))

        call sub_build_HF(m)
        call sub_descale_HF

!         do k = -F, F
!             do j = 1, N
!                 write(11, "(ES15.5, I15, 3ES15.5)") r_grid(j), k, (basis_f(j, k, i), i = 1, 3)
!             end do
!             write(11, "(ES15.5, I15, 3ES15.5)")
!         end do
    end subroutine process_basis_HF ! ::::::::::::::::::::::::::::::::::::::::::
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! END PROCESS ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

end module inner






















































