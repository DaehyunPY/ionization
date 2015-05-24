module inner
<<<<<<< HEAD
<<<<<<< HEAD
    ! use H, E, R, K, A matrix 
    ! not use S matrix 
    use kind_type 
    use global 
    implicit none
    complex(dp),       allocatable, private, protected :: inner_a(:)
    complex(dp), save, allocatable, private, protected :: inner_u(:, :)
contains 
    

! ==================================================
! COEFFICIENT
! ==================================================
! inner coefficient --------------------------------
subroutine inner_coeff(l)
    use math_const,  only: i => math_i 
    use fgsl, only: fgsl_sf_bessel_jsl, fgsl_sf_bessel_ysl
    use hamiltonian, only: coord_E
    integer(i4), intent(in) :: l 
    real   (dp) :: ka, sb_j, sb_y, tmp2
    complex(dp) :: tmp1 
    integer(i4) :: j

    ka = (2.d0*Mass*coord_E(1_i4))**0.5d0*Bound
    sb_j = fgsl_sf_bessel_jsl(l, ka)
    sb_y = fgsl_sf_bessel_ysl(l, ka)
    
    tmp1 = A(l)*(sb_j -K(l)*sb_y)
    do j = 1, N 
        tmp2       = H(j, N)/(2.d0*Mass*(E(j) -coord_E(1_i4)))
        inner_a(j) = tmp1*tmp2/R(l)
    end do     
end subroutine inner_coeff
=======
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
=======
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




>>>>>>> 9f36052... new code




>>>>>>> 9f36052... new code


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

<<<<<<< HEAD
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

<<<<<<< HEAD
! ==================================================
! PROCESS
! ==================================================
! inner achive -------------------------------------
subroutine PROC_inner_achive(l)
    integer(i4), intent(in) :: l 
    complex(qp) :: sum 
    integer(i4) :: i, j

    if(.not. allocated(inner_u)) allocate(inner_u(0:L, 1:N))
    if(.not. allocated(inner_a)) allocate(inner_a(1:N))
    call inner_coeff(l) 
    do i = 1, N 
        sum = 0.d0 
        do j = 1, N 
            sum = sum +inner_a(j)*H(j, i)
        end do 
        inner_u(l, i) = sum 
    end do 
    if(allocated(inner_a)) deallocate(inner_a)
end subroutine PROC_inner_achive
! end inner achive ---------------------------------
! inner plot ---------------------------------------
subroutine PROC_inner_plot 
    use math_const,  only: pi => math_pi, degree => math_degree
    use hamiltonian, only: coord_r, coord_theta
    use fgsl, only: fgsl_sf_legendre_Pl
    integer  (i1), parameter :: file_psi1 = 101, file_psi2 = 102
    character(30), parameter :: form_psi  = '(30ES25.10)'
    real     (dp), parameter :: radian_to_degree = 1.d0/degree 
    real     (dp) :: tmp, unit_theta
    complex  (qp) :: sum 
    integer  (i4) :: i, j, k 

    unit_theta = 1.d0  
    if(op_degree == "Y") unit_theta = radian_to_degree

    open(file_psi1, file = "output/inner_u_0.d")
    sum = 0.d0 
    do i = 1, N
        sum = inner_u(0, i)
        write(file_psi1, form_psi) coord_r(i), dble(abs(sum)**2.d0)
    end do 
    close(file_psi1)

    open(file_psi2, file = "output/inner_psi.d")
    do i = 1, N, N/pr 
        do j = 0, ptheta
            sum = 0.d0 
            do k = 0, L 
                tmp = cos(coord_theta(j))
                sum = sum +inner_u(k, i)/coord_r(i)*fgsl_sf_legendre_Pl(k, tmp)
            end do 
            write(file_psi2, form_psi) coord_r(i), coord_theta(j)*unit_theta, dble(abs(sum)**2.d0)
        end do 
        write(file_psi2, form_psi) 
    end do 
    if(allocated(inner_u)) deallocate(inner_u)
    close(file_psi2)
end subroutine PROC_inner_plot
! end inner plot -----------------------------------
=======
        term_dipole = 0.5d0*amp_grid(m)*r_grid(i)
    end function term_dipole ! :::::::::::::::::::::::::::::::::::::::::::::::::
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! END HAMILTONIAN --------------------------------------------------------------

=======
        term_dipole = 0.5d0*amp_grid(m)*r_grid(i)
    end function term_dipole ! :::::::::::::::::::::::::::::::::::::::::::::::::
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! END HAMILTONIAN --------------------------------------------------------------

>>>>>>> 9f36052... new code









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

<<<<<<< HEAD
>>>>>>> 9f36052... new code
=======
>>>>>>> 9f36052... new code
end module inner






















































