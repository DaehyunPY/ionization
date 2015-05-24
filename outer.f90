module outer
    ! use K, A
    ! not use H, E, R, S
    use kind_type
    use global
    implicit none
contains


! ==================================================
! COEFFICIENT
! ==================================================
! outer coefficient --------------------------------
function outer_u(l, r)
    use gsl_bessel,  only: gsl_sf_bessel_jsl, gsl_sf_bessel_ysl
    use hamiltonian, only: coord_E
    integer(i4), intent(in) :: l
    real   (dp), intent(in) :: r
    real   (dp) :: kr, sb_j, sb_y
    complex(dp) :: outer_u

    kr = (2.d0*Mass*coord_E(1_i4))**0.5d0*r
    sb_j = gsl_sf_bessel_jsl(l, kr)
    sb_y = gsl_sf_bessel_ysl(l, kr)

    outer_u = A(l)*(sb_j -K(l)*sb_y)*r
end function outer_u










! ==================================================
! PROCESS
! ==================================================
! outer plot ---------------------------------------
subroutine PROC_outer_plot
    use math_const,   only: pi => math_pi, degree => math_degree
    use hamiltonian,  only: coord_theta
    use gsl_legendre, only: gsl_sf_legendre_Pl
    integer  (i1), parameter :: file_psi1 = 101, file_psi2 = 102
    character(30), parameter :: form_psi  = '(30ES25.10)'
    real     (dp), parameter :: radian_to_degree = 1.d0/degree
    real     (dp) :: unit_theta, tmp, r, dr
    complex  (qp) :: sum
    integer  (i4) :: i, j, k

    unit_theta = 1.d0
    if(op_degree == 1) unit_theta = radian_to_degree
    dr = Ra/dble(N)

    open(file_psi1, file = "output/outer_u_0.d")
    sum = 0.d0
    do i = 1, N, N/pr
        r = Ra +dr*dble(i)
        sum = outer_u(0_i4, r)
        write(file_psi1, form_psi) r, dble(abs(sum)**2.d0)
    end do
    close(file_psi1)

    open(file_psi2, file = "output/outer_psi.d")
    do i = 1, N, N/pr
        r = Ra +dr*dble(i)
        do j = 0, ptheta
            sum = 0.d0
            do k = 0, L
                tmp = cos(coord_theta(j))
                sum = sum +outer_u(k, r)/r*gsl_sf_legendre_Pl(k, tmp)
            end do
            write(file_psi2, form_psi) r, coord_theta(j)*unit_theta, dble(abs(sum))**2.d0
        end do
        write(file_psi2, form_psi)
    end do
    close(file_psi2)
end subroutine PROC_outer_plot
! end outer plot -----------------------------------
end module outer
