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
    use fgsl, only: fgsl_sf_bessel_jsl, fgsl_sf_bessel_ysl
    use hamiltonian, only: coord_E
    integer(i4), intent(in) :: l 
    real   (dp), intent(in) :: r 
    real   (dp) :: kr, sb_j, sb_y
    complex(dp) :: outer_u

    kr = (2.0_dp*Mass*coord_E(1_i4))**0.5_dp*r
    sb_j = fgsl_sf_bessel_jsl(l, kr)
    sb_y = fgsl_sf_bessel_ysl(l, kr)

    outer_u = A(l)*(sb_j -K(l)*sb_y)*r 
end function outer_u










! ==================================================
! PROCESS
! ==================================================
! outer plot ---------------------------------------
subroutine PROC_outer_plot 
    use math_const,  only: pi => math_pi, degree => math_degree
    use hamiltonian, only: coord_theta
    use fgsl, only: fgsl_sf_legendre_Pl
    integer  (i1), parameter :: file_psi1 = 101, file_psi2 = 102
    character(30), parameter :: form_psi  = '(30ES25.10)'
    real     (dp), parameter :: radian_to_degree = 1.0_dp/degree 
    real     (dp) :: unit_theta, tmp, r, dr 
    complex  (qp) :: sum 
    integer  (i4) :: i, j, k 

    unit_theta = 1.0_dp 
    if(op_degree == "Y") unit_theta = radian_to_degree
    dr = Bound/dble(N)

    open(file_psi1, file = "output/outer_u_0.d")
    sum = 0.0_dp 
    do i = 1, N, N/pr 
        r = Bound +dr*dble(i)
        sum = outer_u(0_i4, r)
        write(file_psi1, form_psi) r, dble(abs(sum)**2.0_dp)
    end do 
    close(file_psi1)

    open(file_psi2, file = "output/outer_psi.d")
    do i = 1, N, N/pr 
        r = Bound +dr*dble(i)
        do j = 0, ptheta
            sum = 0.0_dp 
            do k = 0, L 
                tmp = cos(coord_theta(j))
                sum = sum +outer_u(k, r)/r*fgsl_sf_legendre_Pl(k, tmp)
            end do 
            write(file_psi2, form_psi) r, coord_theta(j)*unit_theta, dble(abs(sum))**2.0_dp 
        end do 
        write(file_psi2, form_psi) 
    end do 
    close(file_psi2)
end subroutine PROC_outer_plot
! end outer plot -----------------------------------
end module outer
