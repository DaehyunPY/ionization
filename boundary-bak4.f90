module boundary
    use kind_type
    use global 
    implicit none
contains 


subroutine mat_R(l)
    character(30), parameter  :: form_out = '(1A15, 1ES15.3)'
    integer  (i4), intent(in) :: l 
    real     (qp) :: sum 
    integer  (i4) :: i 

    sum = 0.d0 
    do i = 1, (2*F +1)*N 
        sum = sum +H(i, N)**2.d0/(2.d0*Mass*(E(i) -Scatt)) 
    end do 
    R(l) = sum/ra 
    write(file_log, form_out) "R: ", R(l)
end subroutine mat_R


subroutine mat_K(l)
    use gsl_special, only: gsl_sf_bessel_jsl, gsl_sf_bessel_ysl
    character(30), parameter  :: form_out = '(1A15, 1ES15.3)'
    integer  (i4), intent(in) :: l 
    real     (dp) :: ka, sb_j, sb_y, diff_j, diff_y 
    real     (dp) :: agamma, tmp1, tmp2 

    ka = (2.d0*Mass*Scatt)**0.5d0*ra
    sb_j = gsl_sf_bessel_jsl(l, ka)
    sb_y = gsl_sf_bessel_ysl(l, ka)
    if(l /= 0_i4) then 
        diff_j = gsl_sf_bessel_jsl(l -1_i4, ka) -dble(l +1_i4)/ka*gsl_sf_bessel_jsl(l, ka)
        diff_y = gsl_sf_bessel_ysl(l -1_i4, ka) -dble(l +1_i4)/ka*gsl_sf_bessel_ysl(l, ka)
    else if(l == 0_i4) then 
        diff_j = -gsl_sf_bessel_jsl(0_i4 +1_i4, ka) -dble(0_i4 +1_i4)/ka*gsl_sf_bessel_jsl(0_i4, ka)
        diff_y = -gsl_sf_bessel_ysl(0_i4 +1_i4, ka) -dble(0_i4 +1_i4)/ka*gsl_sf_bessel_ysl(0_i4, ka)
    end if 

    agamma = 1.d0/R(l) -1.d0
    tmp1   = ka*diff_j -agamma*sb_j 
    tmp2   = ka*diff_y -agamma*sb_y 
    K(l)   = tmp1/tmp2
    write(file_log, form_out) "K: ", K(l)
end subroutine mat_K


subroutine mat_S(l)
    use math_const, only: i => math_i 
    character(60), parameter  :: form_out = '(1A15, 1ES15.3, 1ES15.3, "i")'
    integer  (i4), intent(in) :: l 

    S(l) = (1.d0 +i*K(l))/(1.d0 -i*K(l))
    write(file_log, form_out) "S: ", S(l)
end subroutine mat_S


subroutine mat_A(l)
    use math_const, only: i => math_i 
    character(60), parameter  :: form_out = '(1A15, 1ES15.3, 1ES15.3, "i")'
    integer  (i4), intent(in) :: l 
    complex  (dp) :: sign 
    real     (dp) :: tmp1, tmp2 

    if(mod(l, 4_i4) == 0_i4) then 
        sign = 1.d0 
    else if(mod(l, 4_i4) == 1_i4) then 
        sign = i 
    else if(mod(l, 4_i4) == 2_i4) then 
        sign = -1.d0 
    else if(mod(l, 4_i4) == 3_i4) then 
        sign = -i 
    end if 
    tmp1 = (1.d0 +real(S(l)))/2.d0 
    tmp2 = aimag(S(l))/2.d0 
    A(l) = sign*(2.d0*dble(l) +1.d0)*exp(tmp1 +i*tmp2)
    write(file_log, form_out) "A: ", A(l)
end subroutine mat_A












! ==================================================
! PROCESS
! ==================================================
! boundary matrix ----------------------------------
subroutine PROC_boundary_mat(l) ! It must be called after PROC_input, PROC_H 
    integer(i4), intent(in) :: l 

    call mat_R(l)
    call mat_K(l)
    call mat_S(l)  
    call mat_A(l)
end subroutine PROC_boundary_mat
! end boundary matrix ------------------------------
end module boundary
