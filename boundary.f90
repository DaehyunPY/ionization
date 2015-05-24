module boundary
    use kind_type
    use global 
    implicit none
    integer(i1), save, private, protected :: op_mat(1:6) = -1 
contains 


! ==================================================
! MATRIX 
! ==================================================
! matrix R -----------------------------------------
subroutine mat_R(m, l)
    use hamiltonian, only: coord_E
    character(30), parameter  :: form_out = '(1A15, 1ES15.3)'
    integer  (i4), intent(in) :: m, l 
    real     (qp) :: sum 
    integer  (i4) :: i
    if(.not. allocated(R)) allocate(R(l:l))
    sum = 0.d0 
    do i = 1, N 
        sum = sum +H(i, N)**2.d0/(E(i) -coord_E(m))
    end do 
    R(l) = sum/(2.d0*Mass*Bound)
    write(file_log, form_out) "R: ", R(l)
end subroutine mat_R
! matrix K -----------------------------------------
subroutine mat_K(m, l)
    use hamiltonian, only: coord_E
    use fgsl, only: fgsl_sf_bessel_jsl, fgsl_sf_bessel_ysl
    character(30), parameter  :: form_out = '(1A15, 1ES15.3)'
    integer  (i4), intent(in) :: m, l 
    real     (dp) :: ka, sb_j, sb_y, diff_j, diff_y 
    real     (dp) :: agamma, tmp1, tmp2 
    if(.not. allocated(K)) allocate(K(l:l))
    ka = (2.d0*Mass*coord_E(m))**0.5d0*Bound
    sb_j = fgsl_sf_bessel_jsl(l, ka)
    sb_y = fgsl_sf_bessel_ysl(l, ka)
    if(l /= 0) then 
        diff_j = fgsl_sf_bessel_jsl(l -1_i4, ka) -dble(l +1)/ka*fgsl_sf_bessel_jsl(l, ka)
        diff_y = fgsl_sf_bessel_ysl(l -1_i4, ka) -dble(l +1)/ka*fgsl_sf_bessel_ysl(l, ka)
    else if(l == 0) then 
        diff_j = -fgsl_sf_bessel_jsl(1_i4, ka) -dble(1)/ka*fgsl_sf_bessel_jsl(0_i4, ka)
        diff_y = -fgsl_sf_bessel_ysl(1_i4, ka) -dble(1)/ka*fgsl_sf_bessel_ysl(0_i4, ka)
    end if 
    agamma = 1.d0/R(l) -1.d0
    tmp1   = ka*diff_j -agamma*sb_j 
    tmp2   = ka*diff_y -agamma*sb_y 
    K(l)   = tmp1/tmp2
    write(file_log, form_out) "K: ", K(l)
end subroutine mat_K
! matrix S -----------------------------------------
subroutine mat_S(m, l)
    use math_const, only: i => math_i 
    character(60), parameter  :: form_out = '(1A15, 1ES15.3, 1ES15.3, "i")'
    integer  (i4), intent(in) :: m, l 
    if(.not. allocated(S)) allocate(S(m:m, l:l))
    S(m, l) = (1.d0 +i*K(l))/(1.d0 -i*K(l))
    write(file_log, form_out) "S: ", S(m, l)
end subroutine mat_S
! matrix A -----------------------------------------
subroutine mat_A(m, l)
    use math_const, only: i => math_i 
    character(60), parameter  :: form_out = '(1A15, 1ES15.3, 1ES15.3, "i")'
    integer  (i4), intent(in) :: m, l 
    complex  (dp) :: sign 
    real     (dp) :: tmp1, tmp2 
    if(.not. allocated(A)) allocate(A(l:l))
    if(mod(l, 4_i4) == 0) then 
        sign = 1.d0 
    else if(mod(l, 4_i4) == 1) then 
        sign = i 
    else if(mod(l, 4_i4) == 2) then 
        sign = -1.d0 
    else if(mod(l, 4_i4) == 3) then 
        sign = -i 
    end if 
    tmp1 = (1.d0 +real(S(m, l)))/2.d0 
    tmp2 = aimag(S(m, l))/2.d0 
    A(l) = sign*(2.d0*dble(l) +1.d0)*exp(tmp1 +i*tmp2)
    write(file_log, form_out) "A: ", A(l)
end subroutine mat_A










! ==================================================
! PROCESS
! ==================================================
! boundary matrix ----------------------------------
subroutine PROC_boundary_mat(input_m, input_l)
    integer(i4), intent(in) :: input_m, input_l 

    if(op_mat(1) < 0) then 
        ! H(1), E(2), R(3), K(4), S(5), A(6) matrix 
        op_mat(:) = 0 
        if(op_inner == "Y") then 
            ! use H, E, R, K, A matrix 
            ! not use S matrix 
            op_mat(1) = 1 ! H
            op_mat(2) = 1 ! E
            op_mat(3) = 1 ! R
            op_mat(4) = 1 ! K
            op_mat(6) = 1 ! A
        end if 
        if(op_outer == "Y") then 
            ! use K, A 
            ! not use H, E, R, S
            op_mat(4) = 1 ! K
            op_mat(6) = 1 ! A 
        end if 
        if(op_dcs == "Y" .or. op_tcs == "Y" .or. op_ps == "Y") then 
            ! use only S matrix 
            op_mat(5) = 1 ! S 
        end if 
    end if 
    if(op_mat(3) == 1 .and. (.not. allocated(R))) allocate(R(0:L))
    if(op_mat(4) == 1 .and. (.not. allocated(K))) allocate(K(0:L))
    if(op_mat(5) == 1 .and. (.not. allocated(S))) allocate(S(1:M, 0:L))
    if(op_mat(6) == 1 .and. (.not. allocated(A))) allocate(A(0:L))

    call mat_R(input_m, input_l) ! use only H, E matrix 
    call mat_K(input_m, input_l) ! use onyl R matrix 
    call mat_S(input_m, input_l) ! use only K matrix   
    if(op_mat(6) == 1) call mat_A(input_m, input_l) ! use only S matrix 

    if(op_mat(3) == 0 .and. allocated(R)) deallocate(R)
    if(op_mat(4) == 0 .and. allocated(K)) deallocate(K)
    if(op_mat(5) == 0 .and. allocated(S)) deallocate(S)
    if(op_mat(6) == 0 .and. allocated(A)) deallocate(A)
end subroutine PROC_boundary_mat
! end boundary matrix ------------------------------
! matrix out ---------------------------------------
subroutine PROC_mat_out 
    if(op_mat(1) == 0 .and. allocated(H)) deallocate(H)
    if(op_mat(2) == 0 .and. allocated(E)) deallocate(E)
end subroutine PROC_mat_out 
! end matrix out -----------------------------------
end module boundary
