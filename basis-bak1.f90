module basis
    use kind_type 
    use global 
    implicit none
    real(dp), pointer, private, protected :: mat_H(:, :)
contains


! ==================================================
! FUNCTION
! ==================================================
! kinetic term -------------------------------------
function mat_Kinet(i, j)
    use hamiltonian, only: Delta_rho 
    integer(i4), intent(in) :: i, j 
    real   (dp) :: mat_Kinet, tmp 
    tmp       = -Delta_rho(i, j)/dr_p_drho**2_dp 
    mat_Kinet = tmp/(2_dp*Mass)
end function mat_Kinet 
! potential term -----------------------------------
function mat_Poten(i)
    use hamiltonian, only: coord_r, Poten_r 
    integer(i4), intent(in) :: i
    real   (dp) :: mat_Poten, tmp 
    tmp       = Poten_r(coord_r(i))
    mat_Poten = tmp*Charge
end function mat_Poten 
! angular term -------------------------------------
function mat_angular(l, i)
    use hamiltonian, only: coord_r, angular_r
    integer(i4), intent(in) :: l, i
    real   (dp) :: mat_angular, tmp 
    tmp         = angular_r(l)/coord_r(i)**2_dp
    mat_angular = tmp/(2_dp*Mass)
end function mat_angular


! ==================================================
! CALCULATION 
! ==================================================
! dsyev diagonalization ----------------------------
subroutine diag
    character(1), parameter   :: jobz = 'Vectors', uplo = 'Upper'
    real    (dp), allocatable :: work(:)
    integer (i8), save        :: lwork = -1, n, lda 
    integer (i8) :: info 
    
    if(lwork < 0) then 
        n    = size(mat_H(:, 1))
        lda  = size(mat_H(1, :))
        info = 0 
        if(.not. allocated(work)) allocate(work(1))
        call DSYEV(jobz, uplo, n, mat_H, lda, E, work, lwork, info)
        lwork = int(work(1))
        if(allocated(work)) deallocate(work)
    end if 

    info = 0 
    if(.not. allocated(work)) allocate(work(1:lwork))
    call DSYEV(jobz, uplo, n, mat_H, lda, E, work, lwork, info)
    if(info /= 0) stop "Error #2: subroutine diag"
    if(allocated(work)) deallocate(work)
end subroutine diag


! ==================================================
! TEST
! ==================================================
! check matrix -------------------------------------
! subroutine check_mat
!     integer  (i1), parameter :: file_check = 101
!     character(30), parameter :: & 
!         form_num   = '(20X, 5000(I15, 5X))', &
!         form_check = '(20X, 5000(ES15.5, 5X))', &
!         form_ch    = '(1I10)', &
!         form_ch1   = '(1I15, 5X, ', &
!         form_ch2   = '20X, ', &
!         form_ch3   = '(ES15.5, 5X), 20X, ', &
!         form_ch4   = '(ES15.5, 5X))'
!     character(10) :: ch1, ch2 
!     integer  (i4) :: n, i, j 

!     n = size(mat_H(:, 1))
!     open(file_check, file = "output/check_mat.d")
!     write(file_check, form_num) (j, j = 1, n)
!     do i = 1, n 
!         if(i == 1) then 
!             write(ch1, form_ch) 0 
!             write(ch2, form_ch) n -1
!             write(file_check, form_ch1//form_ch2//ch2//form_ch4) i, (mat_H(i, j), j = 2, n)
!             ch1 = ""
!             ch2 = ""
!         else if(i == n) then 
!             write(ch1, form_ch) n -1 
!             write(ch2, form_ch) 0 
!             write(file_check, form_ch1//ch1//form_ch4) i, (mat_H(i, j), j = 1, n -1) 
!             ch1 = ""
!             ch2 = ""
!         else 
!             write(ch1, form_ch) i -1 
!             write(ch2, form_ch) n -i 
!             write(file_check, form_ch1//ch1//form_ch3//ch2//form_ch4) i, (mat_H(i, j), j = 1, i -1), (mat_H(i, j), j = i +1, n)
!             ch1 = ""
!             ch2 = ""
!         end if 
!     end do 
!     write(file_check, *)
!     write(file_check, form_num) (j, j = 1, n)
!     write(file_check, form_check) (mat_H(j, j), j = 1, n)
!     write(file_check, *)
!     close(file_check)
! end subroutine check_mat










! ==================================================
! PROCESS
! ==================================================
! hamiltonian --------------------------------------
subroutine PROC_H(l) 
    character(30), parameter  :: form_out = '(1A15, 10F10.3)'
    integer  (i4), intent(in) :: l 
    real     (dp), pointer    :: p1_H(:, :, :, :), p2_H(:, :, :), p_H(:)
    real     (dp) :: sign, tmp 
    integer  (i4) :: i1, i2, j1, j2 

    if(associated(mat_H)) nullify(mat_H)
    if(associated(p1_H))  nullify(p1_H)
    if(associated(p2_H))  nullify(p2_H)
    if(associated(p_H))   nullify(p_H)
    allocate(mat_H(1:N*(2*F +1), 1:N*(2*F +1)))
    allocate(p1_H (1:N, -F:F,    1:N, -F:F))
    allocate(p2_H (1:N, -F:F,    1:N*(2*F +1)))
    allocate(p_H  (1:N*(2*F +1)   *N*(2*F +1)))
    p1_H (1:N, -F:F,    1:N, -F:F)    => p_H(1:N*(2*F +1)*N*(2*F +1))
    p2_H (1:N, -F:F,    1:N*(2*F +1)) => p_H(1:N*(2*F +1)*N*(2*F +1))
    mat_H(1:N*(2*F +1), 1:N*(2*F +1)) => p_H(1:N*(2*F +1)*N*(2*F +1))
    p_H(:) = 0_dp
    E(:)   = 0_dp 
    do j2 = -F, F 
    do j1 =  1, N 
        do i2 = max(j2 -1_i4, -F), min(j2 +1_i4, F)
        do i1 = 1, N
            tmp = 0_dp 
            if(i2 == j2) then 
                             tmp = tmp +mat_Kinet(i1, j1)
                if(i1 == j1) tmp = tmp +mat_Poten(i1) +mat_angular(l, i1)
                if(i1 == j1) tmp = tmp +dble(i2)*Freq
            else if(i2 == j2 -1 .or. i2 == j2 +1) then 
                if(i1 == j1) tmp = tmp +0.5_dp*Amp
            end if 
            p1_H(i1, i2, j1, j2) = tmp 
        end do 
        end do
    end do 
    end do
!     call check_mat
    call diag

    H(:, :) = 0_dp
    sign    = 1_dp 
    i1      = 1 
    j1      = apho
!     if(size(H(1, 0, :)) == 1) i1 = N 
    if(size(H(1, :)) == 1) i1 = N 
    do i2 = 1, (2*F +1)*N 
        sign = 1_dp 
        if(mat_H(1, i2) < 0_dp) sign = -1_dp 
        do j2 = i1, N 
!         do j1 = -F, F 
            tmp = sign*p2_H(j2, j1, i2)
            tmp = tmp/(coord_weight(j2)*dr_p_drho)**0.5_dp
!             tmp = tmp*dble(2*F +1)**0.5_dp 
!             H(i2, j1, j2) = tmp
            H(i2, j2) = tmp
!         end do 
        end do 
    end do 
    p_H(:) = 0_dp
    if(associated(mat_H)) nullify(mat_H)
    if(associated(p1_H))  nullify(p1_H)
    if(associated(p2_H))  nullify(p2_H)
    if(associated(p_H))   nullify(p_H)
    write(file_log, form_out) "Energy: ", (E(i1), i1 = 1, 5) 
end subroutine PROC_H
! end hamiltonian ----------------------------------
! basis plot ---------------------------------------
subroutine PROC_basis_plot(num)
    use hamiltonian, only: coord_r
    integer  (i1), parameter  :: file_psi = 101, file_ene = 102
    character(30), parameter  :: form_gen = '(1000ES25.10)'
    integer  (i4), intent(in) :: num 
    integer  (i4) :: i, j
    character (3) :: ch 

    write(ch, '(I3.3)') num 
    open(file_psi, file = "output/basis_u_"//ch//".d")
    open(file_ene, file = "output/basis_energy_"//ch//".d")

    write(file_psi, form_gen) 0_dp, (0_dp, j = 1, 20)
    do i = 1, N 
        write(file_psi, form_gen) coord_r(i), (H(j, i), j = 1, 10), (H(j, i), j = N -9, N)
    end do 
    do i = 1, (2*F +1)*N 
        write(file_ene, form_gen) dble(i), E(i) 
    end do 
    close(file_psi)
    close(file_ene)
end subroutine PROC_basis_plot
! end basis plot -----------------------------------
end module basis
