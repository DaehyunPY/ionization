module basis
    use kind_type 
    use global 
    implicit none
    real(dp), allocatable, private, protected :: mat_H(:, :)
contains


! ==================================================
! FUNCTION
! ==================================================
! kinetic term -------------------------------------
function mat_kinet(i, j)
    use hamiltonian, only: Delta_rho 
    integer(i4), intent(in) :: i, j 
    real   (dp) :: mat_kinet, tmp 
    tmp        = -Delta_rho(i, j)/dr_p_drho**2.d0 
    mat_kinet  = tmp/(2.d0*Mass)
end function mat_kinet 
! potential term -----------------------------------
function mat_poten(i)
    use hamiltonian, only: coord_r, Poten_r 
    integer(i4), intent(in) :: i
    real   (dp) :: mat_poten, tmp 
    tmp        = Poten_r(coord_r(i))
    mat_poten  = tmp*Charge
end function mat_poten 
! angular term -------------------------------------
function mat_angular(l, i)  
    use hamiltonian, only: coord_r, angular_r
    integer(i4), intent(in) :: l, i
    real   (dp) :: mat_angular, tmp 
    tmp         = angular_r(l)/coord_r(i)**2.d0
    mat_angular = tmp/(2.d0*Mass)
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
    if(info /= 0) stop "Error #817: subroutine diag"
    if(allocated(work)) deallocate(work)
end subroutine diag










! ==================================================
! PROCESS
! ==================================================
! hamiltonian --------------------------------------
subroutine PROC_H(l) 
    character(30), parameter  :: form_out = '(1A15, 10F9.3)'
    integer  (i4), intent(in) :: l 
    real     (dp) :: sign, tmp 
    integer  (i4) :: i, i1, j

    if(.not. (op_basis == "Y" .or. op_inner == "Y")) then 
        if(.not. allocated(H)) allocate(H(1:N, N:N))
    else if(op_basis == "Y" .or. op_inner == "Y") then 
        if(.not. allocated(H)) allocate(H(1:N, 1:N))
    end if 
    if(.not. allocated(E)) allocate(E(1:N))
    if(.not. allocated(mat_H)) allocate(mat_H(1:N, 1:N))
    mat_H(:, :) = 0.d0
    E(:)        = 0.d0 
    do i = 1, N 
        do j = 1, N
            mat_H(i, j) = mat_H(i, j) +mat_Kinet(i, j)
        enddo
        mat_H(i, i) = mat_H(i, i) +mat_Poten(i) +mat_angular(l, i)
    enddo
    call diag

    H(:, :) = 0.d0
    sign    = 1.d0 
    i1      = 1 
    if(size(H(1, :)) == 1) i1 = N 
    do j = 1, N 
        sign = 1.d0 
        if(mat_H(1, j) < 0.d0) sign = -1.d0 
        do i = i1, N 
            tmp = sign*mat_H(i, j)
            tmp = tmp/(coord_weight(i)*dr_p_drho)**0.5d0 
            H(j, i) = tmp
        end do 
    end do 
    if(allocated(mat_H)) deallocate(mat_H)
    write(file_log, form_out) "Energy: ", (E(i), i = 1, 5) 
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

    write(file_psi, form_gen) (0.d0, j = 0, 20)
    do i = 1, N 
        write(file_psi, form_gen) coord_r(i), (H(j, i), j = 1, 10), (H(j, i), j = N -9, N)
        write(file_ene, form_gen) dble(i), E(i) 
    end do 
    close(file_psi)
    close(file_ene)
end subroutine PROC_basis_plot
! end basis plot -----------------------------------
end module basis
