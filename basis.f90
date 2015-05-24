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
! dipole term -------------------------------------
function mat_dipole(i)  
    use hamiltonian, only: coord_r
    integer(i4), intent(in) :: i
    real   (dp) :: mat_dipole, tmp, alpha = 1.d-3
!     tmp         = Amp*coord_r(i)
    tmp         = Amp*coord_r(i)*(1.d0 -exp(alpha*(coord_r(i) -Bound)))
!     mat_dipole = tmp*Charge 
    mat_dipole = tmp
end function mat_dipole
! floquet term -------------------------------------
! function mat_floquet(n)
!     use hamiltonian, only: coord_r
!     integer(i4), intent(in) :: n
!     real   (dp) :: mat_floquet, tmp 
!     tmp         = dble(n)*Freq
!     mat_floquet = tmp
! end function mat_floquet
! dipole term --------------------------------------
! function mat_dipole(i)
!     use hamiltonian, only: coord_r
!     integer(i4), intent(in) :: i
!     real   (dp) :: mat_dipole, tmp 
!     tmp         = 0.5d0*Amp*coord_r(i)
!     mat_dipole  = tmp*Charge
! end function mat_dipole
! matrix R -----------------------------------------
function func_R(energy)
    real   (dp), intent(in) :: energy
    real   (dp) :: func_R 
    real   (qp) :: sum 
    integer(i4) :: i
    sum = 0.d0 
    do i = 1, N 
        sum = sum +H(i, N)**2.d0/(E(i) -energy)
    end do 
    func_R = sum/(2.d0*Mass*Bound)
end function func_R
! matrix R -----------------------------------------


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
! check matrix -------------------------------------
! subroutine check_norm
!     real(qp) :: sum, total  
!     integer  :: i1, i2, j 
!     total = 0.d0 
!     do j = 1, (2*F +1)*N
!         sum = 0.d0 
!         do i2 = -F, F
!             do i1 = 1, N 
!                 sum = sum +H(j, i2, i1)*H(j, i2, i1)*coord_weight(i1)*dr_p_drho
!             end do 
!         end do 
!         write(*, *) j, dble(sum)
!         total = total +sum 
!     end do 
!     write(*, *) 'total', dble(total/((2*F +1)*N))
! end subroutine check_norm










! ==================================================
! PROCESS
! ==================================================
! hamiltonian --------------------------------------
subroutine PROC_H(l) 
    use hamiltonian, only: coord_r, Poten_r 
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
        mat_H(i, i) = mat_H(i, i) +mat_Poten(i) +mat_angular(l, i) +mat_dipole(i)
!         mat_H(i, i) = mat_H(i, i) +mat_Poten(i) +mat_angular(l, i) 
    enddo
!     call check_mat ! for test 
    call diag
!     call check_mat ! for test 

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
!     call check_norm ! for test 
    if(allocated(mat_H)) deallocate(mat_H)
    write(file_log, form_out) "Energy: ", (E(i), i = 1, 5) 

    do j = 1, N 
        write(30, *) coord_r(j), mat_poten(j) +mat_dipole(j)
    end do 
    write(30, *)
    write(30, *)
    print *, "hello"
end subroutine PROC_H
! end hamiltonian ----------------------------------
! basis plot ---------------------------------------
subroutine PROC_basis_plot(l)
    use hamiltonian, only: coord_r
    integer  (i1), parameter  :: file_psi = 101, file_ene = 102
    character(30), parameter  :: form_gen = '(25ES25.10)'
    integer  (i4), intent(in) :: l 
    integer  (i4) :: i, j
    character (3) :: ch 

    write(ch, '(I3.3)') l 
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
! bound states -------------------------------------
subroutine PROC_bound(l)
    use fgsl, only: fgsl_sf_bessel_ksl_scaled
    integer  (i1), parameter  :: file_bound = 101
    character(30), parameter  :: form_bound = '(2ES25.10)', form_out = '(1A15, 1I15, 2ES15.3)'
    integer  (i4), intent(in) :: l 
    real     (dp) :: dE, E0, Ei, tmp1, tmp2, d1, d2 
    real     (dp) :: ka, sb_ka, diff_ka 
    integer  (i4) :: i, j, m, num 
    character (3) :: ch 

    write(ch, '(I3.3)') l 
    open(file_bound, file = "output/bound_"//ch//".d")
    m   = 1000
    num = 0 
    do j = 1, N 
        tmp1 = E(j)
        tmp2 = E(j +1)
        if(.not. tmp1 < 0.d0) exit 
        if(.not. tmp2 < 0.d0) tmp2 = 0.d0 
        E0  = tmp1
        dE  = (tmp2 -tmp1)/dble(m)
        d1  = 0.d0 
        do i = 1, m -1 
            Ei = E0 +dble(i)*dE 
            ka = (2.d0*Mass*(-Ei))**0.5d0*Bound
            sb_ka = fgsl_sf_bessel_ksl_scaled(l, ka)/exp(ka)*ka 
            if(.not. l == 0) then 
                tmp1 = dble(l)*fgsl_sf_bessel_ksl_scaled(l, ka)/exp(ka)
                tmp2 = dble(l +1)*fgsl_sf_bessel_ksl_scaled(l +1_i4, ka)/exp(ka)
            else if(l == 0) then 
                tmp1 = 0.d0 
                tmp2 = dble(l +1)*fgsl_sf_bessel_ksl_scaled(l +1_i4, ka)/exp(ka)
            end if 
            diff_ka = -ka**2.d0/dble(2*l +1)*(tmp1 +tmp2)
            d2 = sb_ka -func_R(Ei)*(sb_ka +diff_ka)
            write(file_bound, form_bound) Ei, d2
            if(d1*d2 < 0.d0) then 
                num  = num +1 
                tmp1 = -dE*abs(d2)/(abs(d1) +abs(d2))
                tmp2 = Ei +tmp1 
                write(file_log, form_out) "Bound State: ", num, tmp2
            end if 
            d1 = d2 
        end do 
    end do 
    close(file_bound)
end subroutine PROC_bound
! end bound states ---------------------------------
end module basis
