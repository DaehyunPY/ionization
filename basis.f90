module basis
    use kind_const
    use global
    implicit none
    real(dp), save, allocatable, protected :: H1(:, :), E1(:), EF(:)
    real(dp), save, pointer, protected :: HF(:, :, :)
contains


! ==================================================
! HARMILTONIAN TERMS
! ==================================================
! kinetic term -------------------------------------
function term_kinet(i, j)
    use hamiltonian, only: Delta_grid, dr_pdrho
    integer(i4), intent(in) :: i, j
    real(dp) :: term_kinet
    term_kinet = -Delta_grid(i, j)/dr_pdrho**2.d0/2.d0
end function term_kinet
! potential term -----------------------------------
function term_poten(i)
    use hamiltonian, only: r_grid, poten_r
    integer(i4), intent(in) :: i
    real(dp) :: term_poten
    term_poten = poten_r(r_grid(i))
end function term_poten
! angular term -------------------------------------
function term_angular(i, l)
    use hamiltonian, only: r_grid
    integer(i4), intent(in) :: i, l
    real(dp) :: term_angular
    term_angular = dble(l)*(dble(l) +1.d0)/r_grid(i)**2.d0/2.d0
end function term_angular
! floquet term -------------------------------------
function term_floquet(n)
    use hamiltonian, only: r_grid
    integer(i4), intent(in) :: n
    real(dp) :: term_floquet
    term_floquet = -dble(n)*Freq
end function term_floquet
! dipole term --------------------------------------
function term_dipole(i, l, m)
    use hamiltonian, only: Amp_grid, r_grid, dipole_r
    integer(i4), intent(in) :: i, l, m
    real(dp) :: term_dipole
    term_dipole = 0.5d0*Amp_grid(m)*dipole_r(r_grid(i), l)
end function term_dipole
! end funciton -------------------------------------










! ==================================================
! SUB-CALCULATE
! ==================================================
! single hamiltonian -------------------------------
subroutine SUB_single(l)
    use mylinear, only: diag_sym
    integer(i4), intent(in) :: l
    integer(i4) :: i, j
    H1(:, :) = 0.d0
    E1(:) = 0.d0

    do j = 1, N
        do i = 1, j
            H1(i, j) = term_kinet(i, j)
        end do
        H1(j, j) = H1(j, j) +term_poten(j) +term_angular(j, l)
    end do
    call diag_sym(H1(:, :), E1(:))
end subroutine SUB_single
! floquet hamiltonian ------------------------------
subroutine SUB_floquet(l, m)
    use mylinear, only: diag_sym
    integer(i4), intent(in) :: l, m
    real(dp), pointer :: HF_p0(:), HF_p1(:, :, :, :), HF_p2(:, :)
    real(dp) :: tmp
    integer(i4) :: i, j, k
    nullify(HF_p0)
    nullify(HF_p1)
    nullify(HF_p2)
    allocate(HF_p0(1:N*(2*F +1)*N*(2*F +1)))
    allocate(HF_p1(1:N, -F:F, 1:N, -F:F))
    allocate(HF_p2(1:N*(2*F +1), 1:N*(2*F +1)))
    HF   (1:N, -F:F, 1:N*(2*F +1))    => HF_p0(1:N*(2*F +1)*N*(2*F +1))
    HF_p1(1:N, -F:F, 1:N, -F:F)       => HF_p0(1:N*(2*F +1)*N*(2*F +1))
    HF_p2(1:N*(2*F +1), 1:N*(2*F +1)) => HF_p0(1:N*(2*F +1)*N*(2*F +1))
    HF_p1(:, :, :, :) = 0.d0

    ! note: x1, floquet block, x2, floquet block
    do k = -F, F
        tmp = term_floquet(k)
        do j = 1, N
            if(.not. k == -F) then
                HF_p1(j, k -1, j, k) = term_dipole(j, l, m)
            end if
            do i = 1, j
                HF_p1(i, k, j, k) = term_kinet(i, j)
            end do
            HF_p1(j, k, j, k) = HF_p1(j, k, j, k) +term_poten(j) +term_angular(j, l) +tmp
        end do
    end do
    call diag_sym(HF_p2, EF)
    nullify(HF_p1)
    nullify(HF_p2)
    nullify(HF_p0)
end subroutine SUB_floquet
! descale single wave function ---------------------
subroutine SUB_descale_H1
    use hamiltonian, only: weight_grid, dr_pdrho
    integer(i4) :: j
    do j = 1, N
        H1(1:N, j) = H1(1:N, j)/(weight_grid(1:N)*dr_pdrho)**0.5d0
    end do
end subroutine SUB_descale_H1
! descale floquet wave function --------------------
subroutine SUB_descale_HF
    use hamiltonian, only: weight_grid, dr_pdrho
    integer(i4) :: i, j
    do j = 1, N*(2*F +1)
        do i = -F, F
            HF(1:N, i, j) = HF(1:N, i, j)/(weight_grid(1:N)*dr_pdrho)**0.5d0
        end do
    end do
end subroutine SUB_descale_HF
! end sub-calculate --------------------------------










! ==================================================
! PROCESS
! ==================================================
! hamiltonian --------------------------------------
subroutine PROC_H(l, m)
    character(30), parameter :: form_out1 = '(1A15, 5F9.3)', form_out2 = '(1A15, 1ES15.3, 1ES15.3)'
    integer(i4), intent(in) :: l, m
    integer(i4) :: i

    if(m == 0) then
        if(allocated(H1)) deallocate(H1)
        if(allocated(E1)) deallocate(E1)
        nullify(HF)
        if(allocated(EF)) deallocate(EF)
        allocate(H1(1:N, 1:N))
        allocate(E1(1:N))
        call SUB_single(l)
        call SUB_descale_H1
        write(file_log, form_out1) "Energy: ", (E1(i), i = 1, 5)
        ! todo: matrix size option (op_mat_f)

    else if(m > 0) then
        if(allocated(H1)) deallocate(H1)
        if(allocated(E1)) deallocate(E1)
        nullify(HF)
        if(allocated(EF)) deallocate(EF)
        allocate(HF(1:N, -F:F, 1:N*(2*F +1)))
        allocate(EF(1:N*(2*F +1)))
        call SUB_floquet(l, m)
        call SUB_descale_HF
        write(file_log, form_out2) "Dressed E: ", minval(EF(:)), maxval(EF(:))
        ! todo: matrix size option (op_mat_f)
    end if
end subroutine PROC_H
! break --------------------------------------------
subroutine PROC_basis_break
    if(allocated(H1)) deallocate(H1)
    if(allocated(E1)) deallocate(E1)
    nullify(HF)
    if(allocated(EF)) deallocate(EF)
    ! todo: make break option (module basis)
end subroutine PROC_basis_break
! out ----------------------------------------------
subroutine PROC_basis_out
    if(allocated(H1)) deallocate(H1)
    if(allocated(E1)) deallocate(E1)
    nullify(HF)
    if(allocated(EF)) deallocate(EF)
end subroutine PROC_basis_out
! end process --------------------------------------
end module basis
