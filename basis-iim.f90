module basis
    use kind_type 
    use global 
    implicit none
    real(dp), save, allocatable, protected :: H1(:, :), P1(:, :), E1(:)
    real(dp), save, pointer,     protected :: H(:, :, :), E(:)
contains


! ==================================================
! FUNCTION
! ==================================================
! kinetic term -------------------------------------
function term_kinet(i, j)
    use hamiltonian, only: delta_rho, dr_p_drho
    integer(i4), intent(in) :: i, j 
    real   (dp) :: term_kinet
    term_kinet = -delta_rho(i, j)/dr_p_drho**2.d0 &
                    /(2.d0*Mass)
end function term_kinet 
! potential term -----------------------------------
function term_poten(i)
    use hamiltonian, only: coord_r, poten_r 
    integer(i4), intent(in) :: i
    real   (dp) :: term_poten
    term_poten = poten_r(coord_r(i)) & 
                    *Charge
end function term_poten 
! angular term -------------------------------------
function term_angular(i, l)  
    use hamiltonian, only: coord_r
    integer(i4), intent(in) :: i, l 
    real   (dp) :: term_angular
    term_angular = dble(l)*(dble(l) +1.d0)/coord_r(i)**2.d0 & 
                    /(2.d0*Mass)
end function term_angular
! floquet term -------------------------------------
function term_floquet(n)
    use hamiltonian, only: coord_r
    integer(i4), intent(in) :: n
    real   (dp) :: term_floquet
    term_floquet = -dble(n)*Freq
end function term_floquet
! dipole term --------------------------------------
function term_dipole(i)
    use hamiltonian, only: coord_r
    integer(i4), intent(in) :: i
    real   (dp) :: term_dipole
    term_dipole = 0.5d0*Amp*coord_r(i) & 
                    *Charge
end function term_dipole
! end funciton -------------------------------------










! ==================================================
! SUB-CALCULATE
! ==================================================
! single hamiltonian -------------------------------
subroutine SUB_single(l)
    use linear, only: diag_sym 
    integer(i4), intent(in) :: l 
    integer(i4) :: i, j 
    H1(:, :) = 0.d0 
    E1(:) = 0.d0 
    do j = 1, N 
        do i = 1, N 
            H1(i, j) = term_kinet(i, j)
        end do 
        H1(j, j) = H1(j, j) +term_poten(j) +term_angular(j, l)
    end do 
    P1(:, :) = H1(:, :)
    call diag_sym(P1(:, :), E1(:))
end subroutine SUB_single
! end single hamiltonian ---------------------------
! floquet hamiltonian ------------------------------
subroutine SUB_iim
    use linear, only: solve_band_vec
    real(dp) :: lambda1, lambda2, dlambda
    real(qp) :: sum 
    real(dp), pointer :: & 
        H_p1(:, :, :), H_p2(:), E_p(:, :), &
        HF(:, :, :), HF_p1(:, :), HF_p2(:), &
        H2(:, :, :), H2_p1(:, :), H2_p2(:), & 
        x(:, :), x_p(:), y(:, :), y_p(:)
    integer(i4) :: nb, nf, i, j 

    nullify(H_p1)
    nullify(H_p2)
    nullify(E_p)
    nullify(HF)
    nullify(HF_p1)
    nullify(HF_p2)
    nullify(H2)
    nullify(H2_p1)
    nullify(H2_p2)
    nullify(x)
    nullify(x_p)
    nullify(y)
    nullify(y_p)
    allocate(H_p1 (1:N*(2*F +1), 1:N, -F:F))
    allocate(H_p2 (1:N*(2*F +1)*N*(2*F +1)))
    allocate(E_p  (1:N, -F:F))
    allocate(HF   (-N:N, 1:N, -F:F))
    allocate(HF_p1(-N:N, 1:N*(2*F +1)))
    allocate(HF_p2(1:(2*N +1)*N*(2*F +1)))
    allocate(H2   (-2*N:N, 1:N, -F:F))
    allocate(H2_p1(-2*N:N, 1:N*(2*F +1)))
    allocate(H2_p2(1:(3*N +1)*N*(2*F +1)))
    allocate(x    (1:N, -F:F))
    allocate(x_p  (1:N*(2*F +1)))
    allocate(y    (1:N, -F:F))
    allocate(y_p  (1:N*(2*F +1)))
    H    (1:N, -F:F, 1:N*(2*F +1)) => H_p2 (1:N*(2*F +1)*N*(2*F +1))
    H_p1 (1:N*(2*F +1), 1:N, -F:F) => H_p2 (1:N*(2*F +1)*N*(2*F +1))
    E_p  (1:N, -F:F)               => E    (1:N*(2*F +1))
    HF   (-N:N, 1:N, -F:F)         => HF_p2(1:(2*N +1)*N*(2*F +1))
    HF_p1(-N:N, 1:N*(2*F +1))      => HF_p2(1:(2*N +1)*N*(2*F +1))
    H2   (-2*N:N, 1:N, -F:F)       => H2_p2(1:(3*N +1)*N*(2*F +1))
    H2_p1(-2*N:N, 1:N*(2*F +1))    => H2_p2(1:(3*N +1)*N*(2*F +1))
    x    (1:N, -F:F)               => x_p  (1:N*(2*F +1))
    y    (1:N, -F:F)               => y_p  (1:N*(2*F +1))

    HF(:, :, :) = 0.d0 
    do j = 1, N 
        do i = 1, N 
            HF(i -j, j, -F:F) = H1(i, j)
        end do 
    end do 
    do i = -F, F 
        HF(0, :, i) = HF(0, :, i) +term_floquet(i)
    end do 
    do i = 1, N 
        HF( N, i, -F:F) = term_dipole(i)
        HF(-N, i, -F:F) = term_dipole(i)
    end do 

    E(:) = 0.d0 
    H(:, :, :) = 0.d0 
    do nf = -F, F 
        do nb = 1, N 
            lambda2 = E1(nb) +term_floquet(nf)
            x(:, :) = 0.d0 
            do i = -F, F 
                x(1:N, i) = P1(1:N, nb)
            end do 
            dlambda = 1.d0 
!             do while(dlambda > 1.d-16)
            do while(dlambda > 1.d-12)
                lambda1 = lambda2
                H2(:, :, :) = 0.d0 
                H2(-N:N, :, :) = HF(-N:N, :, :)
                H2(0, :, :) = H2(0, :, :) -lambda1
                y(:, :) = x(:, :)
                call solve_band_vec(H2_p1, y_p)
                sum = 0.d0 
                do i = 1, N*(2*F +1)
                    sum = sum +x_p(i)*y_p(i)
                end do 
                lambda2 = lambda1 +1.d0/sum 
                dlambda = abs(1.d0/sum/lambda1)
                sum = 0.d0 
                do i = 1, N*(2*F +1)
                    sum = sum +y_p(i)**2.d0 
                end do 
                x(:, :) = y(:, :)/sum**0.5d0 
            end do 
            E_p(nb, nf) = lambda2
            H_p1(:, nb, nf) = x_p(:)
        end do 
    end do 

    nullify(H_p1)
    nullify(H_p2)
    nullify(E_p)
    nullify(HF)
    nullify(HF_p1)
    nullify(HF_p2)
    nullify(H2)
    nullify(H2_p1)
    nullify(H2_p2)
    nullify(x)
    nullify(x_p)
    nullify(y)
    nullify(y_p)
end subroutine SUB_iim
! end floquet hamiltonian --------------------------
! descale single wave function ---------------------
subroutine SUB_descale_H1
    use hamiltonian, only: coord_weight, dr_p_drho
    integer(i4) :: i0, i, j 
    H1(:, :) = P1(:, :)
    if(allocated(P1)) deallocate(P1)
    i0 = 1 
    if(size(H1(:, 1)) == 1) i0 = N 
    do j = 1, N 
        do i = i0, N 
            H1(i, j) = H1(i, j) &
                        /(coord_weight(i)*dr_p_drho)**0.5_dp
        end do 
    end do 
end subroutine SUB_descale_H1
! end descale single wave function -----------------
! descale floquet wave function --------------------
subroutine SUB_descale_H
    use hamiltonian, only: coord_weight, dr_p_drho
    integer(i4) :: i0, i1, i2, j 
    i0 = 1 
    if(size(H(:, 0, 1)) == 1) i0 = N 
    do j = 1, N*(2*F +1)
        do i2 = -F, F 
            do i1 = i0, N 
                H(i1, i2, j) = H(i1, i2, j) &
                                    /(coord_weight(i1)*dr_p_drho)**0.5_dp
            end do 
        end do 
    end do 
end subroutine SUB_descale_H
! end descale floquet wave function ----------------










! ==================================================
! PROCESS
! ==================================================
! hamiltonian --------------------------------------
subroutine PROC_H(l) 
    character(30), parameter :: form_out = '(1A15, 5F9.3)'
    integer(i4), intent(in) :: l 
    integer(i4) :: i

    if(allocated(H1)) deallocate(H1)
    if(allocated(P1)) deallocate(P1)
    if(allocated(E1)) deallocate(E1)
    allocate(H1(1:N, 1:N))
    allocate(P1(1:N, 1:N))
    allocate(E1(1:N))
    call SUB_single(l)
    write(file_log, form_out) "Energy: ", (E1(i), i = 1, 5) 

    nullify(H)
    nullify(E)
    allocate(H(1:N, -F:F, 1:N*(2*F +1)))
    allocate(E(1:N*(2*F +1)))
    call SUB_iim
    write(file_log, form_out) "Dressed E: ", E(1), E(N*(2*F +1))

!     if(.not. allocated(H)) then 
!         if(op_mat_f == 1) then 
!             allocate(H(1:N*(2*F +1), -F:F, 1:N))
!         else if(op_mat_f == 0) then 
!             allocate(H(1:N*(2*F +1), -F:F, N:N))
!         end if 
!     end if 
    call SUB_descale_H1
    call SUB_descale_H

!     if(allocated(H1)) deallocate(H1)
!     if(allocated(P1)) deallocate(P1)
!     if(allocated(E1)) deallocate(E1)
!     nullify(H)
!     nullify(E)
end subroutine PROC_H
! end hamiltonian ----------------------------------
! basis plot ---------------------------------------
subroutine PROC_basis_plot(l)
    use hamiltonian, only: coord_r
    integer  (i1), parameter  :: file_psi = 101, file_ene = 102
    character(30), parameter  :: & 
        form_gen = '(2X, 25ES25.10)', & 
        form_tit = '("# ", 1A75)', &
        form_sub = '("# ", 2A25, 1A50)', &
        form_lin = '("# ", 2A25, 20ES25.10)'
    integer  (i4), intent(in) :: l 
    integer  (i4) :: i1, i2, j1, j2, loc(1)
    character (3) :: ch 

    write(ch, '(I3.3)') l 

        open(file_psi, file = "output/basis_u_"//ch//".d")
        open(file_ene, file = "output/basis_energy_"//ch//".d")
        write(file_psi, form_tit) "==========================================================================="
        write(file_psi, form_tit) "BASIS FUNCTION OF INNER REIGON IN THE ABSENCE OF THE FIELD"
        write(file_psi, form_tit) "==========================================================================="  
        write(file_psi, form_sub) "", " RADIAL COORDINATE ", " WAVE FUNCTION PER R (ENERGY) "
        write(file_psi, form_lin) " ----------------------- ", " ----------------------- ", &
                                    (E1(j1), j1 = 1, 10), (E1(j1), j1 = N -9, N)
        write(file_psi, *)
        write(file_psi, form_gen) 0.d0, 0.d0, (0.d0, j1 = 1, 20)
        do i1 = 1, N 
            write(file_psi, form_gen) 0.d0, coord_r(i1), & 
                (H1(i1, j1), j1 = 1, 10), (H1(i1, j1), j1 = N -9, N)
        end do 
        write(file_ene, form_tit) "==========================================================================="
        write(file_ene, form_tit) "ENERGY LEVEL OF INNER REGION IN THE ABSENCE OF THE FIELD"
        write(file_ene, form_tit) "==========================================================================="
        write(file_ene, form_sub) " STATE NUMBER ", " ENERGY "
        write(file_ene, form_sub) " ----------------------- ", " ----------------------- "
        write(file_ene, *)
        do j1 = 1, N 
            write(file_ene, form_gen) dble(j1), E(j1)
        end do 
        close(file_psi)
        close(file_ene)
    
    if(.not. F == 0) then 
        loc(:) = minloc(E(:))
        j1 = loc(1)
        loc(:) = maxloc(E(:))
        j2 = loc(1)
        open(file_psi, file = "output/dressed_u_"//ch//".d")
        open(file_ene, file = "output/dressed_energy_"//ch//".d")
        write(file_psi, form_tit) "==========================================================================="
        write(file_psi, form_tit) "DRESSED BASIS FUNCTION OF INNER REIGON"
        write(file_psi, form_tit) "==========================================================================="  
        write(file_psi, form_sub) " PHOTON SPACE ", " RADIAL COORDINATE ", " WAVE FUNCTION PER R (ENERGY) "
        write(file_psi, form_lin) " ----------------------- ", " ----------------------- ", &
                                    E(j1), E(j2)
        write(file_psi, *)
        do i2 = -F, F 
            write(file_psi, form_gen) dble(i2), 0.d0, 0.d0, 0.d0 
            do i1 = 1, N 
                write(file_psi, form_gen) dble(i2), coord_r(i1), & 
                    H(i1, i2, j1), H(i1, i2, j2)
            end do 
            write(file_psi, *)
        end do 
        write(file_ene, form_tit) "==========================================================================="
        write(file_ene, form_tit) "DRESSED ENERGY LEVEL OF INNER REGION"
        write(file_ene, form_tit) "==========================================================================="
        write(file_ene, form_sub) " STATE NUMBER ", " ENERGY "
        write(file_ene, form_sub) " ----------------------- ", " ----------------------- "
        write(file_ene, *)
        do j1 = 1, (2*F +1)*N 
            write(file_ene, form_gen) dble(j1), E(j1)
        end do 
        close(file_psi)
        close(file_ene)
    end if 
end subroutine PROC_basis_plot
! end basis plot -----------------------------------
! break --------------------------------------------
subroutine PROC_basis_break
!     if(allocated(H1)) deallocate(H1)
!     if(allocated(P1)) deallocate(P1)
!     if(allocated(E1)) deallocate(E1)
!     if(op_mat_h == 0) nullify(H)
!     if(op_mat_e == 0) nullify(E)
end subroutine PROC_basis_break
! end break ----------------------------------------
! out ----------------------------------------------
subroutine PROC_basis_out
    if(allocated(H1)) deallocate(H1)
    if(allocated(P1)) deallocate(P1)
    if(allocated(E1)) deallocate(E1)
    nullify(H)
    nullify(E)
end subroutine PROC_basis_out
! end out ------------------------------------------
end module basis
