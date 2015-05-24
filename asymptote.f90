module asymptote
    ! use only S matrix 
    use kind_type
    use global 
    implicit none
    complex(dp),    allocatable, private, protected :: outer_f(:)
    real(dp), save, allocatable, private, protected :: PS(:, :)
contains 
    

! ==================================================
! MATRIX
! ==================================================
! matrix f -----------------------------------------
subroutine mat_f(j)
    use hamiltonian, only: coord_E
    use math_const, only: i => math_i
    integer(i4), intent(in) :: j 
    real   (dp) :: k, tmp1, tmp2 
    integer(i4) :: h 

    k = (2.d0*Mass*coord_E(j))**0.50
    do h = 0, L 
        tmp1 = aimag(S(j, h))/2.d0 
        tmp2 = (1.d0 -real(S(j, h)))/2.d0 
        outer_f(h) = (2.d0*dble(h) +1.d0)/k*(tmp1 +i*tmp2)
    end do 
end subroutine mat_f 










! ==================================================
! PROCESS
! ==================================================
! cs plot ------------------------------------------
subroutine PROC_CS_plot 
    use math_const,  only: pi => math_pi, degree => math_degree
    use unit_const,  only: au_bohr
    use hamiltonian, only: coord_r, coord_theta, coord_E
    use fgsl, only: fgsl_sf_legendre_Pl
    integer  (i1), parameter :: file_dcs = 101, file_tcs = 102 
    character(30), parameter :: form_cs  = '(30ES25.10)'
    character(30), parameter :: form_out = '(1A15, 5X, 1ES25.10)'
    real     (dp), parameter :: radian_to_degree = 1.d0/degree 
    real     (dp), parameter :: au_to_AA = au_bohr*10.d0**10.d0 
    complex  (dp) :: tmp1 
    real     (dp) :: unit_theta, unit_cs, k, tmp2 
    complex  (qp) :: sum
    integer  (i4) :: i, j 

    unit_theta = 1.d0 
    if(op_degree == "Y") unit_theta = radian_to_degree
    unit_cs    = 1.d0 
    if(op_aa == "Y") unit_cs = (au_to_AA)**2.d0

    if(.not. allocated(outer_f)) allocate(outer_f(0:L))
    call mat_f(1_i4)

    open(file_tcs, file = "output/total_cs.d")
    sum = 0.d0 
    k   = (2.d0*Mass*coord_E(1_i4))**0.50
    do i = 0, L 
        tmp1 = 4.d0*pi/k*outer_f(i)
        sum  = sum +tmp1 
        write(file_tcs, form_cs) dble(i), aimag(tmp1)
    end do 
    tmp1 = sum 
    write(file_log, form_out) "total sigma: ", aimag(tmp1)*unit_cs
    close(file_tcs)

    open(file_dcs, file = "output/diff_cs.d")
    do j = 0, ptheta 
        sum = 0.d0 
        do i = 0, L 
            tmp2 = cos(coord_theta(j))
            sum  = sum +outer_f(i)*fgsl_sf_legendre_Pl(i, tmp2)
        end do 
        write(file_dcs, form_cs) coord_theta(j)*unit_theta, abs(sum)**2.d0*unit_cs
    end do 
    close(file_dcs)
    if(allocated(outer_f)) deallocate(outer_f)
end subroutine PROC_CS_plot
! end cs plot --------------------------------------
! e vs cs plot -------------------------------------
subroutine PROC_E_vs_CS_plot
    use math_const,  only: pi => math_pi
    use unit_const,  only: other_e_eV, au_hartree, au_bohr
    use hamiltonian, only: coord_E
    integer  (i1), parameter   :: file_cs  = 101
    character(30), parameter   :: form_cs  = '(30ES25.10)', form_out = '(1A15, 5X, 1ES25.10)'
    real     (dp), parameter   :: au_to_eV = au_hartree/other_e_ev
    real     (dp), parameter   :: au_to_AA = au_bohr*10.d0**10.d0 
    complex  (dp) :: tmp1
    complex  (qp) :: sum1  
    real     (dp) :: unit_e, unit_cs, k, tmp2
    real     (qp) :: sum2
    integer  (i4) :: i, j

    unit_e  = 1.d0 
    if(op_ev == "Y") unit_e  = au_to_eV
    unit_cs = 1.d0 
    if(op_aa == "Y") unit_cs = (au_to_AA)**2.d0

    open(file_cs, file = "output/energy_vs_cs.d")
    if(.not. allocated(outer_f)) allocate(outer_f(0:L))
    do j = 1, M 
        k = (2.d0*Mass*coord_E(j))**0.50
        outer_f(:) = 0.d0 
        call mat_f(j)
        sum1 = 0.d0 
        sum2 = 0.d0 
        do i = 0, L 
            tmp1 = 4.d0*pi/k*outer_f(i)
            sum1 = sum1 +tmp1 
            tmp2 = 4.d0*pi/(2.d0*Mass*coord_E(j))*(2.d0*dble(i) +1.d0)
            sum2 = sum2 +tmp2 
        end do 
        tmp1  = sum1 
        tmp2  = sum2 
        write(file_cs, form_cs) coord_E(j)*unit_e, aimag(tmp1)*unit_cs, tmp2*unit_cs 
    end do 
    if(allocated(outer_f)) deallocate(outer_f)
    close(file_cs)
end subroutine PROC_E_vs_CS_plot
! end e vs cs plot ---------------------------------
! e vs phase & lifttime plot -----------------------
subroutine PROC_E_vs_PS_plot
    use math_const,  only: pi => math_pi
    use unit_const,  only: other_e_eV, au_hartree
    use hamiltonian, only: coord_E, diff_E 
    integer  (i1), parameter   :: file_ph  = 101, file_lt = 102 
    character(30), parameter   :: form_gen = '(1I5, 2ES25.10)', form_out = '(1A15, 1I15, 2ES15.3)'
    real     (dp), parameter   :: au_to_eV = au_hartree/other_e_ev
    real     (dp), allocatable :: work(:)
    integer  (i4), allocatable :: r(:)
    complex  (dp) :: tmp1 
    real     (dp) :: unit_e, tmp2
    integer  (i4) :: i, j, num 

    unit_e  = 1.d0  
    if(op_ev == "Y") unit_e  = au_to_eV

    if(.not. allocated(PS)) allocate(PS(1:M, 0:L))
    do i = 0, L 
        do j = 1, M 
            tmp1     = 0.5d0*log(S(j, i))
            PS(j, i) = aimag(tmp1)
        end do 
    end do 

    if(.not. allocated(work)) allocate(work(1:M))
    do i = 0, L 
        work(:)  = PS(:, i)
        PS(:, i) = 0.d0  
        num = 0_i4 
        do j = 1, M 
            if(j /= 1) then 
                if(work(j -1)*work(j) < 0.d0) then 
                    if(work(j) > 0.d0) then 
                        tmp2 = work(j) -pi
                        if(abs(tmp2 -work(j -1)) < abs(work(j) -work(j -1))) num = num -1_i4 
                    else if(work(j) < 0.d0) then 
                        tmp2 = work(j) +pi
                        if(abs(tmp2 -work(j -1)) < abs(work(j) -work(j -1))) num = num +1_i4 
                    end if 
                end if
            end if 
            PS(j, i) = work(j) +dble(num)*pi 
        end do 
        tmp2 = (maxval(PS(:, i)) +minval(PS(:, i)))/2.d0/pi
        num  = int(tmp2)
        if(tmp2 < 0.d0) num = num -1_i4
        PS(:, i) = PS(:, i) -dble(num)*pi 
    end do 

    open(file_ph, file = "output/energy_vs_ps.d")
    do i = 0, L 
        do j = 1, M 
            write(file_ph, form_gen) i, coord_E(j)*unit_e, PS(j, i)
        end do 
        write(file_ph, form_gen)
        write(file_ph, form_gen)
    end do 
    close(file_ph)
    open(file_lt, file = "output/energy_vs_lt.d")
    if(.not. allocated(r)) allocate(r(1))
    do i = 0, L 
        work(:) = 0.d0 
        do j = 1, M 
            tmp2 = 0.d0 
            if(j < 3) then 
                tmp2 = diff_E(u0 = PS(j, i), up1 = PS(j +1, i), up2 = PS(j +2, i))
            else if(j > M -2) then 
                tmp2 = diff_E(um2 = PS(j -2, i), um1 = PS(j -1, i), u0 = PS(j, i))
            else if(j >= 3 .and. j <= M -2) then 
                tmp2 = diff_E(um2 = PS(j -2, i), um1 = PS(j -1, i), u0 = PS(j, i), up1 = PS(j +1, i), up2 = PS(j +2, i))
            end if 
            tmp2 = tmp2*2.d0 
            write(file_lt, form_gen) i, coord_E(j)*unit_e, tmp2 
            work(j) = tmp2 
        end do 
        write(file_lt, form_gen)
        write(file_lt, form_gen)
        tmp2 = maxval(work(:))
        r(:) = maxloc(work(:))
        num  = r(1)
        write(file_log, form_out) "Resonance: ", i, coord_E(num)*unit_e, 4.d0/tmp2 
    end do 
    close(file_lt)
    if(allocated(r))    deallocate(r)
    if(allocated(work)) deallocate(work)
    if(allocated(PS))   deallocate(PS)
end subroutine PROC_E_vs_PS_plot
! end e vs phase & lifttime plot -------------------
end module asymptote
