module inf
    use kind_type
    use global 
    implicit none
    complex(dp),    allocatable, private, protected :: outer_f(:)
    real(dp), save, allocatable, private, protected :: CS(:), PS(:, :)
contains 
    

subroutine mat_f 
    use math_const, only: i => math_i
    real   (dp) :: k, tmp1, tmp2 
    integer(i4) :: j 

    k = (2.0_dp*Mass*Scatt)**0.50
    do j = 0, L 
        tmp1 = aimag(S(j))/2.0_dp 
        tmp2 = (1.0_dp -real(S(j)))/2.0_dp 
        outer_f(j) = (2.0_dp*dble(j) +1.0_dp)/k*(tmp1 +i*tmp2)
    end do 
end subroutine mat_f 










! ==================================================
! PROCESS
! ==================================================
! cs plot ------------------------------------------
subroutine PROC_CS_plot 
    use math_const,  only: pi => math_pi, degree => math_degree
    use unit_const,  only: au_bohr
    use hamiltonian, only: coord_r, coord_theta
    use gsl_special, only: gsl_sf_legendre_Pl
    integer  (i1), parameter :: file_dcs = 101, file_tcs = 102 
    character(30), parameter :: form_cs  = '(30ES25.10)'
    character(30), parameter :: form_out = '(1A15, 5X, 1ES25.10)'
    real     (dp), parameter :: radian_to_degree = 1.0_dp/degree 
    real     (dp), parameter :: au_to_AA = au_bohr*10.0_dp**10.0_dp 
    complex  (dp) :: tmp1 
    real     (dp) :: unit_theta, unit_cs, k, tmp2 
    complex  (qp) :: sum
    integer  (i4) :: i, j 

    unit_theta = 1.0_dp 
    if(op_degree == "Y") unit_theta = radian_to_degree
    unit_cs    = 1.0_dp 
    if(op_aa == "Y") unit_cs = (au_to_AA)**2_dp

    if(.not. allocated(outer_f)) allocate(outer_f(0:L))
    call mat_f 

    open(file_tcs, file = "output/total_cs.d")
    sum = 0.0_dp 
    k   = (2.0_dp*Mass*Scatt)**0.50
    do i = 0, L 
        tmp1 = 4.0_dp*pi/k*outer_f(i)
        sum  = sum +tmp1 
        write(file_tcs, form_cs) dble(i), aimag(tmp1)
    end do 
    tmp1 = sum 
    write(file_log, form_out) "total sigma: ", aimag(tmp1)*unit_cs
    close(file_tcs)

    open(file_dcs, file = "output/diff_cs.d")
    do j = 0, ptheta 
        sum = 0.0_dp 
        do i = 0, L 
            tmp2 = cos(coord_theta(j))
            sum  = sum +outer_f(i)*gsl_sf_legendre_Pl(i, tmp2)
        end do 
        write(file_dcs, form_cs) coord_theta(j)*unit_theta, abs(sum)**2.0_dp*unit_cs
    end do 
    close(file_dcs)
    if(allocated(outer_f)) deallocate(outer_f)
end subroutine PROC_CS_plot
! end cs plot --------------------------------------
! cs achive ----------------------------------------
subroutine PROC_CS_achive(j)
    use math_const, only: pi => math_pi
    character(30), parameter  :: form_out = '(1A15, 5X, 1ES25.10)'
    integer  (i4), intent(in) :: j 
    real     (dp) :: k 
    complex  (qp) :: sum 
    complex  (dp) :: tmp 
    integer  (i4) :: i

    if(.not. allocated(CS))      allocate(CS(1:M))
    if(.not. allocated(outer_f)) allocate(outer_f(0:L))
    call mat_f 

    sum = 0.0_dp 
    k   = (2.0_dp*Mass*Scatt)**0.50
    do i = 0, L 
        tmp = 4.0_dp*pi/k*outer_f(i)
        sum = sum +tmp 
    end do 
    tmp   = sum 
    CS(j) = aimag(tmp)
    write(file_log, form_out) "total sigma: ", aimag(tmp)
    if(allocated(outer_f)) deallocate(outer_f)
end subroutine PROC_CS_achive
! end cs achive ------------------------------------
! e vs cs plot -------------------------------------
subroutine PROC_E_vs_CS_plot
    use math_const,  only: pi => math_pi
    use unit_const,  only: other_e_eV, au_hartree, au_bohr
    use hamiltonian, only: coord_E
    integer  (i1), parameter   :: file_cs  = 101
    character(30), parameter   :: form_cs  = '(30ES25.10)'
    real     (dp), parameter   :: au_to_eV = au_hartree/other_e_ev
    real     (dp), parameter   :: au_to_AA = au_bohr*10.0_dp**10.0_dp 
    real     (dp) :: unit_e, unit_cs, tmp 
    real     (qp) :: sum 
    integer  (i4) :: i, j

    unit_e  = 1.0_dp 
    if(op_ev == "Y") unit_e  = au_to_eV
    unit_cs = 1.0_dp 
    if(op_aa == "Y") unit_cs = (au_to_AA)**2_dp

    open(file_cs, file = "output/energy_vs_cs.d")
    do j = 1, M 
        sum = 0.0_dp 
        do i = 0, L 
            tmp = 4_dp*pi/(2_dp*Mass*coord_E(j))*(2_dp*dble(i) +1_dp)
            sum = sum +tmp 
        end do 
        tmp = sum 
        write(file_cs, form_cs) coord_E(j)*unit_e, CS(j)*unit_cs, tmp*unit_cs 
    end do 
    if(allocated(CS)) deallocate(CS)
    close(file_cs)
end subroutine PROC_E_vs_CS_plot
! end e vs cs plot ---------------------------------
! phase achive -------------------------------------
subroutine PROC_PS_achive(j)
    use math_const, only: pi => math_pi
    integer(i4), intent(in) :: j 
    complex(dp) :: tmp 
    integer(i4) :: i

    if(.not. allocated(PS)) allocate(PS(0:L, 1:M))
    do i = 0, L 
        tmp      = 0.5_dp*log(S(i))
        PS(i, j) = aimag(tmp)
    end do 
end subroutine PROC_PS_achive
! end phase achive ---------------------------------
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
    real     (dp) :: unit_e, tmp 
    integer  (i4) :: i, j, num 

    unit_e  = 1.0_dp 
    if(op_ev == "Y") unit_e  = au_to_eV
    if(.not. allocated(work)) allocate(work(1:M))
    do i = 0, L 
        work(:)  = PS(i, :)
        PS(i, :) = 0.0_dp 
        num = 0_i4 
        do j = 1, M 
            if(j /= 1) then 
                if(work(j -1)*work(j) < 0.0_dp) then 
                    if(work(j) > 0.0_dp) then 
                        tmp = work(j) -pi
                        if(abs(tmp -work(j -1)) < abs(work(j) -work(j -1))) num = num -1_i4 
                    else if(work(j) < 0.0_dp) then 
                        tmp = work(j) +pi
                        if(abs(tmp -work(j -1)) < abs(work(j) -work(j -1))) num = num +1_i4 
                    end if 
                end if
            end if 
            PS(i, j) = work(j) +dble(num)*pi 
        end do 
    end do 

    open(file_ph, file = "output/energy_vs_ps.d")
    do i = 0, L 
        do j = 1, M 
            write(file_ph, form_gen) i, coord_E(j)*unit_e, PS(i, j)
        end do 
        write(file_ph, form_gen)
    end do 
    close(file_ph)
    open(file_lt, file = "output/energy_vs_lt.d")
    if(.not. allocated(r)) allocate(r(1))
    do i = 0, L 
        work(:) = 0.0_dp 
        do j = 1, M 
            tmp = 0.0_dp 
            if(j < 3) then 
                tmp = diff_E(u0 = PS(i, j), up1 = PS(i, j +1), up2 = PS(i, j +2))
            else if(j > M -2) then 
                tmp = diff_E(um2 = PS(i, j -2), um1 = PS(i, j -1), u0 = PS(i, j))
            else if(j >= 3 .and. j <= M -2) then 
                tmp = diff_E(um2 = PS(i, j -2), um1 = PS(i, j -1), u0 = PS(i, j), up1 = PS(i, j +1), up2 = PS(i, j +2))
            end if 
            tmp = tmp*2.0_dp 
            write(file_lt, form_gen) i, coord_E(j)*unit_e, tmp 
            work(j) = tmp 
        end do 
        write(file_lt, form_gen)
        tmp  = maxval(work(:))
        r(:) = maxloc(work(:))
        num  = r(1)
        write(file_log, form_out) "Resonance: ", i, coord_E(num)*unit_e, 4.0_dp/tmp 
    end do 
    close(file_lt)
    if(allocated(r))    deallocate(r)
    if(allocated(work)) deallocate(work)
    if(allocated(PS))   deallocate(PS)
end subroutine PROC_E_vs_PS_plot
! end e vs phase & lifttime plot -------------------
end module inf
