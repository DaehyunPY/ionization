module hamiltonian 
    use kind_type
    use global 
    implicit none
    ! coordnation
    real(dp), save, protected :: dr_p_drho
    real(dp), save, protected, private :: dtheta, dE 
    real(dp), save, allocatable, protected :: coord_weight(:)
    real(dp), save, allocatable, protected, private :: coord_rho(:), coord_dshape(:, :)
    ! potiontial 
    real(dp), save, private, protected :: &
        Z, alphab, cutoff, &   ! general 
        z0, z1, z2, &          ! II-1 GELTMAN
        pD, pH, pdelta, &      ! II-2 GREEN & YUKAWA
        pA, pB, alpha, beta, & ! II-3 KYLSTRA & BUCKINGHAM
        alpha1, alpha2, alpha3 ! III  CHEN & NOUS 
    integer(i1), save, private, protected :: op_poten 
contains


! ==================================================
! COORDNATION
! ==================================================
! r ------------------------------------------------
function coord_r(i) 
    integer(i4), intent(in) :: i
    real   (dp) :: coord_r
    coord_r = dr_p_drho*coord_rho(i) +Ra/2.d0 
end function coord_r
! theta --------------------------------------------
function coord_theta(i)
    integer(i4), intent(in) :: i
    real   (dp) :: coord_theta
    coord_theta = dtheta*dble(i)
end function coord_theta
! E ------------------------------------------------
function coord_E(i)
    integer(i4), intent(in) :: i
    real   (dp) :: coord_E 
    coord_E = dE*dble(i)
end function coord_E
! end coordination ---------------------------------


! ==================================================
! OPERATOR 
! ==================================================
! diff ---------------------------------------------
function diff_E(um2, um1, u0, up1, up2)
    use math_const, only: pi => math_pi 
    real(dp), intent(in) :: u0 
    real(dp), intent(in), optional :: um2, um1, up1, up2 
    real(dp) :: diff_E, tmp 
    tmp = 0.d0 
    if(present(um2) .and. present(up2)) then 
        tmp = +1.d0/12.d0 * um2 &
              -8.d0/12.d0 * um1 &
              +8.d0/12.d0 * up1 &
              -1.d0/12.d0 * up2 
    else if(present(um2) .and. (.not. present(up2))) then 
        tmp = +1.d0/2.d0  * um2 &
              -4.d0/2.d0  * um1 &
              +3.d0/2.d0  * u0 
    else if((.not. present(um2)) .and. present(up2)) then 
        tmp = -3.d0/2.d0  * u0  &
              +4.d0/2.d0  * up1 &
              -1.d0/2.d0  * up2 
    else if((.not. present(um2)) .and. (.not. present(up2))) then 
        stop "Error #999: function diff_E"
    end if 
    diff_E = tmp/dE 
end function diff_E
! delta --------------------------------------------
function delta_rho(i, j)
    use math_const, only: pi => math_pi 
    integer(i4), intent(in) :: i, j 
    real   (qp) :: sum 
    real   (dp) :: delta_rho
    integer(i4) :: k, n  
    n   = dble(size(coord_rho(:))) -1 
    sum = 0.d0 
    do k = 0, n 
        sum = sum -coord_dshape(i, k)*coord_dshape(j, k)
    end do 
    delta_rho = sum 
end function delta_rho
! potential ----------------------------------------
function poten_r(r)
    real(dp), intent(in) :: r
    real(dp) :: poten_r, stat, pol, tmp 
    if(op_poten == 0) then 
!         tmp     = 7.75d0*r**2.d0*exp(-r) ! example 2.6.1
        tmp     = -2.5d0 ! example P G Burke 
        poten_r = tmp/Charge
    else if(op_poten == 1) then 
        poten_r = Z/r 
    else if(op_poten == 2) then 
        stat    = -exp(-2.d0*z0*r)*(z1 +z2/r)
        pol     = -alphab/(2.d0*r**4.d0)*(1.d0 -exp(-r))**6.d0 
        tmp     = stat +pol 
        poten_r = tmp/(-1.d0)
    else if(op_poten == 3) then 
        stat    = -(Z/pH)*(exp(-r/pD)/r)*(1.d0 +(pH -1.d0)*exp(-r*pH/pD))
        pol     = -alphab/(2.d0*(r**2.d0 +pdelta**2.d0)**2.d0)
        tmp     = stat +pol 
        poten_r = tmp/(-1.d0)
    else if(op_poten == 4) then 
        stat    = -Z*( &
                    +pA**2.d0/(4.d0*alpha**3.d0)*exp(-2.d0*alpha*r)*(alpha +1.d0/r) &
                    +4.d0*pA*pB/(alpha +beta)**3.d0*exp(-(alpha +beta)*r)*((alpha +beta)/2.d0 +1.d0/r) & 
                    +pB**2.d0/(4.d0*beta**3.d0)*exp(-2.d0*beta*r)*(beta +1.d0/r))
        pol     = -alphab/(2.d0*(cutoff**2.d0 +r**2.d0)**2.d0)
        tmp     = stat +pol 
        poten_r = tmp/(-1.d0)
    else if(op_poten == 5) then 
        stat    = -(Z/r)*exp(-alpha1*r) -alpha2*exp(-alpha3*r)
        pol     = -alphab/(2.d0*r**4.d0) * (1.d0 -exp(-(r/cutoff)**3.d0))**2.d0 
        tmp     = stat +pol 
        poten_r = tmp/(-1.d0)
    end if 
end function poten_r
! end operator -------------------------------------


! ==================================================
! TEST 
! ==================================================
! test potential -----------------------------------
! subroutine TEST_poten
!     integer  (i1), parameter :: file_poten = 101
!     character(30), parameter :: form_poten = '(30ES25.10)'
!     integer  (i4) :: i 
!     open(file_poten, file = "output/poten.d")
!     do i = 1, N
!         write(file_poten, form_poten) coord_r(i), poten_r(coord_r(i))*Charge
!     end do
!     close(file_poten)
! end subroutine TEST_poten
! test coordination system -------------------------
! subroutine TEST_coord
!     use fgsl, only: fgsl_sf_legendre_Pl
!     integer  (i1), parameter :: file_coord = 101
!     character(30), parameter :: form_coord = '(30ES25.10)'
!     integer  (i4) :: i 
!     open(file_coord, file = "output/coord.d")
!     write(file_coord, form_coord) coord_r(0_i4), coord_rho(0), coord_weight(0)
!     do i = 1, N -1 
!         write(file_coord, form_coord) coord_r(i), coord_rho(i), coord_weight(i), & 
!             N*(coord_rho(i)*fgsl_sf_legendre_Pl(N, coord_rho(i)) -fgsl_sf_legendre_Pl(N -1_i4, coord_rho(i))) &
!                 /(coord_rho(i)**2.d0 -1.d0)
!     end do 
!     write(file_coord, form_coord) coord_r(N), coord_rho(N), coord_weight(N)
!     close(file_coord)
! end subroutine TEST_coord
! end test -----------------------------------------










! ==================================================
! PROCESS
! ==================================================
! input --------------------------------------------
subroutine PROC_input
    use math_const, only: pi => math_pi
    use unit_const, only: other_e_eV, au_hartree
    character(60),  parameter :: & 
        form_laser = '(4/, 2(45X, 1F15.8, /), /)', &
        form_part  = '(4/, 3(45X, 1F15.8, /), 1/, 2(45X, 1F15.8, /), /)', &
        form_poten = '(4/, 1(45X, 1I15, /), 3(45X, 1F15.8, /))'
    character(30),  parameter :: &
        form_p1 = '(20/)', &
        form_p2 = '( 2/, 3(45X, 1F15.8, /), 15/)', &
        form_p3 = '( 7/, 3(45X, 1F15.8, /), 10/)', &
        form_p4 = '(11/, 4(45X, 1F15.8, /),  5/)', &
        form_p5 = '(16/, 3(45X, 1F15.8, /),  1/)'
    character(60),  parameter :: &
        form_cal = '(4/, 2(45X, 1F15.8, /), 6(45X, 1I15, /), /)'
    character(120), parameter :: &
        form_opt   = & 
            '(6/, 3(45X, 6X, 1I15, /), /, 2(45X, 6X, 1I15, /), /, 3(45X, 6X, 1I15, /), 3/, 2(45X, 6X, 1I15, /))'
    real(dp),  parameter :: eV_to_au = other_e_ev/au_hartree
    real(dp) :: tmp1, tmp2 

    open(file_input, file = "input.d")
    ! laser field ----------------------------------
    read(file_input, form_laser) Amp, Freq
    if(.not. Amp  >  0.d0) stop "Error #101: check 'LASER FIELD - AMPLITUDE'"
    if(.not. Freq >  0.d0) stop "Error #102: check 'LASER FIELD - FREQUENCY'"
    ! particle -------------------------------------
    read(file_input, form_part ) Mass, Charge, Spin, tmp1, tmp2 
    if(.not. Mass >  0.d0)  stop "Error #201: check 'PARTICLE - MASS'"    
    if(.not. Spin == 0.5d0) stop "Error #202: check 'PARTICLE - SPIN'"
    if(tmp1 >= 0.d0 .and. tmp2 >= 0.d0) stop "Error #203: check 'PARTICLE - SCATTERING ENERGY'"
    if(tmp1 <  0.d0 .and. tmp2 <  0.d0) stop "Error #204: check 'PARTICLE - SCATTERING ENERGY'"
    if(tmp1 >= 0.d0) Scatt = tmp1 
    if(tmp2 >= 0.d0) Scatt = tmp2*eV_to_au
    ! potential ------------------------------------
    read(file_input, form_poten) op_poten, Z, alphab, cutoff
    if(.not. (op_poten >= 0 .and. op_poten <= 6)) stop "Error #301: check 'POTENTIAL - TYPE'"
    if(op_poten == 0) then 
        read(file_input, form_p1)         
    else if(op_poten == 1) then 
        read(file_input, form_p1) 
    else if(op_poten == 2) then 
        read(file_input, form_p2) z0, z1, z2
    else if(op_poten == 3) then 
        read(file_input, form_p3) pD, pH, pdelta
        if(pD < 0.d0 .and. pH < 0.d0) stop "Error #302: check 'POTENTIAL - COEFFICIENT D, H'"
        if(pD < 0.d0)     pD     = pH/Z**0.4d0 
        if(pH < 0.d0)     pH     = pD*Z**0.4d0 
        if(pdelta < 0.d0) pdelta = (alphab/(2.d0*Z**(1.d0/3.d0)))**(0.25d0)
    else if(op_poten == 4) then 
        read(file_input, form_p4) pA, pB, alpha, beta
    else if(op_poten == 5) then 
        read(file_input, form_p5) alpha1, alpha2, alpha3
    end if 
    ! calculation ----------------------------------
    read (file_input, form_cal) Ra, Rap, L, N, F, M, pr, ptheta
    if(.not. Ra  >  0.d0) stop "Error #401: check 'CALCULATION - INNER REGION SIZE'"
    if(.not. Rap >= Ra)   stop "Error #402: check 'CALCULATION - OUTER REGION SIZE'"
    if(.not. L  >= 0) stop "Error #403: check 'CALCULATION - MAXIUM OF ANGULAR MOMANTUM L'"
    if(.not. N  >  0) stop "Error #404: check 'CALCULATION - GRID NUMBER OF r COORDINATES'"
    if(.not. F  >= 0) stop "Error #405: check 'CALCULATION - FLOQUET NUMBER'"
    if(.not. M  >  0) stop "Error #406: check 'CALCULATION - GRID NUMBER OF Energy COORDINATES'"
    if(.not. pr >  0) stop "Error #407: check 'CALCULATION - 3D PLOT NUMBER OF r COORDINATES'"
    if(.not. ptheta > 0) stop "Error #408: check 'CALCULATION - 3D PLOT NUMBER OF theta COORDINATES'"
    if(pr > N) pr = N 
    dtheta = pi/dble(ptheta)
    dE     = Scatt/dble(M) 
    ! option ---------------------------------------
    read (file_input, form_opt) & 
        op_ev, op_degree, op_aa, & 
        op_basis, op_bound, & 
        op_dcs, op_inner, op_outer, &
        op_tcs, op_ps
    if(.not.(op_ev     == 1 .or. op_ev     == 0)) stop "Error #501: check 'OPTION - USE ENERGY UNIT eV'"
    if(.not.(op_degree == 1 .or. op_degree == 0)) stop "Error #502: check 'OPTION - USE ANGULAR UNIT degree'"
    if(.not.(op_aa     == 1 .or. op_aa     == 0)) stop "Error #503: check 'OPTION - USE CROSS SECTION UNIT A^2'"
    if(.not.(op_basis  == 1 .or. op_basis  == 0)) stop "Error #504: check 'OPTION - BASIS FUNCTION'"
    if(.not.(op_bound  == 1 .or. op_bound  == 0)) stop "Error #505: check 'OPTION - BOUND STATES ENERGY & FUNCTION'"
    if(.not.(op_dcs    == 1 .or. op_dcs    == 0)) stop "Error #506: check 'OPTION - CROSS SECTION FUNCTION'"
    if(.not.(op_inner  == 1 .or. op_inner  == 0)) stop "Error #507: check 'OPTION - INNER REGION WAVE FUNCTION'"
    if(.not.(op_outer  == 1 .or. op_outer  == 0)) stop "Error #508: check 'OPTION - OUTER REGION WAVE FUNCTION'"
    if(.not.(op_tcs    == 1 .or. op_tcs    == 0)) stop "Error #509: check 'OPTION - KINETIC ENERGY vs CROSS SECTION'"
    if(.not.(op_ps     == 1 .or. op_ps     == 0)) stop "Error #510: check 'OPTION - KINETIC ENERGY vs PHASE'"
    close(file_input) 
    open (file_log, file = "output/log.d")
    if(.not. (op_tcs == 1 .or. op_ps == 1)) then 
        M  = 1
        dE = Scatt
    else if(op_tcs == 1 .or. op_ps == 1) then 
        op_basis = 0
        op_bound = 0
        op_dcs   = 0
        op_inner = 0 
        op_outer = 0
    end if 
    ! matrix ---------------------------------------
    if(op_basis == 1) then 
        op_mat_f = 1 ! use full size H 
    end if 
    if(op_inner == 1) then 
        op_mat_f = 1 ! use full size H 
        op_mat_h = 1 ! use H, E, R, K, A matrix 
        op_mat_e = 1 ! not use S matrix 
        op_mat_r = 1
        op_mat_k = 1
        op_mat_a = 1
    end if 
    if(op_outer == 1) then 
        op_mat_k = 1 ! use K, A 
        op_mat_a = 1 ! not use H, E, R, S
    end if 
    if(op_dcs == 1 .or. op_tcs == 1 .or. op_ps == 1) then 
        op_mat_s = 1 ! use only S matrix 
    end if 
end subroutine PROC_input
! end input ----------------------------------------
! information --------------------------------------
subroutine PROC_inform
    use unit_const, only: other_e_eV, au_hartree
    character(30), parameter :: form_out = '(1A15, 1ES15.3)'
    real     (dp), parameter :: au_to_eV = au_hartree/other_e_ev

    write(file_log, *) "================================================================="
    write(file_log, *) "LASER FIELD"
    write(file_log, *) "================================================================="
    write(file_log, *) " -------------------------------------------  -------------------"
    write(file_log, *) " AMPLITUDE                              [au] ", Amp 
    write(file_log, *) " FREQUENCY                              [au] ", Freq
    write(file_log, *) " - "
    write(file_log, *) " - "
    write(file_log, *) "================================================================="
    write(file_log, *) "PARTICLE: ELECTRON"
    write(file_log, *) "================================================================="
    write(file_log, *) " -------------------------------------------  ------------------ "
    write(file_log, *) " MASS                                   [au] ", Mass 
    write(file_log, *) " CHARGE                                 [au] ", Charge
    write(file_log, *) " SPIN                                    [1] ", Spin 
    write(file_log, *) " KINETIC ENERGY                         [au] ", Scatt
    write(file_log, *) "                                        [eV] ", Scatt*au_to_eV
    write(file_log, *) " - "
    write(file_log, *) " - "
    if(op_poten == 0) then 
        write(file_log, *) "================================================================="
        write(file_log, *) "POTENTIAL: EXAMPLE"
        write(file_log, *) "================================================================="
        write(file_log, *) " -------------------------------------------  ------------------ "
        write(file_log, *) " - "
        write(file_log, *) " - "
    else if(op_poten == 1) then 
        write(file_log, *) "================================================================="
        write(file_log, *) "POTENTIAL: COULOMB"
        write(file_log, *) "================================================================="
        write(file_log, *) " -------------------------------------------  ------------------ "
        write(file_log, *) " ELEMENT NUMBER                          [1] ", Z 
        write(file_log, *) " - "
        write(file_log, *) " - "
    else if(op_poten == 2) then 
        write(file_log, *) "================================================================="
        write(file_log, *) "POTENTIAL: GELTMAN"
        write(file_log, *) "================================================================="
        write(file_log, *) " -------------------------------------------  ------------------ "
        write(file_log, *) " POLARIZABILITY                         [au] ", alphab
        write(file_log, *) " COEFFICIENT z0                         [au] ", z0
        write(file_log, *) " COEFFICIENT z1                         [au] ", z1
        write(file_log, *) " COEFFICIENT z2                         [au] ", z2
        write(file_log, *) " - "
        write(file_log, *) " - "
    else if(op_poten == 3) then 
        write(file_log, *) "================================================================="
        write(file_log, *) "POTENTIAL: GREEN & YUKAWA"
        write(file_log, *) "================================================================="
        write(file_log, *) " -------------------------------------------  ------------------ "
        write(file_log, *) " ELEMENT NUMBER                          [1] ", Z 
        write(file_log, *) " COEFFICIENT D                          [au] ", pD 
        write(file_log, *) " COEFFICIENT H                          [au] ", pH 
        write(file_log, *) " COEFFICIENT delta                      [au] ", pdelta
        write(file_log, *) " - "
        write(file_log, *) " - "
    else if(op_poten == 4) then 
        write(file_log, *) "================================================================="
        write(file_log, *) "POTENTIAL: KYLSTRA & BUCKINGHAM"
        write(file_log, *) "================================================================="
        write(file_log, *) " -------------------------------------------  ------------------ "
        write(file_log, *) " ELEMENT NUMBER                          [1] ", Z 
        write(file_log, *) " POLARIZABILITY                         [au] ", alphab
        write(file_log, *) " CUTOFF                                 [au] ", cutoff
        write(file_log, *) " COEFFICIENT A                          [au] ", pA 
        write(file_log, *) " COEFFICIENT B                          [au] ", pB 
        write(file_log, *) " COEFFICIENT alpha                      [au] ", alpha 
        write(file_log, *) " COEFFICIENT beta                       [au] ", beta 
        write(file_log, *) " - "
        write(file_log, *) " - "
    else if(op_poten == 5) then 
        write(file_log, *) "================================================================="
        write(file_log, *) "POTENTIAL: CHEN & NOUS"
        write(file_log, *) "================================================================="
        write(file_log, *) " -------------------------------------------  ------------------ "
        write(file_log, *) " ELEMENT NUMBER                          [1] ", Z 
        write(file_log, *) " POLARIZABILITY                         [au] ", alphab
        write(file_log, *) " CUTOFF                                 [au] ", cutoff
        write(file_log, *) " COEFFICIENT alpha1                     [au] ", alpha1
        write(file_log, *) " COEFFICIENT alpha2                     [au] ", alpha2
        write(file_log, *) " COEFFICIENT alpha3                     [au] ", alpha3 
        write(file_log, *) " - "
        write(file_log, *) " - "
    end if 
    write(file_log, *) "================================================================="
    write(file_log, *) "UNIT"
    write(file_log, *) "================================================================="
    write(file_log, *) " -------------------------------------------  ------------------ "
    write(file_log, *) " DEFAULT                                               au        "
    if(op_ev     == 1) write(file_log, *) " ENERGY                                                eV        "
    if(op_degree == 1) write(file_log, *) " ANGULAR                                           degree        "
    if(op_aa     == 1) write(file_log, *) " CROSS SECTION                                        A^2        "
    write(file_log, *) " - "
    write(file_log, *) " - "
    write(file_log, *) "================================================================="
    write(file_log, *) "CALCULATION"
    write(file_log, *) "================================================================="
    write(file_log, *) " -------------------------------------------  ------------------ "
    write(file_log, *) " INNER REGION SIZE                      [au] ", Ra 
    write(file_log, *) " OUTER REGION SIZE                      [au] ", Rap 
    write(file_log, *) " MAXIUM OF ANGULAR MOMANTUM L            [1] ", L 
    write(file_log, *) " GRID NUMBER OF r COORDINATES            [1] ", N 
    write(file_log, *) " FLOQUET NUMBER                          [1] ", F 
    write(file_log, *) " GRID NUMBER OF Energy COORDINATES       [1] ", M  
    write(file_log, *) " - "
    write(file_log, *) " - "
end subroutine PROC_inform
! end information ----------------------------------
! coordination -------------------------------------
subroutine PORC_coord
    use fgsl,   only: fgsl_sf_legendre_Pl
    use linear, only: diag_sym_band
    real(dp), allocatable :: X(:, :)
    real(dp) :: tmp
    real(qp) :: sum 
    integer(i4) :: i, j, k 
    
    if(.not. allocated(coord_rho))    allocate(coord_rho   (0:N))
    if(.not. allocated(coord_weight)) allocate(coord_weight(0:N))
    if(.not. allocated(coord_dshape)) allocate(coord_dshape(0:N, 0:N))
    if(.not. allocated(X)) allocate(X(1:2, 1:N -1))
    X = 0.d0 
    do i = 1, N -2
        tmp = (dble(i)*dble(i +2))/(dble(2*i +1)*dble(2*i +3))
        tmp = tmp**0.5d0  
        X(1, i +1) = tmp 
    end do 
    call diag_sym_band(X)
    
    tmp = 2.d0/(dble(N +1)*dble(N))
    coord_rho   (0) = -1.d0
    coord_weight(0) = tmp 
    do i = 1, N -1
        coord_rho   (i) = X(1, i)
        coord_weight(i) = tmp/(fgsl_sf_legendre_Pl(N, coord_rho(i)))**2.d0 
    end do 
    coord_rho   (N) = 1.d0
    coord_weight(N) = tmp 
    dr_p_drho = Ra/(coord_rho(N) -coord_rho(0))
    if(allocated(X)) deallocate(X) 

    do i = 0, N 
        do j = 0, N 
            if(i == j) then 
                sum = 0.d0 
                do k = 0, N 
                    if(k /= i) then 
                        sum = sum +1.d0/(coord_rho(i) -coord_rho(k))
                    end if 
                end do 
            else if(i /= j) then 
                sum = 1.d0/(coord_rho(i) -coord_rho(j))
                do k = 0, N 
                    if(k /= i .and. k /= j) then 
                        sum = sum*(coord_rho(j) -coord_rho(k))/(coord_rho(i) -coord_rho(k))
                    end if 
                end do 
            end if 
            sum = sum*(coord_weight(j)/coord_weight(i))**0.5d0
            coord_dshape(i, j) = sum 
        end do 
    end do 
!     call TEST_poten
!     call TEST_coord
end subroutine PORC_coord
! end coordination ---------------------------------
! out ----------------------------------------------
subroutine PROC_hamiltonian_out
    if(allocated(coord_rho))    deallocate(coord_rho)
    if(allocated(coord_weight)) deallocate(coord_weight)
    if(allocated(coord_dshape)) deallocate(coord_dshape)
    close(file_log)
end subroutine PROC_hamiltonian_out
! end out ------------------------------------------
end module hamiltonian 
