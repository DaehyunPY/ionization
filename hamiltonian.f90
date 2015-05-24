module hamiltonian 
    use kind_type
    use global 
    implicit none
    real(dp), allocatable, private, protected :: X(:, :)
    real(dp), save, private, protected :: dtheta, dE 
    real(dp), save, private, protected :: &
        Z, alphab, cutoff, &   ! general 
        z0, z1, z2, &          ! II-1 GELTMAN
        pD, pH, pdelta, &      ! II-2 GREEN & YUKAWA
        pA, pB, alpha, beta, & ! II-3 KYLSTRA & BUCKINGHAM
        alpha1, alpha2, alpha3 ! III  CHEN & NOUS 
    integer(i1), save, private, protected :: ty 
contains


! ==================================================
! COORDNATION
! ==================================================
! r ------------------------------------------------
function coord_r(i) 
    integer(i4), intent(in) :: i
    real   (dp) :: coord_r
    coord_r = dr_p_drho*coord_rho(i) +Bound/2.d0 
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
        stop "Error #271: function diff_E"
    end if 
    diff_E = tmp/dE 
end function diff_E
! delta --------------------------------------------
function Delta_rho(i, j)
    use math_const, only: pi => math_pi 
    integer(i4), intent(in) :: i, j 
    real   (qp) :: sum 
    real   (dp) :: Delta_rho
    integer(i4) :: k, n  
    n   = dble(size(coord_rho(:))) -1 
    sum = 0.d0 
    do k = 0, n 
        sum = sum -coord_dshape(i, k)*coord_dshape(j, k)
    end do 
    Delta_rho = sum 
end function Delta_rho
! potential ----------------------------------------
function Poten_r(r)
    real(dp), intent(in) :: r
    real(dp) :: Poten_r, stat, pol, tmp 
    if(ty == 0) then 
!         tmp     = 7.75d0*r**2.d0*exp(-r) ! example 2.6.1
        tmp     = -2.5d0 ! example P G Burke 
        Poten_r = tmp/Charge
    else if(ty == 1) then 
        Poten_r = Z/r 
    else if(ty == 2) then 
        stat    = -exp(-2.d0*z0*r)*(z1 +z2/r)
        pol     = -alphab/(2.d0*r**4.d0)*(1.d0 -exp(-r))**6.d0 
        tmp     = stat +pol 
        Poten_r = tmp/(-1.d0)
    else if(ty == 3) then 
        stat    = -(Z/pH)*(exp(-r/pD)/r)*(1.d0 +(pH -1.d0)*exp(-r*pH/pD))
        pol     = -alphab/(2.d0*(r**2.d0 +pdelta**2.d0)**2.d0)
        tmp     = stat +pol 
        Poten_r = tmp/(-1.d0)
    else if(ty == 4) then 
        stat    = -Z*( &
                    +pA**2.d0/(4.d0*alpha**3.d0)*exp(-2.d0*alpha*r)*(alpha +1.d0/r) &
                    +4.d0*pA*pB/(alpha +beta)**3.d0*exp(-(alpha +beta)*r)*((alpha +beta)/2.d0 +1.d0/r) & 
                    +pB**2.d0/(4.d0*beta**3.d0)*exp(-2.d0*beta*r)*(beta +1.d0/r))
        pol     = -alphab/(2.d0*(cutoff**2.d0 +r**2.d0)**2.d0)
        tmp     = stat +pol 
        Poten_r = tmp/(-1.d0)
    else if(ty == 5) then 
        stat    = -(Z/r)*exp(-alpha1*r) -alpha2*exp(-alpha3*r)
        pol     = -alphab/(2.d0*r**4.d0) * (1.d0 -exp(-(r/cutoff)**3.d0))**2.d0 
        tmp     = stat +pol 
        Poten_r = tmp/(-1.d0)
    end if 
end function Poten_r
! angular ------------------------------------------
function angular_r(l)
    integer(i4), intent(in) :: l 
    real   (dp) :: angular_r 
    angular_r = dble(l)*(dble(l) +1.d0)
end function angular_r


! ==================================================
! CALCULATION 
! ==================================================
! dsbev diagonalization ----------------------------
subroutine diag
    character(1), parameter :: jobz = 'No vectors', uplo = 'Upper'
    integer (i8), parameter :: kd   = 1 
    integer (i8) :: n, ldab, ldz, info 
    real    (dp), allocatable :: W(:), Z(:, :), work(:)
    n    = size(X(1, :))
    ldab = size(X(:, 1))
    if(.not. allocated(W))    allocate(W(1:n))
    if(.not. allocated(Z))    allocate(Z(1, 1))
    if(.not. allocated(work)) allocate(work(1:3*n -2))
    ldz  = size(Z (:, 1)) 
    info = 0
    call DSBEV(jobz, uplo, n, kd, X, ldab, W, Z, ldz, work, info)
    if(info /= 0) stop "Error #153: subroutine diag"
    X(:, :) = 0.d0 
    X(1, :) = W(:)
    if(allocated(W))    deallocate(W)
    if(allocated(Z))    deallocate(Z)
    if(allocated(work)) deallocate(work)
end subroutine diag


! ==================================================
! TEST 
! ==================================================
! check potential ----------------------------------
! subroutine check_poten
!     integer  (i1), parameter :: file_poten = 101
!     character(30), parameter :: form_poten = '(30ES25.10)'
!     integer  (i4) :: i 

!     open(file_poten, file = "output/poten.d")
!     do i = 1, N
!         write(file_poten, form_poten) coord_r(i), Poten_r(coord_r(i))*Charge
!     end do
!     close(file_poten)
! end subroutine check_poten
! end check potential ------------------------------
! check coordination system ------------------------
! subroutine check_coord
!     use gsl_special, only: gsl_sf_legendre_Pl
!     integer  (i1), parameter :: file_coord = 101
!     character(30), parameter :: form_coord = '(30ES25.10)'
!     integer  (i4) :: i 

!     open(file_coord, file = "output/coord.d")
!     write(file_coord, form_coord) coord_r(0_i4), coord_rho(0), coord_weight(0)
!     do i = 1, N -1 
!         write(file_coord, form_coord) coord_r(i), coord_rho(i), coord_weight(i), & 
!             N*(coord_rho(i)*gsl_sf_legendre_Pl(N, coord_rho(i)) -gsl_sf_legendre_Pl(N -1_i4, coord_rho(i))) &
!                 /(coord_rho(i)**2.d0 -1.d0)
!     end do 
!     write(file_coord, form_coord) coord_r(N), coord_rho(N), coord_weight(N)
!     close(file_coord)
! end subroutine check_coord
! end check coordination system --------------------








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
        form_cal = '(4/, 1(45X, 1F15.8, /), 6(45X, 1I15, /), /)'
    character(120), parameter :: &
        form_opt   = & 
            '(6/, 3(45X, 6X, 1A1, /), /, 2(45X, 6X, 1A1, /), /, 3(45X, 6X, 1A1, /), 3/, 2(45X, 6X, 1A1, /))'
    real(dp),  parameter :: eV_to_au = other_e_ev/au_hartree
    real(dp) :: tmp1, tmp2 

    open(file_input, file = "input.d")
    ! laser field ----------------------------------
    read(file_input, form_laser) Amp, Freq
    if(.not. Amp  >  0.d0) stop "Error #051: check 'LASER FIELD - AMPLITUDE'"
    if(.not. Freq >  0.d0) stop "Error #052: check 'LASER FIELD - FREQUENCY'"
    ! particle -------------------------------------
    read(file_input, form_part ) Mass, Charge, Spin, tmp1, tmp2 
    if(.not. Mass >  0.d0)  stop "Error #053: check 'PARTICLE - MASS'"    
    if(.not. Spin == 0.5d0) stop "Error #054: check 'PARTICLE - SPIN'"
    if(tmp1 >= 0.d0 .and. tmp2 >= 0.d0) stop "Error #055: check 'PARTICLE - SCATTERING ENERGY'"
    if(tmp1 <  0.d0 .and. tmp2 <  0.d0) stop "Error #056: check 'PARTICLE - SCATTERING ENERGY'"
    if(tmp1 >= 0.d0) Scatt = tmp1 
    if(tmp2 >= 0.d0) Scatt = tmp2*eV_to_au
    ! potential ------------------------------------
    read(file_input, form_poten) ty, Z, alphab, cutoff
    if(.not. (ty >= 0 .and. ty <= 6)) stop "Error #051: check 'POTENTIAL - TYPE'"
    if(ty == 0) then 
        read(file_input, form_p1)         
    else if(ty == 1) then 
        read(file_input, form_p1) 
    else if(ty == 2) then 
        read(file_input, form_p2) z0, z1, z2
    else if(ty == 3) then 
        read(file_input, form_p3) pD, pH, pdelta
        if(pD < 0.d0 .and. pH < 0.d0) stop "Error #057: check 'POTENTIAL - COEFFICIENT D, H'"
        if(pD < 0.d0)     pD     = pH/Z**0.4d0 
        if(pH < 0.d0)     pH     = pD*Z**0.4d0 
        if(pdelta < 0.d0) pdelta = (alphab/(2.d0*Z**(1.d0/3.d0)))**(0.25d0)
    else if(ty == 4) then 
        read(file_input, form_p4) pA, pB, alpha, beta
    else if(ty == 5) then 
        read(file_input, form_p5) alpha1, alpha2, alpha3
    end if 
    ! calculation ----------------------------------
    read (file_input, form_cal) Bound, L, N, F, M, pr, ptheta
    if(.not. Bound  >  0.d0) stop "Error #058: check 'CALCULATION - BOUNDARY SIZE'"
    if(.not. L      >= 0)    stop "Error #059: check 'CALCULATION - MAXIUM OF ANGULAR MOMANTUM L'"
    if(.not. N      >  0)    stop "Error #060: check 'CALCULATION - GRID NUMBER OF r COORDINATES'"
    if(.not. F      >= 0)    stop "Error #061: check 'CALCULATION - FLOQUET NUMBER'"
    if(.not. M      >  0)    stop "Error #062: check 'CALCULATION - GRID NUMBER OF Energy COORDINATES'"
    if(.not. pr     >  0)    stop "Error #063: check 'CALCULATION - 3D PLOT NUMBER OF r COORDINATES'"
    if(.not. ptheta >  0)    stop "Error #064: check 'CALCULATION - 3D PLOT NUMBER OF theta COORDINATES'"
    if(pr > N) pr = N 
    dtheta = pi/dble(ptheta)
    dE     = Scatt/dble(M) 
    ! option ---------------------------------------
    read (file_input, form_opt) & 
        op_ev, op_degree, op_aa, & 
        op_basis, op_bound, & 
        op_dcs, op_inner, op_outer, &
        op_tcs, op_ps
    if(.not.(op_ev     == "Y" .or. op_ev     == "N")) stop "Error #065: check 'OPTION - USE ENERGY UNIT eV'"
    if(.not.(op_degree == "Y" .or. op_degree == "N")) stop "Error #066: check 'OPTION - USE ANGULAR UNIT degree'"
    if(.not.(op_aa     == "Y" .or. op_aa     == "N")) stop "Error #067: check 'OPTION - USE CROSS SECTION UNIT A^2'"
    if(.not.(op_basis  == "Y" .or. op_basis  == "N")) stop "Error #069: check 'OPTION - BASIS FUNCTION'"
    if(.not.(op_bound  == "Y" .or. op_bound  == "N")) stop "Error #415: check 'OPTION - BOUND STATES ENERGY & FUNCTION'"
    if(.not.(op_dcs    == "Y" .or. op_dcs    == "N")) stop "Error #070: check 'OPTION - CROSS SECTION FUNCTION'"
    if(.not.(op_inner  == "Y" .or. op_inner  == "N")) stop "Error #071: check 'OPTION - INNER REGION WAVE FUNCTION'"
    if(.not.(op_outer  == "Y" .or. op_outer  == "N")) stop "Error #072: check 'OPTION - OUTER REGION WAVE FUNCTION'"
    if(.not.(op_tcs    == "Y" .or. op_tcs    == "N")) stop "Error #073: check 'OPTION - KINETIC ENERGY vs CROSS SECTION'"
    if(.not.(op_ps     == "Y" .or. op_ps     == "N")) stop "Error #074: check 'OPTION - KINETIC ENERGY vs PHASE'"
    close(file_input) 
    open (file_log, file = "output/log.d")
    if(.not. (op_tcs == "Y" .or. op_ps == "Y")) then 
        M  = 1
        dE = Scatt
    else if(op_tcs == "Y" .or. op_ps == "Y") then 
        op_basis = "N"
        op_bound = "N"
        op_dcs   = "N"
        op_inner = "N" 
        op_outer = "N"
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

    if(ty == 0) then 
        write(file_log, *) "================================================================="
        write(file_log, *) "POTENTIAL: EXAMPLE"
        write(file_log, *) "================================================================="
        write(file_log, *) " -------------------------------------------  ------------------ "
        write(file_log, *) " - "
        write(file_log, *) " - "

    else if(ty == 1) then 
        write(file_log, *) "================================================================="
        write(file_log, *) "POTENTIAL: COULOMB"
        write(file_log, *) "================================================================="
        write(file_log, *) " -------------------------------------------  ------------------ "
        write(file_log, *) " ELEMENT NUMBER                          [1] ", Z 
        write(file_log, *) " - "
        write(file_log, *) " - "

    else if(ty == 2) then 
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

    else if(ty == 3) then 
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

    else if(ty == 4) then 
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

    else if(ty == 5) then 
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
    if(op_ev     == "Y") write(file_log, *) " ENERGY                                                eV        "
    if(op_degree == "Y") write(file_log, *) " ANGULAR                                           degree        "
    if(op_aa     == "Y") write(file_log, *) " CROSS SECTION                                        A^2        "
    write(file_log, *) " - "
    write(file_log, *) " - "

    write(file_log, *) "================================================================="
    write(file_log, *) "CALCULATION"
    write(file_log, *) "================================================================="
    write(file_log, *) " -------------------------------------------  ------------------ "
    write(file_log, *) " BOUNDARY SIZE                          [au] ", Bound 
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
    use gsl_special, only: gsl_sf_legendre_Pl
    real   (dp) :: tmp
    real   (qp) :: sum 
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
    call diag
    
    tmp = 2.d0/(dble(N +1)*dble(N))
    coord_rho   (0) = -1.d0
    coord_weight(0) = tmp 
    do i = 1, N -1
        coord_rho   (i) = X(1, i)
        coord_weight(i) = tmp/(gsl_sf_legendre_Pl(N, coord_rho(i)))**2.d0 
    end do 
    coord_rho   (N) = 1.d0
    coord_weight(N) = tmp 
    dr_p_drho = Bound/(coord_rho(N) -coord_rho(0))
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
!     call check_poten ! for test 
!     call check_coord ! for test 
end subroutine PORC_coord
! end coordination ---------------------------------
! out ----------------------------------------------
subroutine PROC_prog_out 
    if(allocated(coord_rho))    deallocate(coord_rho)
    if(allocated(coord_weight)) deallocate(coord_weight)
    if(allocated(coord_dshape)) deallocate(coord_dshape)
    if(allocated(H)) deallocate(H)
    if(allocated(E)) deallocate(E)
    if(allocated(R)) deallocate(R)
    if(allocated(K)) deallocate(K)
    if(allocated(S)) deallocate(S)
    if(allocated(A)) deallocate(A)
    close(file_log)
end subroutine PROC_prog_out 
! end out ------------------------------------------
end module hamiltonian 
