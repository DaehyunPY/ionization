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
    coord_r = dr_p_drho*coord_rho(i) +ra/2.0_dp 
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
! nabla --------------------------------------------
! function Nabla_E(i, j)
!     use math_const, only: pi => math_pi 
!     integer(i4), intent(in) :: i, j 
!     real   (dp) :: Nabla_E, tmp, dbi, dbj 
!     dbi = dble(i)
!     dbj = dble(j)
!     if(i /= j) then 
!         tmp = ((dbi -dbj)**2.0_dp*sin(pi*(dbj +dbi))+pi*(-(dbi -dbj)**2.0_dp*(dbj +dbi)*cos(pi*(dbj +dbi)) &
!             -cos(pi*(dbi -dbj))*(dbi -dbj)*(dbj +dbi)**2.0_dp)+sin(pi*(dbi -dbj))*(dbj +dbi)**2.0_dp) &
!             /(pi*dE*(dbi -dbj)**2.0_dp*(dbj +dbi)**2.0_dp)
!     else if(i == j) then 
!         tmp = (sin(pi*(dbj +dbi))-pi*(dbj +dbi)*cos(pi*(dbj +dbi)))/(pi*dE*(dbj +dbi)**2.0_dp)
!     end if 
!     Nabla_E = -tmp 
! end function Nabla_E
! diff ---------------------------------------------
function diff_E(um2, um1, u0, up1, up2)
    use math_const, only: pi => math_pi 
    real(dp), intent(in) :: u0 
    real(dp), intent(in), optional :: um2, um1, up1, up2 
    real(dp) :: diff_E, tmp 
    tmp = 0.0_dp 
    if(present(um2) .and. present(up2)) then 
        tmp = +1.0_dp/12.0_dp * um2 &
              -8.0_dp/12.0_dp * um1 &
              +8.0_dp/12.0_dp * up1 &
              -1.0_dp/12.0_dp * up2 
    else if(present(um2) .and. (.not. present(up2))) then 
        tmp = +1.0_dp/2.0_dp  * um2 &
              -4.0_dp/2.0_dp  * um1 &
              +3.0_dp/2.0_dp  * u0 
    else if((.not. present(um2)) .and. present(up2)) then 
        tmp = -3.0_dp/2.0_dp  * u0  &
              +4.0_dp/2.0_dp  * up1 &
              -1.0_dp/2.0_dp  * up2 
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
    sum = 0.0_dp 
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
        tmp     = 7.75_dp*r**2.0_dp*exp(-r)
!         tmp     = -2.5_dp 
        Poten_r = tmp/Charge
    else if(ty == 1) then 
        Poten_r = Z/r 
    else if(ty == 2) then 
        stat    = -exp(-2.0_dp*z0*r)*(z1 +z2/r)
        pol     = -alphab/(2.0_dp*r**4.0_dp)*(1.0_dp -exp(-r))**6.0_dp 
        tmp     = stat +pol 
        Poten_r = tmp/(-1.0_dp)
    else if(ty == 3) then 
        stat    = -(Z/pH)*(exp(-r/pD)/r)*(1.0_dp +(pH -1.0_dp)*exp(-r*pH/pD))
        pol     = -alphab/(2.0_dp*(r**2.0_dp +pdelta**2.0_dp)**2.0_dp)
        tmp     = stat +pol 
        Poten_r = tmp/(-1.0_dp)
    else if(ty == 4) then 
        stat    = -Z*( &
                    +pA**2.0_dp/(4.0_dp*alpha**3.0_dp)*exp(-2.0_dp*alpha*r)*(alpha +1.0_dp/r) &
                    +4.0_dp*pA*pB/(alpha +beta)**3.0_dp*exp(-(alpha +beta)*r)*((alpha +beta)/2.0_dp +1.0_dp/r) & 
                    +pB**2.0_dp/(4.0_dp*beta**3.0_dp)*exp(-2.0_dp*beta*r)*(beta +1.0_dp/r))
        pol     = -alphab/(2.0_dp*(cutoff**2.0_dp +r**2.0_dp)**2.0_dp)
        tmp     = stat +pol 
        Poten_r = tmp/(-1.0_dp)
    else if(ty == 5) then 
        stat    = -(Z/r)*exp(-alpha1*r) -alpha2*exp(-alpha3*r)
        pol     = -alphab/(2.0_dp*r**4.0_dp) * (1.0_dp -exp(-(r/cutoff)**3.0_dp))**2.0_dp 
        tmp     = stat +pol 
        Poten_r = tmp/(-1.0_dp)
    end if 
end function Poten_r
! angular ------------------------------------------
function angular_r(l)
    integer(i4), intent(in) :: l 
    real   (dp) :: angular_r 
    angular_r = dble(l)*(dble(l) +1.0_dp)
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
    if(info /= 0) stop "Error #1: subroutine diag"

    X(:, :) = 0.0_dp 
    X(1, :) = W(:)
    if(allocated(W))    deallocate(W)
    if(allocated(Z))    deallocate(Z)
    if(allocated(work)) deallocate(work)
end subroutine diag










! ==================================================
! PROCESS
! ==================================================
! input --------------------------------------------
subroutine PROC_input
    use math_const, only: pi => math_pi
    use unit_const, only: other_e_eV, au_hartree
    character(60), parameter :: & 
        form_part  = '(4/, 3(45X, 1F15.8, /), /, 2(45X, 1F15.8, /), /)', &
        form_poten = '(4/, 1(45X, 1I15, /), 3(45X, 1F15.8, /))'
    character(30), parameter :: &
        form_p1    = '(20/)', &
        form_p2    = '( 2/, 3(45X, 1F15.8, /), 15/)', &
        form_p3    = '( 7/, 3(45X, 1F15.8, /), 10/)', &
        form_p4    = '(11/, 4(45X, 1F15.8, /),  5/)', &
        form_p5    = '(16/, 3(45X, 1F15.8, /),  1/)'
    character(60), parameter :: &
        form_cal   = '(4/, 1(45X, 1F15.8, /), 5(45X, 1I15, /), /)'
    character(90), parameter :: &
        form_opt   = '(6/, 3(45X, 6X, 1A1, /), /, 5(45X, 6X, 1A1, /), 3/, 2(45X, 6X, 1A1, /))'
    real     (dp), parameter :: eV_to_au = other_e_ev/au_hartree
    real     (dp) :: tmp1, tmp2 

    open(file_input, file = "input.d")
    read(file_input, form_part) Mass, Charge, Spin, tmp1, tmp2 
    if(Mass <  0.0_dp) stop "Input value is not valid: check 'PARTICLE - MASS'."
    if(Spin /= 0.5_dp) stop "Input value is not valid: check 'PARTICLE - SPIN'."
    if(tmp1 >= 0.0_dp .and. tmp2 >= 0.0_dp) stop "Input value is not valid: check 'PARTICLE - SCATTERING ENERGY'."
    if(tmp1 <  0.0_dp .and. tmp2 <  0.0_dp) stop "Input value is not valid: check 'PARTICLE - SCATTERING ENERGY'."
    if(tmp1 >= 0.0_dp) Scatt = tmp1 
    if(tmp2 >= 0.0_dp) Scatt = tmp2*eV_to_au
    read(file_input, form_poten) ty, Z, alphab, cutoff
    if(ty < 0 .or. ty > 6) stop "Input value is not valid: check 'POTENTIAL - TYPE'."
    if(ty == 0) then 
        read(file_input, form_p1)         
    else if(ty == 1) then 
        read(file_input, form_p1) 
    else if(ty == 2) then 
        read(file_input, form_p2) z0, z1, z2
    else if(ty == 3) then 
        read(file_input, form_p3) pD, pH, pdelta
        if(pD < 0.0_dp .and. pH < 0.0_dp) stop "Input value is not valid: check 'POTENTIAL - COEFFICIENT D, H'."
        if(pD < 0.0_dp)     pD     = pH/Z**0.4_dp 
        if(pH < 0.0_dp)     pH     = pD*Z**0.4_dp 
        if(pdelta < 0.0_dp) pdelta = (alphab/(2.0_dp*Z**(1.0_dp/3.0_dp)))**(0.25_dp)
    else if(ty == 4) then 
        read(file_input, form_p4) pA, pB, alpha, beta
    else if(ty == 5) then 
        read(file_input, form_p5) alpha1, alpha2, alpha3
    end if 
    read (file_input, form_cal) ra, N, M, L, pr, ptheta
    if(ra < 0.0_dp)  stop "Input value is not valid: check 'CALCULATION - BOUNDARY SIZE'."
    if(N  < 0)     stop "Input value is not valid: check 'CALCULATION - GRID NUMBER OF r COORDINATES'."
    if(M  < 0)     stop "Input value is not valid: check 'CALCULATION - GRID NUMBER OF Energy COORDINATES'."
    if(L  < 0)     stop "Input value is not valid: check 'CALCULATION - MAXIUM OF QUANTUM NUMBER L'."
    if(pr < 0)     stop "Input value is not valid: check 'CALCULATION - 3D PLOT NUMBER OF r COORDINATES'."
    if(ptheta < 0) stop "Input value is not valid: check 'CALCULATION - 3D PLOT NUMBER OF theta COORDINATES'."
    read (file_input, form_opt) & 
        op_ev, op_degree, op_aa, &
        op_poten, op_basis, op_dcs, op_inner, op_outer, &
        op_tcs, op_ps
    if(op_ev     /= "Y" .and. op_ev     /= "N") stop "Input value is not valid: check 'OPTION - USE ENERGY UNIT eV'."
    if(op_degree /= "Y" .and. op_degree /= "N") stop "Input value is not valid: check 'OPTION - USE ANGULAR UNIT degree'."
    if(op_aa     /= "Y" .and. op_aa     /= "N") stop "Input value is not valid: check 'OPTION - USE CROSS SECTION UNIT A^2'."
    if(op_poten  /= "Y" .and. op_poten  /= "N") stop "Input value is not valid: check 'OPTION - POTENTIAL FUNCTION'."
    if(op_basis  /= "Y" .and. op_basis  /= "N") stop "Input value is not valid: check 'OPTION - BASIS FUNCTION'."
    if(op_dcs    /= "Y" .and. op_dcs    /= "N") stop "Input value is not valid: check 'OPTION - CROSS SECTION FUNCTION'."
    if(op_inner  /= "Y" .and. op_inner  /= "N") stop "Input value is not valid: check 'OPTION - INNER REGION WAVE FUNCTION'."
    if(op_outer  /= "Y" .and. op_outer  /= "N") stop "Input value is not valid: check 'OPTION - OUTER REGION WAVE FUNCTION'."
    if(op_tcs    /= "Y" .and. op_tcs    /= "N") stop "Input value is not valid: check 'OPTION - KINETIC ENERGY vs CROSS SECTION'."
    if(op_ps     /= "Y" .and. op_ps     /= "N") stop "Input value is not valid: check 'OPTION - KINETIC ENERGY vs PHASE'."
    close(file_input) 
    open (file_log, file = "output/log.d")
    
    if(pr > N) pr = N 
    dtheta = pi/dble(ptheta)
    dE     = Scatt/dble(M) 

    if(.not. allocated(coord_rho))    allocate(coord_rho   (0:N))
    if(.not. allocated(coord_weight)) allocate(coord_weight(0:N))
    if(.not. allocated(coord_dshape)) allocate(coord_dshape(0:N, 0:N))
    if(op_basis == "N" .and. op_inner == "N") then 
        if(.not. allocated(H)) allocate(H(1:N, N:N))
    else if(op_basis == "Y" .or. op_inner == "Y") then 
        if(.not. allocated(H)) allocate(H(1:N, 1:N))
    end if 
    if(.not. allocated(E)) allocate(E(1:N))
    if(.not. allocated(R)) allocate(R(0:L))
    if(.not. allocated(K)) allocate(K(0:L))
    if(.not. allocated(S)) allocate(S(0:L))
    if(.not. allocated(A)) allocate(A(0:L))

    if(op_tcs == "N" .and. op_ps == "N") then 
        M  = 1
        dE = Scatt
    else if(op_tcs == "Y" .or. op_ps == "Y") then 
        op_basis = "N"
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
        write(file_log, *) "POTENTIAL: EXAMPLE 2.6.1"
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
    write(file_log, *) " BOUNDARY SIZE                          [au] ", ra 
    write(file_log, *) " GRID NUMBER OF r COORDINATES            [1] ", N 
    write(file_log, *) " MAXIUM OF QUANTUM NUMBER L              [1] ", L 
    write(file_log, *) " - "
    write(file_log, *) " - "
end subroutine PROC_inform
! end information ----------------------------------
! coordination -------------------------------------
subroutine PORC_coord
    use gsl_special, only: gsl_sf_legendre_Pl
    real   (dp) :: tmp
    real   (qp) :: sum 
    integer(i4) :: n, i, j, k 
    
    n = size(coord_rho(:)) -2 
    if(.not. allocated(X)) allocate(X(1:2, 1:n))
    X = 0.0_dp 
    do i = 1, n -1 
        tmp = (dble(i)*dble(i +2))/(dble(2*i +1)*dble(2*i +3))
        tmp = tmp**0.5_dp  
        X(1, i +1) = tmp 
    end do 
    call diag
    
    n   = size(coord_rho(:)) -1 
    tmp = 2.0_dp/(dble(n +1)*dble(n))
    coord_rho   (0) = -1.0_dp
    coord_weight(0) = tmp 
    do i = 1, n -1
        coord_rho   (i) = X(1, i)
        coord_weight(i) = tmp/(gsl_sf_legendre_Pl(n, coord_rho(i)))**2.0_dp 
    end do 
    coord_rho   (n) = 1.0_dp
    coord_weight(n) = tmp 
    dr_p_drho = ra/(coord_rho(n) -coord_rho(0))
    if(allocated(X)) deallocate(X) 

    do i = 0, n 
        do j = 0, n 
            if(i == j) then 
                sum = 0.0_dp 
                do k = 0, n 
                    if(k /= i) then 
                        sum = sum +1.0_dp/(coord_rho(i) -coord_rho(k))
                    end if 
                end do 
            else if(i /= j) then 
                sum = 1.0_dp/(coord_rho(i) -coord_rho(j))
                do k = 0, n 
                    if(k /= i .and. k /= j) then 
                        sum = sum*(coord_rho(j) -coord_rho(k))/(coord_rho(i) -coord_rho(k))
                    end if 
                end do 
            end if 
            sum = sum*(coord_weight(j)/coord_weight(i))**0.5_dp
            coord_dshape(i, j) = sum 
        end do 
    end do 
end subroutine PORC_coord
! end coordination ---------------------------------
! potential plot -----------------------------------
subroutine PROC_Poten_plot
    use gsl_special, only: gsl_sf_legendre_Pl
    integer  (i1), parameter :: file_poten = 101, file_coord = 102 
    character(30), parameter :: form_gen   = '(30ES25.10)'
    integer  (i4) :: i 

    open(file_poten, file = "output/poten.d")
    do i = 1, N
        write(file_poten, form_gen) coord_r(i), Poten_r(coord_r(i))
    end do
    close(file_poten)

    open(file_coord, file = "output/coord.d")
    write(file_coord, form_gen) coord_r(0_i4), coord_rho(0), coord_weight(0)
    do i = 1, N -1 
        write(file_coord, form_gen) coord_r(i), coord_rho(i), coord_weight(i), & 
            N*(coord_rho(i)*gsl_sf_legendre_Pl(N, coord_rho(i)) -gsl_sf_legendre_Pl(N -1_i4, coord_rho(i))) &
                /(coord_rho(i)**2.0_dp -1.0_dp)
    end do 
    write(file_coord, form_gen) coord_r(N), coord_rho(N), coord_weight(N)
    close(file_coord)
end subroutine PROC_Poten_plot
! end potential plot -------------------------------
! out ----------------------------------------------
subroutine PROC_out 
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
end subroutine PROC_out 
! end out ------------------------------------------
end module hamiltonian 
