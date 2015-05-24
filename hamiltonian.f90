module hamiltonian
    use kind_const
    use global
    implicit none
    ! coordnation
    real(dp), save, protected :: dr_pdrho
    real(dp), save, protected, private :: dtheta, aAmp, bAmp, dAmp
    real(dp), save, allocatable, protected :: weight_grid(:)
    real(dp), save, allocatable, protected, private :: rho_grid(:), dshape_grid(:, :)
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
function r_grid(i)
    integer(i4), intent(in) :: i
    real(dp) :: r_grid
    r_grid = dr_pdrho*rho_grid(i) +Ba/2.d0
end function r_grid
! theta --------------------------------------------
function theta_grid(i)
    integer(i4), intent(in) :: i
    real(dp) :: theta_grid
    theta_grid = dtheta*dble(i)
end function theta_grid
! amplitude ----------------------------------------
function Amp_grid(i)
    integer(i4), intent(in) :: i
    real(dp) :: Amp_grid
    if(i == 0) then
        Amp_grid = 0.d0
    else if(i > 0) then
        Amp_grid = dAmp*dble(i) +bAmp
    end if
end function Amp_grid
! end coordination ---------------------------------


! ==================================================
! OPERATOR
! ==================================================
! delta --------------------------------------------
function Delta_grid(i, j)
    use math_const, only: pi => math_pi
    integer(i4), intent(in) :: i, j
    real   (qp) :: sum
    real   (dp) :: Delta_grid
    integer(i4) :: k, n
    n   = dble(size(rho_grid(:))) -1
    sum = 0.d0
    do k = 0, n
        sum = sum -dshape_grid(i, k)*dshape_grid(j, k)
    end do
    Delta_grid = sum
end function Delta_grid
! potential ----------------------------------------
function poten_r(r)
    real(dp), intent(in) :: r
    real(dp) :: poten_r, stat, pol
    if(op_poten == 0) then
!         poten_r = 7.75d0*r**2.d0*exp(-r) ! example 2.6.1
        poten_r = -2.5d0 ! example P G Burke
    else if(op_poten == 1) then
        poten_r = -Z/r
    else if(op_poten == 2) then
        stat    = -exp(-2.d0*z0*r)*(z1 +z2/r)
        pol     = -alphab/(2.d0*r**4.d0)*(1.d0 -exp(-r))**6.d0
        poten_r = stat +pol
    else if(op_poten == 3) then
        stat    = -(Z/pH)*(exp(-r/pD)/r)*(1.d0 +(pH -1.d0)*exp(-r*pH/pD))
        pol     = -alphab/(2.d0*(r**2.d0 +pdelta**2.d0)**2.d0)
        poten_r = stat +pol
    else if(op_poten == 4) then
        stat    = -Z*( &
                    +pA**2.d0/(4.d0*alpha**3.d0)*exp(-2.d0*alpha*r)*(alpha +1.d0/r) &
                    +4.d0*pA*pB/(alpha +beta)**3.d0*exp(-(alpha +beta)*r)*((alpha +beta)/2.d0 +1.d0/r) &
                    +pB**2.d0/(4.d0*beta**3.d0)*exp(-2.d0*beta*r)*(beta +1.d0/r))
        pol     = -alphab/(2.d0*(cutoff**2.d0 +r**2.d0)**2.d0)
        poten_r = stat +pol
    else if(op_poten == 5) then
        stat    = -(Z/r)*exp(-alpha1*r) -alpha2*exp(-alpha3*r)
        pol     = -alphab/(2.d0*r**4.d0) * (1.d0 -exp(-(r/cutoff)**3.d0))**2.d0
        poten_r = stat +pol
    end if
end function poten_r
! dipole -------------------------------------------
function dipole_r(r, l)
    real(dp), intent(in) :: r
    integer(i4), intent(in) :: l
    real(dp) :: dipole_r
    dipole_r = r
end function dipole_r
! end operator -------------------------------------










! ==================================================
! PROCESS
! ==================================================
! input --------------------------------------------
subroutine PROC_input
    use math_const, only: pi => math_pi
    use unit_const, only: other_e_eV, au_hartree
    character(60),  parameter :: &
        form_laser = '(4/, 3(45X, 1F15.8, /), /)', &
        form_poten = '(4/, 1(45X, 1I15, /), 3(45X, 1F15.8, /))'
    character(30),  parameter :: &
        form_p1 = '(20/)', &
        form_p2 = '( 2/, 3(45X, 1F15.8, /), 15/)', &
        form_p3 = '( 7/, 3(45X, 1F15.8, /), 10/)', &
        form_p4 = '(11/, 4(45X, 1F15.8, /),  5/)', &
        form_p5 = '(16/, 3(45X, 1F15.8, /),  1/)'
    character(60),  parameter :: &
        form_cal = '(4/, 3(45X, 1F15.8, /), 5(45X, 1I15, /), /)', &
        form_opt = '(6/, 3(45X, 6X, 1I15, /), /, /)'
    real(dp),  parameter :: eV_to_au = other_e_ev/au_hartree

    open(file_input, file = "input.data")

    ! laser field
    read(file_input, form_laser) Freq, aAmp, bAmp
    if(.not. Freq > 0.d0) stop "Error hamiltonian #139: Check 'FREQUENCY'. "
    if(.not. aAmp > 0.d0) stop "Error hamiltonian #140: Check 'MAXIMUM AMPLITUDE'. "
    if(.not. (bAmp >= 0.d0 .and. bAmp < aAmp)) stop "Error hamiltonian #140: Check 'MINIMUM AMPLITUDE'. "

    ! potential
    read(file_input, form_poten) op_poten, Z, alphab, cutoff
    if(.not. (op_poten >= 0 .and. op_poten <= 6)) stop "Error hamiltonian #301: Check 'POTENTIAL TYPE'. "
    if(op_poten == 0) then
        read(file_input, form_p1)
    else if(op_poten == 1) then
        read(file_input, form_p1)
    else if(op_poten == 2) then
        read(file_input, form_p2) z0, z1, z2
    else if(op_poten == 3) then
        read(file_input, form_p3) pD, pH, pdelta
        if(pD < 0.d0 .and. pH < 0.d0) stop "Error hamiltonian #302: Check '[TYPE 3] COEFFICIENT D, H'. "
        if(pD < 0.d0)     pD     = pH/Z**0.4d0
        if(pH < 0.d0)     pH     = pD*Z**0.4d0
        if(pdelta < 0.d0) pdelta = (alphab/(2.d0*Z**(1.d0/3.d0)))**(0.25d0)
    else if(op_poten == 4) then
        read(file_input, form_p4) pA, pB, alpha, beta
    else if(op_poten == 5) then
        read(file_input, form_p5) alpha1, alpha2, alpha3
    end if

    ! calculation
    read (file_input, form_cal) Ba, Bap, Bbeta, L, N, F, M, p3d
    if(.not. Ba > 0.d0)    stop "Error hamiltonian #166: Check 'INNER REGION SIZE'. "
    if(.not. Bap >= Ba)    stop "Error hamiltonian #167: Check 'OUTER REGION SIZE'. "
    if(.not. Bbeta > 0.d0) stop "Error hamiltonian #168: Check 'PROPAGATE PARAMETER'. "
    if(.not. L >= 0)  stop "Error hamiltonian #403: Check 'MAXIMUM OF ANGULAR MOMANTUM L'. "
    if(.not. N >  0)  stop "Error hamiltonian #404: Check 'GRID NUMBER OF r COORDINATES'. "
    if(.not. F >= 0)  stop "Error hamiltonian #405: Check 'FLOQUET BLOCK NUMBER'. "
    if(.not. M >  0)  stop "Error hamiltonian #406: Check 'GRID NUMBER OF Amplitude COORDINATES'. "
    if(.not. p3d > 0) stop "Error hamiltonian #407: Check 'GRID NUMBER OF 3D PLOT'. "
    dtheta = pi/dble(p3d)
    dAmp = (aAmp -bAmp)/dble(M)

    ! option
    read (file_input, form_opt) op_ev, op_degree, op_aa
    if(.not.(op_ev     == 1 .or. op_ev     == 0)) stop "Error hamiltonian #501: Check '[UNIT] USE ENERGY UNIT eV'. "
    if(.not.(op_degree == 1 .or. op_degree == 0)) stop "Error hamiltonian #502: Check '[UNIT] USE ANGULAR UNIT degree'. "
    if(.not.(op_aa     == 1 .or. op_aa     == 0)) stop "Error hamiltonian #503: Check '[UNIT] USE CROSS SECTION UNIT A^2'. "
    close(file_input)
    open (file_log, file = "output/log.data")
end subroutine PROC_input
! information --------------------------------------
subroutine PROC_inform
    use unit_const, only: other_e_eV, au_hartree
    character(30), parameter :: form_out = '(1A15, 1ES15.3)'
    real     (dp), parameter :: au_to_eV = au_hartree/other_e_ev

    write(file_log, *) "================================================================="
    write(file_log, *) "LASER FIELD SCATTERING"
    write(file_log, *) "================================================================="
    write(file_log, *) " -------------------------------------------  ------------------ "
    write(file_log, *) " FREQUENCY                              [au] ", Freq
    write(file_log, *) " MAXIMUM AMPLITUDE                      [au] ", aAmp
    write(file_log, *) " MINIMUM AMPLITUDE                      [au] ", bAmp
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
    write(file_log, *) " INNER REGION SIZE                      [au] ", Ba
    write(file_log, *) " OUTER REGION SIZE                      [au] ", Bap
    write(file_log, *) " MAXIMUM OF ANGULAR MOMANTUM L           [1] ", L
    write(file_log, *) " GRID NUMBER OF r COORDINATES            [1] ", N
    write(file_log, *) " FLOQUET BLOCK NUMBER                    [1] ", F
    write(file_log, *) " GRID NUMBER OF Amplitude COORDINATES    [1] ", M
    write(file_log, *) " GRID NUMBER OF 3D PLOT                  [1] ", p3d
    write(file_log, *) " - "
    write(file_log, *) " - "
end subroutine PROC_inform
! coordination -------------------------------------
subroutine PORC_coord
    use fgsl,     only: fgsl_sf_legendre_Pl
    use mylinear, only: diag_sym_band
    real(dp), allocatable :: X(:, :)
    real(dp) :: tmp
    real(qp) :: sum
    integer(i4) :: i, j, k

    if(.not. allocated(rho_grid))    allocate(rho_grid   (0:N))
    if(.not. allocated(weight_grid)) allocate(weight_grid(0:N))
    if(.not. allocated(dshape_grid)) allocate(dshape_grid(0:N, 0:N))
    if(.not. allocated(X)) allocate(X(1:2, 1:N -1))
    X = 0.d0
    do i = 1, N -2
        tmp = (dble(i)*dble(i +2))/(dble(2*i +1)*dble(2*i +3))
        tmp = tmp**0.5d0
        X(1, i +1) = tmp
    end do
    call diag_sym_band(X)

    tmp = 2.d0/(dble(N +1)*dble(N))
    rho_grid   (0) = -1.d0
    weight_grid(0) = tmp
    do i = 1, N -1
        rho_grid   (i) = X(1, i)
        weight_grid(i) = tmp/(fgsl_sf_legendre_Pl(N, rho_grid(i)))**2.d0
    end do
    rho_grid   (N) = 1.d0
    weight_grid(N) = tmp
    dr_pdrho = Ba/(rho_grid(N) -rho_grid(0))
    if(allocated(X)) deallocate(X)

    do i = 0, N
        do j = 0, N
            if(i == j) then
                sum = 0.d0
                do k = 0, N
                    if(k /= i) then
                        sum = sum +1.d0/(rho_grid(i) -rho_grid(k))
                    end if
                end do
            else if(i /= j) then
                sum = 1.d0/(rho_grid(i) -rho_grid(j))
                do k = 0, N
                    if(k /= i .and. k /= j) then
                        sum = sum*(rho_grid(j) -rho_grid(k))/(rho_grid(i) -rho_grid(k))
                    end if
                end do
            end if
            sum = sum*(weight_grid(j)/weight_grid(i))**0.5d0
            dshape_grid(i, j) = sum
        end do
    end do
end subroutine PORC_coord
! out ----------------------------------------------
subroutine PROC_hamiltonian_out
    if(allocated(rho_grid))    deallocate(rho_grid)
    if(allocated(weight_grid)) deallocate(weight_grid)
    if(allocated(dshape_grid)) deallocate(dshape_grid)
    close(file_log)
end subroutine PROC_hamiltonian_out
! end process --------------------------------------
end module hamiltonian
