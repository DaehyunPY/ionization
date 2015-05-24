module boundary
    use kind_type
    use global
    implicit none
    complex(dp), save, allocatable, protected, private :: k(:), R(:, :), B(:, :), Scatt(:, :)
    integer(i4), parameter, private :: &
        file_cmplx = 201, file_ion = 202, file_branch1 = 203, file_branch2 = 204, &
        addblock = 10, addcount = 3
    real(dp), save, protected, private :: A0, alpha0, Ep
contains


! ==================================================
! FUNCTIONS
! ==================================================
! single R matrix ----------------------------------
function matching_R1(energy)
    use basis, only: H1, E1
    integer(i4) :: i
    real(dp) :: matching_R1
    real(dp), intent(in) :: energy
    real(qp) :: sum

    sum = 0.d0
    do i = 1, N
        sum = sum +H1(N, i)*H1(N, i)/(E1(i) -energy)
    end do
    matching_R1 = sum/(2.d0*Ba)
end function matching_R1
! matching parameter d1 ----------------------------
function matching_d1(energy, l)
    use fgsl, only: fgsl_sf_bessel_ksl_scaled
    integer(i4), intent(in) :: l
    real(dp) :: ka, sb_ka, diff_ka, tmp1, tmp2, matching_d1
    real(dp), intent(in) :: energy

    ka = (2.d0*(-energy))**0.5d0*Ba
    sb_ka = fgsl_sf_bessel_ksl_scaled(l, ka)/exp(ka)*ka
    if(.not. l == 0) then
        tmp1 = dble(l)*fgsl_sf_bessel_ksl_scaled(l, ka)/exp(ka)
        tmp2 = dble(l +1)*fgsl_sf_bessel_ksl_scaled(l +1_i4, ka)/exp(ka)
    else if(l == 0) then
        tmp1 = 0.d0
        tmp2 = dble(l +1)*fgsl_sf_bessel_ksl_scaled(l +1_i4, ka)/exp(ka)
    end if
    diff_ka = -ka**2.d0/dble(2*l +1)*(tmp1 +tmp2)
    matching_d1 = sb_ka -matching_R1(energy)*(sb_ka +diff_ka)

    if(1.d0 +abs(matching_d1) == 1.d0) then
        matching_d1 = 0.d0
    end if
end function matching_d1
! floquet R matrix ---------------------------------
function matching_RF(f1, f2, energy)
    use basis, only: HF, EF
    complex(dp) :: matching_RF
    complex(dp), intent(in) :: energy
    complex(qp) :: sum
    integer(i4) :: i
    integer(i4), intent(in) :: f1, f2

    sum = 0.d0
    do i = 1, N*(2*F +1)
        sum = sum +HF(N, f1, i)*HF(N, f2, i)/(EF(i) -energy)
    end do
    matching_RF = sum/(2.d0*Ba)
end function matching_RF
! transfer pphi function ---------------------------
function transfer_pphi(l, x)
    integer(i4) :: m, m1, m2, dm, j, count
    integer(i4), intent(in) :: l
    real(dp), intent(in) :: x
    real(dp) :: transfer_pphi, inp1, inp2, tmp1, tmp2
    real(qp) :: sum, dsum1 ,dsum2

    inp1 = A0**2.d0/Freq/8.d0
    inp2 = A0*x
    if(mod(l, 2) == 0) then
        m1 = min(0, l/2)
        m2 = max(0, l/2)
    else if(mod(l, 2) == 1 .or. mod(l, 2) == -1) then
        if(l > 0) then
            m1 = 0
            m2 = (l +1)/2
        else if(l < 0) then
            m1 = (l -1)/2
            m2 = 0
        end if
    end if

    sum  = 0.d0

        dsum2 = 0.d0
        do m = m1, m2
            tmp2 = bessel_jn(m, inp1)*bessel_jn(l -2*m, inp2)
            if(m == m1) tmp1 = tmp2
            dsum2 = dsum2 +tmp2
        end do
        sum = sum +dsum2

    count = 0
    do j = 0, maxloop
        dsum1 = 0.d0
        dsum2 = 0.d0
        do dm = 1, addblock
            m = m1 -j*addblock -dm
            tmp1 = bessel_jn(m, inp1)*bessel_jn(l -2*m, inp2)
            dsum1 = dsum1 +tmp1

            m = m2 +j*addblock +dm
            tmp2 = bessel_jn(m, inp1)*bessel_jn(l -2*m, inp2)
            dsum2 = dsum2 +tmp2
        end do

        tmp1 = sum +dsum1
        tmp2 = sum +dsum2
        if(real(tmp1, kind = dp) == real(sum, kind = dp) &
            .and. real(tmp2, kind = dp) == real(sum, kind = dp)) then
                count = count +1
                sum = sum +dsum1 +dsum2
                if(count == addcount) exit
                cycle
        end if

        count = 0
        sum = sum +dsum1 +dsum2
        if(j == maxloop) print *, "Warning boundary #122: Do-loop arrive at 'maxloop'. "
    end do
    transfer_pphi = sum
end function transfer_pphi
! transfer phi function ----------------------------
function transfer_phi(m, n, x)
    use math_const, only: i => math_i
    use mod_zbes, only: bessel1
    complex(dp) :: transfer_phi, cml, inp1, inp2, sign, tmp0, tmp1, tmp2
    complex(qp) :: sum, dsum1, dsum2
    integer(i4) :: l, l1, l2, dl, j, info, count
    integer(i4), intent(in) :: m, n
    real(dp), intent(in) :: x

    inp1 = exp(i*k(m)*x)
    inp2 = k(m)*alpha0
    l1 = 0
    l2 = n -m
    if(l2 < 0) then
        l1 = n -m
        l2 = 0
    end if

    sum = 0.d0

        dsum2 = 0.d0
        do l = l1, l2
            if(l >= 0) then
                cml = dble(l)
                sign = 1.d0
                call bessel1(cml, inp2, tmp0, info)
                tmp0 = sign*tmp0
            else
                cml = dble(-l)
                if(mod(l, 2) == 0) then
                    sign = 1.d0
                else if(mod(l, 2) == 1 .or. mod(l, 2) == -1) then
                    sign = -1.d0
                end if
                call bessel1(cml, inp2, tmp0, info)
                tmp0 = sign*tmp0
            end if
            if(info == 20 .or. info == 10) then
                tmp0 = 0.d0
            end if
            if(mod(l, 4) == 0) then
                sign = 1.d0
            else if(mod(l, 4) == 1 .or. mod(l, 4) == -3) then
                sign = i
            else if(mod(l, 4) == 2 .or. mod(l, 4) == -2) then
                sign = -1.d0
            else if(mod(l, 4) == 3 .or. mod(l, 4) == -1) then
                sign = -i
            end if
            tmp2 = sign*tmp0*transfer_pphi(l -n +m, x)
            if(l == l1) tmp1 = tmp2
            dsum2 = dsum2 +tmp2
        end do
        sum = sum +dsum2

    count = 0
    do j = 0, maxloop
        dsum1 = 0.d0
        dsum2 = 0.d0
        do dl = 1, addblock
            l = l1 -j*addblock -dl
            if(l >= 0) then
                cml = dble(l)
                sign = 1.d0
                call bessel1(cml, inp2, tmp0, info)
                tmp0 = sign*tmp0
            else
                cml = dble(-l)
                if(mod(l, 2) == 0) then
                    sign = 1.d0
                else if(mod(l, 2) == 1 .or. mod(l, 2) == -1) then
                    sign = -1.d0
                end if
                call bessel1(cml, inp2, tmp0, info)
                tmp0 = sign*tmp0
            end if
            if(info == 20 .or. info == 10) then
                tmp0 = 0.d0
            end if
            if(mod(l, 4) == 0) then
                sign = 1.d0
            else if(mod(l, 4) == 1 .or. mod(l, 4) == -3) then
                sign = i
            else if(mod(l, 4) == 2 .or. mod(l, 4) == -2) then
                sign = -1.d0
            else if(mod(l, 4) == 3 .or. mod(l, 4) == -1) then
                sign = -i
            end if
            tmp1 = sign*tmp0*transfer_pphi(l -n +m, x)
            dsum1 = dsum1 +tmp1

            l = l2 +j*addblock +dl
            if(l >= 0) then
                cml = dble(l)
                sign = 1.d0
                call bessel1(cml, inp2, tmp0, info)
                tmp0 = sign*tmp0
            else
                cml = dble(-l)
                if(mod(l, 2) == 0) then
                    sign = 1.d0
                else if(mod(l, 2) == 1 .or. mod(l, 2) == -1) then
                    sign = -1.d0
                end if
                call bessel1(cml, inp2, tmp0, info)
                tmp0 = sign*tmp0
            end if
            if(info == 20 .or. info == 10) then
                tmp0 = 0.d0
            end if
            if(mod(l, 4) == 0) then
                sign = 1.d0
            else if(mod(l, 4) == 1 .or. mod(l, 4) == -3) then
                sign = i
            else if(mod(l, 4) == 2 .or. mod(l, 4) == -2) then
                sign = -1.d0
            else if(mod(l, 4) == 3 .or. mod(l, 4) == -1) then
                sign = -i
            end if
            tmp2 = sign*tmp0*transfer_pphi(l -n +m, x)
            dsum2 = dsum2 +tmp2
        end do

        tmp1 = sum +dsum1
        tmp2 = sum +dsum2
        if(cmplx(tmp1, kind = dp) == cmplx(sum, kind = dp) &
            .and. cmplx(tmp2, kind = dp) == cmplx(sum, kind = dp)) then
                count = count +1
                sum = sum +dsum1 +dsum2
                if(count == addcount) exit
                cycle
        end if

        count = 0
        sum = sum +dsum1 +dsum2
        if(j == maxloop) print *, "Warning boundary #239: Do-loop arrive at 'maxloop'. "
    end do
    transfer_phi = inp1*sum
end function transfer_phi
! differential transfer phi funciton ---------------
function transfer_dphi(m, n, x)
    use math_const, only: i => math_i
    integer(i4), intent(in) :: m, n
    real(dp), intent(in) :: x
    complex(dp) :: transfer_dphi
    complex(qp) :: sum

    sum = 0.d0
    sum = sum +i*k(m)*transfer_phi(m, n, x)
    sum = sum +A0*0.5d0*transfer_phi(m, n +1, x)
    sum = sum -A0*0.5d0*transfer_phi(m, n -1, x)
    transfer_dphi = sum
end function transfer_dphi
! matching parameter dF ----------------------------
function matching_dF(m, energy)
    use linear, only: det_cmplx
    complex(dp) :: matching_dF
    complex(dp), intent(in) :: energy
    integer(i4), intent(in) :: m

    if(allocated(k)) deallocate(k)
    if(allocated(R)) deallocate(R)
    if(allocated(B)) deallocate(B)
    allocate(k(-F:F))
    allocate(R(-F:F, -F:F))
    allocate(B(-F:F, -F:F))
    call ssub_channel_k(m, energy)
    call ssub_lngth_R(energy)
    call ssub_matching_B
    matching_dF = det_cmplx(B)
    if(allocated(k)) deallocate(k)
    if(allocated(R)) deallocate(R)
    if(allocated(B)) deallocate(B)
end function matching_dF
! end functions ------------------------------------










! ==================================================
! SUB-SUB-CALCULATION
! ==================================================
! channel p ----------------------------------------
subroutine ssub_channel_k(m, energy)
    use math_const, only: i => math_i, pi => math_pi
    use hamiltonian, only: Amp_grid
    complex(dp) :: tmp
    complex(dp), intent(in) :: energy
    integer(i4) :: j1
    integer(i4), intent(in) :: m

    A0 = -Amp_grid(m)/Freq
    alpha0 = Amp_grid(m)/Freq**2.d0
    Ep = A0**2.d0/4.d0

    k(:) = 0.d0
    do j1 = -F, F
        tmp = 2.d0*(energy -Ep +dble(j1)*Freq)
!         k(j1) = tmp**0.5d0
        k(j1) = -tmp**0.5d0
    end do
end subroutine ssub_channel_k
! length gauge R -----------------------------------
subroutine ssub_lngth_R(energy)
    complex(dp), intent(in) :: energy
    integer(i4) :: i, j

    R(:, :) = 0.d0
    do j = -F, F
        do i = -F, F
            R(i, j) = matching_RF(i, j, energy)
        end do
    end do
end subroutine ssub_lngth_R
! matching boundary condition -----------------------
subroutine ssub_matching_B
    integer(i4) :: n1, n2, m
    complex(dp), allocatable :: dB(:, :)
    complex(qp) :: sum

    if(allocated(dB)) deallocate(dB)
    allocate(dB(-F:F, -F:F))
    dB(:, :) = 0.d0
    do m = -F, F
        do n2 = -F, F
            dB(m, n2) = transfer_dphi(m, n2, Ba)
        end do
    end do
    B(:, :) = 0.d0
    do m = -F, F
        do n1 = -F, F
!             sum = 0.d0
!             do n2 = -F, F
!                 sum = sum -Ba*R(n1, n2)*dB(m, n2)
!             end do
!             B(n1, m) = transfer_phi(m, n1, Ba) +sum
            B(n1, m) = transfer_phi(m, n1, Ba) -Ba*R(n1, n1)*dB(m, n1)
        end do
    end do
    if(allocated(dB)) deallocate(dB)
end subroutine ssub_matching_B
! solve a eigen-vector of zero eigen-value ---------
subroutine ssub_eigen_vector(A, x)
    use linear, only: solve_cmplx_vec, det_cmplx
    character(30), parameter :: form_out = "(1I15, 9999ES15.3)"
    complex(dp), allocatable :: work(:, :)
    complex(dp), intent(in) :: A(1:, 1:)
    complex(dp), intent(out) :: x(1:)
    integer(i8) :: i1, j1, i2, j2, n, n0, n1, k, info
    real(qp) :: norm

    n = size(A(1:, 1))
    n0 = int(n/2) +1
    n1 = n0
    if(allocated(work)) deallocate(work)
    allocate(work(1:(n -1), 1:(n -1)))
    if(.not. (size(A(1, 1:)) == n .and. size(x(1:)) == n)) then
        stop "Error boundary #208: A and x have different dimention. "
    end if

    do k = 1, 2
        x(:) = 0.d0
        work(:, :) = 0.d0
        do j1 = 1, n
            if(j1 < n1) then
                j2 = j1
            else if(j1 == n1) then
                do i1 = 1, n
                    if(i1 < n1) then
                        i2 = i1
                    else if(i1 == n1) then
                        cycle
                    else if(i1 > n1) then
                        i2 = i1 -1
                    end if
                    x(i2) = -A(i1, j1)
                end do
                cycle
            else if(j1 > n1) then
                j2 = j1 -1
            end if
            do i1 = 1, n
                if(i1 < n1) then
                    i2 = i1
                else if(i1 == n1) then
                    cycle
                else if(i1 > n1) then
                    i2 = i1 -1
                end if
                work(i2, j2) = A(i1, j1)
            end do
        end do
        x(n) = 1.d0

        info = 0
        call solve_cmplx_vec(work, x(1:(n -1)), info)

        if(info == 0) then
            exit
        else if(info > 0) then
            if(k == 1) then
                if(info < n1) then
                    n1 = info
                else
                    n1 = info +1
                end if
            else
                stop "Error boundary #257: matrix-A has two or more zero eigen-value. "
            end if
        end if
    end do
    if(k == 2) print *, "Warning boundary #261: You can ignore the warning message. "
    if(allocated(work)) deallocate(work)

    do i1 = n, 1, -1
        if(i1 <= n1) exit
        x(i1) = x(i1 -1)
    end do
    x(n1) = 1.d0
!     if(n1 /= n0) then
!         x(:) = x(:)*abs(x(n0))/x(n0)
!     end if
    norm = 0.d0
    do i1 = 1, n
        norm = norm +abs(x(i1))**2.d0
    end do
    x(:) = x(:)/norm**0.5d0
end subroutine ssub_eigen_vector
! end sub-sub-calculation --------------------------










! ==================================================
! SUB-CALCULATION
! ==================================================
! find initial bound state -------------------------
subroutine SUB_initial_bound(l, energy)
    use basis, only: H1, E1
    integer(i4), intent(in) :: l
    real(dp), intent(out) :: energy
    real(dp) :: x1, x2, dx, dxold, f, fx1, fx2, df, sign, tmp
    integer(i4) :: i

    x1 = E1(1)
    x2 = E1(2)
    if(real(1.d0 +abs(H1(N, 1)), kind = sp) == real(1.d0, kind = sp)) then
        energy = x1

    else
        energy = 0.d0
        if(.not. x1 < 0.d0) stop "Error boundary #700: There is no bound state. "
        if(x2 > 0.d0) x2 = 0.d0
        dxold = abs(x2 -x1)

        dx = abs(x1*tiny3)
        if(.not. dx < 0.5d0*dxold) stop "Error boundary #705: 'tiny3' is too big. "
        x1 = x1 +dx
        x2 = x2 -dx

        fx1 = matching_d1(x1, l)
        fx2 = matching_d1(x2, l)
        sign = fx1*fx2
        if(sign > 0.d0) then
            stop "Error boundary #713: Bound state must be bracketed. "
        else if(fx1 == 0.d0) then
            energy = x1
            return
        else if(fx2 == 0.d0) then
            energy = x2
            return
        end if

        tmp = 1.d0
        energy = (x1 +tmp*x2)/(1.d0 +tmp)

    !     note: do-loop code
        do i = 1, maxloop
            f = matching_d1(energy, l)
            sign = fx1*f
            if(sign > 0.d0) then
                x1 = energy
                fx1 = f
                dx = abs(x1*tiny3)
                energy = x1 +dx
                if(.not. energy < x2) then
                    tmp = abs(fx1)/abs(fx2)
                    energy = (x1 +tmp*x2)/(1.d0 +tmp)
                    return
                end if
                f = matching_d1(energy, l)
                df = (f -fx1)/dx
                sign = fx1*f
                if(sign < 0.d0) then
                    x2 = energy
                    fx2 = f
                    tmp = abs(fx1)/abs(fx2)
                    energy = (x1 +tmp*x2)/(1.d0 +tmp)
                    return
                else if(sign == 0.d0) then
                    return
                end if

            else if(sign < 0.d0) then
                x2 = energy
                fx2 = f
                dx = abs(x1*tiny3)
                energy = x2 -dx
                if(.not. energy > x1) then
                    tmp = abs(fx1)/abs(fx2)
                    energy = (x1 +tmp*x2)/(1.d0 +tmp)
                    return
                end if
                f = matching_d1(energy, l)
                df = (fx2 -f)/dx
                sign = fx1*f
                if(sign < 0.d0) then
                    x1 = energy
                    fx1 = f
                    tmp = abs(fx1)/abs(fx2)
                    energy = (x1 +tmp*x2)/(1.d0 +tmp)
                    return
                else if(sign == 0.d0) then
                    return
                end if
            else if(sign == 0.d0) then
                return
            end if

            dx = -f/df
            tmp = energy +dx
            if(energy == tmp) then
                return
            else if(1.d0 +abs(dx/energy) == 1.d0) then
                energy = tmp
                return
            else if(tmp < x2 .and. tmp > x1) then
                energy = tmp
            else
                tmp = abs(fx1)/abs(fx2)
                energy = (x1 +tmp*x2)/(1.d0 +tmp)
            end if
            if(i == maxloop) print *, "Warning boundary #727: Do-loop arrive at 'maxloop'. "
        end do
    end if
end subroutine SUB_initial_bound
! find other bound states --------------------------
subroutine SUB_floquet_bound(m, energy)
    use math_const, only: i => math_i
    character(60), parameter :: form_out = "(1A15, 3ES15.3)"
    complex(dp) :: z1, z2, fz1, df, ddf, zm, zp, dr, fzm, fzp, dz
    complex(dp), intent(inout) :: energy
    integer(i4) :: j
    integer(i4), intent(in) :: m
    real(dp) :: tmp

    z1 = energy
    dz = 1.d0*i
    do j = 1, maxloop
        fz1 = matching_dF(m, z1)
        write(file_log, form_out) "Assumed E: ", z1, abs(fz1)

        dr = abs(z1)*dz/abs(dz)*tiny3
        zm = z1 -dr
        zp = z1 +dr
        fzm = matching_dF(m, zm)
        fzp = matching_dF(m, zp)
        df = (-1.d0/2.d0*fzm +1.d0/2.d0*fzp)/dr
        ddf = (fzm -2.d0*fz1 +fzp)/dr**2.d0
        dz = -fz1/(df -ddf/2.d0/df*fz1)
        tmp = abs(fz1*ddf/df**2.d0)

        z2 = z1 +dz
        if(real(1.d0 +abs(fz1), kind = dp) == real(1.d0, kind = dp)) then
            if(cmplx(z2, kind = sp) == cmplx(z1, kind = sp)) then
                energy = z2
                exit
            end if
        end if

        z1 = z2
        if(j == maxloop) then
            energy = z2
            print *, "Warning boundary #769: Do-loop arrive at 'maxloop'. "
        end if
    end do
end subroutine SUB_floquet_bound
! find bound vector --------------------------------
subroutine SUB_floquet_vector(m, energy, y)
    complex(dp), intent(in) :: energy
    complex(dp), intent(out) :: y(1:)
    integer(i4), intent(in) :: m

    if(allocated(k)) deallocate(k)
    if(allocated(R)) deallocate(R)
    if(allocated(B)) deallocate(B)
    allocate(k(-F:F))
    allocate(R(-F:F, -F:F))
    allocate(B(-F:F, -F:F))
    call ssub_channel_k(m, energy)
    call ssub_lngth_R(energy)
    call ssub_matching_B
    call ssub_eigen_vector(B, y)
    if(allocated(k)) deallocate(k)
    if(allocated(R)) deallocate(R)
    if(allocated(B)) deallocate(B)
end subroutine SUB_floquet_vector
! end sub-calculation ------------------------------










! ==================================================
! TEST
! ==================================================
! end test -----------------------------------------


! ==================================================
! PROCESS
! ==================================================
! matching boundary condition ----------------------
subroutine PROC_matching(mi, li)
    use math_const, only: i => math_i, pi => math_pi
    use hamiltonian, only: Amp_grid
    character(30), parameter :: &
        form_gen    = "(1I15, 9999ES25.10)", &
        form_label  = "('# ', 13X, 1A25, 9999I25)", &
        form_num    = "(1I4.4)", &
        form_out    = "(1A15, 1ES15.3, 1ES15.3)"
    complex(dp) :: x
    complex(dp), allocatable :: y(:)
    integer(i4) :: j
    integer(i4), intent(in) :: mi, li
    integer(i4), save :: f0
    real(dp) :: tmp
    real(qp) :: norm

    if(mi == 0) then
        open(file_cmplx, file = "output/bound.data")
        if(allocated(Scatt)) deallocate(Scatt)
        allocate(Scatt(0:M, 0:L))
        call SUB_initial_bound(li, tmp)
!         todo: find multi bound state
        Scatt(mi, li) = tmp
        write(file_log, form_out) "Bound State: ", tmp
        write(file_cmplx, form_gen) li, Amp_grid(mi), tmp, 0.d0
        write(file_cmplx, form_gen)
        write(file_cmplx, form_gen)
        do j = 1, maxloop
            tmp = real(Scatt(0, li)) +dble(j)*Freq
            if(tmp > 0.d0) then
                f0 = j
                exit
            end if
        end do

    else if(mi > 0) then
        x = Scatt(mi -1, li)
        if(mi == 1) then
!             x = x -i*9.d-11
            open(file_branch1, file = "output/branch1.data")
            open(file_branch2, file = "output/branch2.data")
            open(file_ion, file = "output/ion.data")
            write(file_branch1, form_label) "line number", (j +2, j = 1, 2*F +1)
            write(file_branch1, form_label) "floquet number", (j, j = -F, F)
            write(file_branch2, form_label) "line number", (j +2, j = 1, F -f0 +1)
            write(file_branch2, form_label) "floquet number", (j, j = f0, F)
        end if

!         note: do-loop code
        do j = 1, maxloop
            call SUB_floquet_bound(mi, x)
            if(aimag(x) > 0.d0) then
                x = conjg(x)
            else
                exit
            end if
            if(j == maxloop) print *, "Warning boundary #115: Do-loop arrive at 'maxloop'. "
        end do
        if(.not. aimag(x) < 0.d0) stop "Error boundary #673: Ionization rate is minus. There is somthing wrong. "
        Scatt(mi, li) = x
        write(file_log, form_out) "Bound State: ", x
        write(file_cmplx, form_gen) li, Amp_grid(mi), real(x), aimag(x)
        write(file_ion, form_gen) li, Amp_grid(mi), real(x) -real(Scatt(0, li)), -2.d0*aimag(x)

        if(allocated(y)) deallocate(y)
        allocate(y(-F:F))
        call SUB_floquet_vector(mi, x, y)
        write(file_branch1, form_gen) li, Amp_grid(mi), (abs(y(j))**2.d0, j = -F, F)
        norm = 0.d0
        do j = f0, F
            norm = norm +abs(y(j))**2.d0
        end do
        write(file_branch2, form_gen) li, Amp_grid(mi), (abs(y(j))**2.d0/norm, j = f0, F)
        if(allocated(y)) deallocate(y)
    end if
end subroutine PROC_matching
! out ----------------------------------------------
subroutine PROC_boundary_out
    if(allocated(Scatt)) deallocate(Scatt)
    close(file_cmplx)
    close(file_branch1)
    close(file_branch2)
    close(file_ion)
end subroutine PROC_boundary_out
! end out ------------------------------------------
end module boundary
