module boundary
    use kind_const
    use global
    implicit none
    complex(dp), save, allocatable, protected :: Scatt(:, :)
    complex(dp), save, allocatable, protected, private :: &
        R(:, :), CH(:, :), U1(:, :), U2(:, :), W0(:, :), A(:, :)
    integer(i4), parameter, private :: file_bound = 201, addblock = 10, addcount = 3
    real(dp), save, allocatable, protected, private :: d(:)
    real(dp), save, protected, private :: A0, alpha0, Ep
contains


! ==================================================
! FUNCTIONS (SINGLE)
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
! end functions (single) ---------------------------








! ==================================================
! FUNCTIONS (FLOQUET)
! ==================================================
! floquet R matrix ---------------------------------
function matching_RF(f1, f2, energy)
    use basis, only: HF, EF
    complex(dp) :: matching_RF
    complex(dp), intent(in) :: energy
    complex(qp) :: sum
    integer(i4) :: i
    integer(i4), intent(in) :: f1, f2

    sum = 0.d0
    do i = 1, (2*F +1)*N
        sum = sum +HF(N, f1, i)*HF(N, f2, i)/(EF(i) -energy)
    end do
    matching_RF = sum/(2.d0*Ba)
end function matching_RF
! transfer C ---------------------------------------
function transfer_f(l, x)
    integer(i4) :: m, m1, m2, dm, j, count
    integer(i4), intent(in) :: l
    real(dp), intent(in) :: x
    real(dp) :: transfer_f, inp1, inp2, tmp1, tmp2
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
    transfer_f = sum
end function transfer_f
! determinant matrix WD1 ---------------------------
function det_WD1(x)
    use mylinear, only: det_cmplx
    complex(dp) :: det_WD1, tmp1, tmp2
    complex(dp), allocatable :: work(:, :)
    complex(dp), intent(in)  :: x
    integer(i4) :: i
    if(allocated(work)) deallocate(work)
    allocate(work(-F:F, -F:F))
    work(:, :) = -W0(:, :)
    tmp1 = x**2.d0
    do i = -F, F
        tmp2 = 2.d0*d(i)*x
        work(i, i) = work(i, i) +tmp1 +tmp2
    end do
    det_WD1 = det_cmplx(work)
    if(allocated(work)) deallocate(work)
end function det_WD1
! solve a eigen-vector of zero eigen-value ---------
subroutine eigen_vector(A, v, dA)
    use mylinear, only: solve_cmplx_vec
    complex(dp), allocatable :: work(:, :)
    complex(dp), intent(in) :: A(1:, 1:)
    complex(dp), intent(out) :: v(1:), dA
    complex(qp) :: sum
    integer(i8) :: i1, j1, i2, j2, n, n0, n1, k, info
    real(qp) :: norm

    n = size(A(1:, 1))
    n0 = int(n/2) +1
    n1 = n0
    if(allocated(work)) deallocate(work)
    allocate(work(1:(n -1), 1:(n -1)))
    if(.not. (size(A(1, 1:)) == n .and. size(v(1:)) == n)) then
        stop "Error boundary #208: A and v have different dimention. "
    end if

    do k = 1, 2
        v(:) = 0.d0
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
                    v(i2) = -A(i1, j1)
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
        v(n) = 1.d0

        info = 0
        call solve_cmplx_vec(work, v(1:(n -1)), info)

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
        v(i1) = v(i1 -1)
    end do
    v(n1) = 1.d0
    if(n1 /= n0) then
        dA = 0.d0
    else
        sum = 0.d0
        do i1 = 1, n
            sum = sum +A(n0, i1)*v(i1)
        end do
        if(1.d0 +abs(sum) == 1.d0) then
            dA = 0.d0
        else
            dA = sum
        end if
    end if

    norm = 0.d0
    do i1 = 1, n
        norm = norm +abs(v(i1))**2.d0
    end do
    v(:) = v(:)/norm**0.5d0

    if(n1 /= n0) then
        v(:) = v(:)*abs(v(n0))/v(n0)
    end if
end subroutine eigen_vector
! outer region wave function -----------------------
function wave_G(r, j1, j2)
    use math_const, only: i => math_i
    complex(dp) :: wave_G, p
    integer(i4), intent(in) :: j1, j2
    real(dp), intent(in) :: r
    p = CH(j2, 2)
    wave_G = exp(i*p*r)*A(j1, j2)
end function wave_G
! outer region differcial wave function ------------
function wave_dG(r, j1, j2)
    use math_const, only: i => math_i
    complex(dp) :: wave_dG, p
    integer(i4), intent(in) :: j1, j2
    real(dp), intent(in) :: r
    p = CH(j2, 2)
    wave_dG = i*p*exp(i*p*r)*A(j1, j2)
end function wave_dG
! matching parameter dF ----------------------------
function matching_dF(m, l, energy)
    use fgsl, only: fgsl_sf_bessel_ksl_scaled
    character(60), parameter :: form_out = '(1A15, 3ES15.3)'
    complex(dp) :: matching_dF
    complex(dp), intent(in) :: energy
    integer(i4), intent(in) :: m, l

    if(allocated(U1)) deallocate(U1)
    if(allocated(U2)) deallocate(U2)
    if(allocated(W0)) deallocate(W0)
    if(allocated(d))  deallocate(d)
!     allocate(U1(-F:F, -F:F))
    allocate(U2(-F:F, -F:F))
    allocate(W0(-F:F, -F:F))
    allocate( d(-F:F))
    call ssub_channel(m, energy)

    if(allocated(R))  deallocate(R)
    allocate(R(-F:F, -F:F))
    ! todo: matrix size option
    ! todo: angular momantum elements
    call ssub_lngth_R(energy)
    call ssub_vlcty_R(l)
    call ssub_propa_R
    call ssub_modif_R

    if(allocated(A))  deallocate(A)
    allocate(A(-F:F, -F:F))
    call ssub_A
    if(allocated(U1)) deallocate(U1)
    if(allocated(U2)) deallocate(U2)
    if(allocated(W0)) deallocate(W0)

    call ssub_matching(matching_dF)
    if(allocated(R))  deallocate(R)
    if(allocated(d))  deallocate(d)
    if(allocated(A))  deallocate(A)
!     write(file_log, form_out) "Assumed E: ", energy, abs(matching_dF)
end function matching_dF
! end functions (floquet) --------------------------










! ==================================================
! SUB-SUB-CALCULATION
! ==================================================
! channel p ----------------------------------------
subroutine ssub_channel(m, energy)
    use math_const, only: i => math_i, pi => math_pi
    use mylinear, only: diag_her
    use hamiltonian, only: Amp_grid
    character(30), parameter :: &
        form_ch  = "(1I15, 2ES25.10)", &
        form_num = "(1I4.4)", &
        form_out = "(1A15, 1ES15.3, 1ES15.3, 'i')"
    complex(dp) :: tmp1, tmp2, dx, fx
    complex(dp), allocatable :: ksq(:), work(:, :), fwork(:, :)
    complex(dp), intent(in) :: energy
    complex(qp) :: sum
    integer(i4) :: j1, j2, k1, k2, n, check
    integer(i4), intent(in) :: m
    real(dp) :: tmp3

    A0 = -Amp_grid(m)/Freq
    alpha0 = Amp_grid(m)/Freq**2.d0
    Ep = A0**2.d0/4.d0

    if(allocated(ksq)) deallocate(ksq)
    allocate(ksq(-F:F))
    tmp1 = energy -A0**2.d0/4.d0
    do j2 = -F, F
        tmp2 = 2.d0*(tmp1 +dble(j2)*Freq)
        ksq(j2) = tmp2
    end do
    U2(:, :) = 0.d0
    do j2 = -F +1, F
        U2(j2 -1, j2) = i/2.d0*A0
    end do
    call diag_her(U2, d)
    do j2 = -F, F
        do j1 = -F, F
            sum = 0.d0
            do k1 = -F, F
                sum = sum +conjg(U2(k1, j1))*ksq(k1)*U2(k1, j2)
            end do
            W0(j1, j2) = sum
        end do
    end do

    if(m == 1) then
        CH(:, :) = 0.d0
        do j1 = -F, F
            tmp1 = ksq(j1)**0.5d0
            tmp3 = aimag(tmp1)
            if(tmp3 > 0.d0) then
                CH(j1, 1) = tmp1
                CH(j1, 2) = -tmp1
            else if(tmp3 < 0.d0) then
                CH(j1, 1) = -tmp1
                CH(j1, 2) = tmp1
            end if
        end do
    end if

    if(allocated(ksq))   deallocate(ksq)
    if(allocated(work))  deallocate(work)
    if(allocated(fwork)) deallocate(fwork)
    allocate(work(-F:F, 1:2))
    allocate(fwork(-F:F, 1:2))
    work(:, :) = CH(:, :)
    do j2 = 1, 2
        do j1 = -F, F
            fwork(j1, j2) = det_WD1(CH(j1, j2))
        end do
    end do
!     note: do-loop code
    do n = 1, maxloop
        check = 1
        do j2 = 1, 2
            do j1 = -F, F
                fx = fwork(j1, j2)
                if(1.d0 +abs(fx) == 1.d0) cycle
                sum = 1.d0
                do k2 = 1, 2
                    do k1 = -F, F
                        if(k2 == j2 .and. k1 == j1) cycle
                        sum = sum*(work(j1, j2) -work(k1, k2))
                    end do
                end do
                dx = -fx/sum
                CH(j1, j2) = work(j1, j2) +dx
                if(.not. CH(j1, j2) == work(j1, j2)) check = check*0
                fwork(j1, j2) = det_WD1(CH(j1, j2))
            end do
        end do
        if(check == 1) exit
        if(n == maxloop) then
            tmp1 = maxval(abs(fwork(:, :)))
!             print *, "Warning boundary #412: Do-loop arrive at 'maxloop'. "
!             print *, "  > Following value should be zero. ", abs(tmp1)
!             note: for convinience.
            tmp1 = tmp1**2.d0
            if(1.d0 +tmp1 /= 1.d0) then
                print *, "Warning boundary #412: Do-loop arrive at 'maxloop'. "
                print *, "  > Following value should be zero. ", abs(tmp1)
            end if
        end if
        work(:, :) = CH(:, :)
    end do
    if(allocated(work))  deallocate(work)
    if(allocated(fwork)) deallocate(fwork)
!     write(file_log, form_out)   "channel p: ", CH(0, 1)
!     write(file_log, form_out)   "",            CH(0, 2)
end subroutine ssub_channel
! length gauge R -----------------------------------
subroutine ssub_lngth_R(energy)
    complex(dp), intent(in) :: energy
    integer(i4) :: i, j

    R(:, :) = 0.d0
    do i = -F, F
        do j = -F, F
            R(i, j) = matching_RF(i, j, energy)
        end do
!         R(j, j) = matching_RF(j, j, energy)
    end do
end subroutine ssub_lngth_R
! velocity gauge R ---------------------------------
subroutine ssub_vlcty_R(l)
    use mylinear, only: inverse_cmplx
    complex(dp), allocatable :: rR(:, :), Y(:, :)
    complex(qp) :: sum
    integer(i4) :: i, j, k
    integer(i4), intent(in) :: l
    real(dp), allocatable :: C(:, :), dC(:, :)

    if(allocated(rR)) deallocate(rR)
    if(allocated(Y))  deallocate(Y)
    if(allocated(C))  deallocate(C)
    if(allocated(dC)) deallocate(dC)
    allocate(rR(-F:F, -F:F))
    allocate( Y(-F:F, -F:F))
    allocate( C(-F:F, -F:F))
    allocate(dC(-F:F, -F:F))

    rR(:, :) = R(:, :)
    R(:, :)  = 0.d0
    call inverse_cmplx(rR)
    do j = -F, F
        do i = -F, F
            C (i, j) = transfer_f(i -j, Ba)
            dC(i, j) = A0*0.5d0*(transfer_f(i -j -1, Ba) -transfer_f(i -j +1, Ba))
        end do
    end do
    do j = -F, F
        do i = -F, F
            sum = 0.d0
            do k = -F, F
                sum = sum +C(i, k)*rR(k, j)/Ba
            end do
            Y(i, j) = sum
        end do
    end do
    Y(:, :) = Y(:, :) +dC(:, :)
    call inverse_cmplx(Y)
    do i = -F, F
        do j = -F, F
            sum = 0.d0
            do k = -F, F
                sum = sum +C(i, k)*Y(k, j)/Ba
            end do
            R(i, j) = sum
        end do
    end do

    if(allocated(rR)) deallocate(rR)
    if(allocated(Y))  deallocate(Y)
    if(allocated(C))  deallocate(C)
    if(allocated(dC)) deallocate(dC)
end subroutine ssub_vlcty_R
! propagate R-matrix -------------------------------
subroutine ssub_propa_R
    use math_const, only: i => math_i
    use mylinear, only: diag_her, inverse_cmplx
    complex(dp) :: tmp1, tmp2, tmp3, tmp4
    complex(dp), allocatable :: &
        Y1(:, :), Y2(:, :), V(:, :), tildeG(:, :), G(:, :, :)
    complex(qp) :: sum
    integer(i4) :: j1, j2, k1, k2, l
    real(dp) :: beta, Bh, Bl, Br, lambda, tmp0
    real(dp), allocatable :: lsq(:), lsqold(:)

    ! modify R-matrix to bar R-matrix
    if(allocated(Y1)) deallocate(Y1)
    if(allocated(Y2)) deallocate(Y2)
    allocate(Y1(-F:F, -F:F))
    allocate(Y2(-F:F, -F:F))
    do j2 = -F, F
    tmp4 = exp(-i*d(j2)*Ba)
        do j1 = -F, F
            tmp1 = exp(i*d(j1)*Ba)
            sum  = 0.d0
                do k2 = -F, F
                tmp3 = U2(k2, j2)
                    do k1 = -F, F
                        tmp2 = conjg(U2(k1, j1))*R(k1, k2)
                        sum  = sum +tmp1*tmp2*tmp3*tmp4
                    end do
                end do
            Y2(j1, j2) = sum
        end do
    end do
    do j2 = -F, F
        Y1(:, j2) = i*Ba*Y2(:, j2)*d(j2)
    end do
    do j2 = -F, F
        Y1(j2, j2) = Y1(j2, j2) +1.d0
    end do
    call inverse_cmplx(Y1)
    do j2 = -F, F
        do j1 = -F, F
            sum = 0.d0
            do k1 = -F, F
                sum  = sum +Y1(j1, k1)*Y2(k1, j2)
            end do
            R(j1, j2) = sum
        end do
    end do
    if(allocated(Y1)) deallocate(Y1)
    if(allocated(Y2)) deallocate(Y2)

    if(Ba == Bap) return
    if(allocated(lsq))    deallocate(lsq)
    if(allocated(lsqold)) deallocate(lsqold)
    allocate(lsq(-F:F))
    allocate(lsqold(-F:F))
    beta = Ba*Bbeta
    Bh = Ba
    Bl = 0.d0
    Br = Bl +Bh
    lsqold(:) = 0.d0
    lsq(:) = Bh

!     note: do-loop code
    do l = 1, maxloop
        ! assumed sector size
        sum = 0.d0
        do j1 = -F, F
            sum = sum +abs((lsq(j1) -lsqold(j1))/Bh)
        end do
        sum = sum/dble(2*F +1)
        Bh = beta*sum**(-1.d0/3.d0)
        Bl = Br
        Br = Bl +Bh
        if(Br > Bap) Br = Bap
        lsqold(:) = lsq(:)
        lsq(:) = 0.d0

        ! green-function G-matrix
        if(allocated(V))   deallocate(V)
        allocate(V(-F:F, -F:F))
        do j2 = -F, F
            tmp1 = exp(-i*d(j2)*Bl)
            tmp2 = exp(-i*d(j2)*Br)
            do j1 = -F, F
                V(j1, j2) = 0.5d0*(exp(d(j1)*Bl)*W0(j1, j2)*tmp1 +exp(d(j1)*Br)*W0(j1, j2)*tmp2)
            end do
        end do
        do j1 = -F, F
            V(j1, j1) = V(j1, j1) +d(j1)**2.d0
        end do
        call diag_her(V, lsq)

    !     note: 1[ll], 2[lr], 3[rl], 4[rr]
        if(allocated(tildeG)) deallocate(tildeG)
        allocate(tildeG(1:4, -F:F))
        tildeG(:, :) = 0.d0
        do j1 = -F, F
            if(lsq(j1) > 0.d0) then
                lambda = lsq(j1)**0.5d0
                tmp1 = -lambda*sinh(lambda*(Br -Bl))
                tildeG(1, j1) = 1.d0/tmp1
                tildeG(2, j1) = cosh(lambda*(Bl -Br))/tmp1
                tildeG(3, j1) = cosh(lambda*(Br -Bl))/tmp1
                tildeG(4, j1) = 1.d0/tmp1
            else if(lsq(j1) < 0.d0) then
                lambda = (-lsq(j1))**0.5d0
                tmp0 = lambda*sin(lambda*(Br -Bl))
                tildeG(1, j1) = 1.d0/tmp0
                tildeG(2, j1) = cos(lambda*(Bl -Br))/tmp0
                tildeG(3, j1) = cos(lambda*(Br -Bl))/tmp0
                tildeG(4, j1) = 1.d0/tmp0
            end if
        end do
        if(allocated(G)) deallocate(G)
        allocate(G(1:4, -F:F, -F:F))
        G(:, :, :) = 0.d0
        do j2 = -F, F
            do j1 = -F, F
                do k2 = 1, 4
                    sum = 0.d0
                    do k1 = -F, F
                        sum = sum +V(j1, k1)*tildeG(k2, k1)*conjg(V(j2, k1))
                    end do
                    G(k2, j1, j2) = sum
                end do
            end do
        end do
        if(allocated(V))      deallocate(V)
        if(allocated(tildeG)) deallocate(tildeG)

        ! propagate R-matrix
        if(allocated(Y1)) deallocate(Y1)
        if(allocated(Y2)) deallocate(Y2)
        allocate(Y1(-F:F, -F:F))
        allocate(Y2(-F:F, -F:F))
        Y1(:, :) = R(:, :) -G(1, :, :)
        call inverse_cmplx(Y1)
        Y2(:, :) = 0.d0
        do j2 = -F, F
            do j1 = -F, F
                sum = 0.d0
                do k2 = -F, F
                    tmp1 = G(2, k2, j2)
                    do k1 = -F, F
                        sum = sum +G(3, j1, k1)*Y1(k1, k2)*tmp1
                    end do
                end do
                Y2(j1, j2) = sum
            end do
        end do
        R(:, :) = -(G(4, :, :) +Y2(:, :))
        if(allocated(Y1)) deallocate(Y1)
        if(allocated(Y2)) deallocate(Y2)
        if(Br == Bap) exit
        if(l == maxloop) print *, "Warning boundary #658: Do-loop arrive at 'maxloop'. "
    end do
    if(allocated(lsqold)) deallocate(lsqold)
    if(allocated(lsq))    deallocate(lsq)
end subroutine ssub_propa_R
! modif R ------------------------------------------
subroutine ssub_modif_R
    use math_const, only: i => math_i
    character(60), parameter :: form_out = "(1A15, 1ES15.3, 1ES15.3, 'i')"
    complex(dp) :: tmp
    complex(dp), allocatable :: barR(:, :)
    integer(i4) :: j1, j2

    ! modify bar R-matrix to tilde R-matrix
    if(allocated(barR)) deallocate(barR)
    allocate(barR(-F:F, -F:F))
    barR(:, :) = R(:, :)
    R(:, :) = 0.d0
    do j2 = -F, F
        tmp = exp(i*d(j2)*Bap)
        do j1 = -F, F
            R(j1, j2) = exp(-i*d(j1)*Bap)*barR(j1, j2)*tmp
        end do
    end do
    if(allocated(barR)) deallocate(barR)
!     write(file_log, form_out) "(modified) R: ", R(0, 0)
end subroutine ssub_modif_R
! matrix A -----------------------------------------
subroutine ssub_A
    complex(dp) :: x, tmp1, tmp2
    complex(dp), allocatable :: work(:, :), y(:)
    integer(i4) :: i, j
    real(dp) :: tmp3

    if(allocated(work)) deallocate(work)
    if(allocated(y))    deallocate(y)
    allocate(work(-F:F, -F:F))
    allocate(y(-F:F))
    do j = -F, F
        work(:, :) = -W0(:, :)
        x = CH(j, 2)
        tmp1 = x**2.d0
        do i = -F, F
            tmp2 = 2.d0*d(i)*x
            work(i, i) = work(i, i) +tmp1 +tmp2
        end do
        call eigen_vector(work, y, tmp2)
        tmp3 = abs(tmp2)
        if(.not. real(1.d0 +tmp3, kind = sp) == real(1.d0, kind = sp)) then
!             print *, "Warning boundary #760: dA of ssub_A is not zero. "
!             print *, "  > Following value should be zero. ", tmp3
!             note: for convinience
            tmp3 = tmp3**2.d0
            if(1.d0 +tmp3 == 1.d0) then
                print *, "Warning boundary #760: dA of ssub_A is not zero. "
                print *, "  > Following value should be zero. ", tmp3
            end if
        end if
        A(:, j) = y(:)
    end do
    if(allocated(work)) deallocate(work)
    if(allocated(y))    deallocate(y)
end subroutine ssub_A
! matching boundary condition -----------------------
subroutine ssub_matching(dF)
    use math_const, only: i => math_i
    use mylinear, only: det_cmplx
    character(60), parameter :: form_out = "(1A15, 1ES15.3, 1ES15.3, 'i')"
    complex(dp), allocatable :: Gap1(:, :), Gap2(:, :), work(:, :), y(:)
    complex(dp), intent(out) :: dF
    complex(qp) :: sum
    integer(i4) :: j1, j2, k

    if(allocated(Gap1)) deallocate(Gap1)
    if(allocated(Gap2)) deallocate(Gap2)
    if(allocated(work)) deallocate(work)
    if(allocated(y))    deallocate(y)
    allocate(Gap1(-F:F, -F:F))
    allocate(Gap2(-F:F, -F:F))
    allocate(work(-F:F, -F:F))
    allocate(y(-F:F))

    do j2 = -F, F
        do j1 = -F, F
            Gap1(j1, j2) = wave_G(Bap, j1, j2)
            Gap2(j1, j2) = wave_dG(Bap, j1, j2)
        end do
    end do
    do j1 = -F, F
        Gap2(j1, :) = Gap2(j1, :) +i*d(j1)*Gap1(j1, :)
    end do
    do j2 = -F, F
        do j1 = -F, F
            sum = 0.d0
            do k = -F, F
                sum = sum +R(j1, k)*Bap*Gap2(k, j2)
            end do
            work(j1, j2) = sum
        end do
    end do
    work(:, :) = Gap1(:, :) -work(:, :)

    dF = det_cmplx(work)
!     call eigen_vector(work, y, dF)
    if(allocated(Gap1)) deallocate(Gap1)
    if(allocated(Gap2)) deallocate(Gap2)
    if(allocated(work)) deallocate(work)
    if(allocated(y))    deallocate(y)
!     write(file_log, form_out) "matching: ", dF
end subroutine ssub_matching
! end sub-sub-calculation --------------------------









! ==================================================
! SUB-CALCULATION
! ==================================================
! find initial bound state -------------------------
subroutine SUB_initial_bound(l, energy)
    use basis, only: E1
    integer(i4), intent(in) :: l
    real(dp), intent(out) :: energy
    real(sp) :: sign
    real(dp) :: x1, x2, dx, dxold, f, fx1, fx2, df, tmp
    integer(i4) :: i

    x1 = E1(1)
    x2 = E1(N)
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
end subroutine SUB_initial_bound
! find other bound states --------------------------
subroutine SUB_floquet_bound(m, l, energy)
    use math_const, only: i => math_i
    character(60), parameter :: form_out = '(1A15, 3ES15.3)'
    complex(dp) :: z1, z2, fz, df, ddf, zm2, zm1, zp1, zp2, dr, fzm2, fzm1, fzp1, fzp2, dz
    complex(dp), intent(inout) :: energy
    integer(i4) :: j
    integer(i4), intent(in) :: m, l

    z1 = energy
    dz = 1.d0
    do j = 1, maxloop
        fz = matching_dF(m, l, z1)
        write(file_log, form_out) "Assumed E: ", z1, abs(fz)

        dr = abs(z1)*dz/abs(dz)*tiny3
        zm2 = z1 -2.d0*dr
        zm1 = z1 -dr
        zp1 = z1 +dr
        zp2 = z1 +2.d0*dr
        fzm2 = matching_dF(m, l, zm2)
        fzm1 = matching_dF(m, l, zm1)
        fzp1 = matching_dF(m, l, zp1)
        fzp2 = matching_dF(m, l, zp2)
        df = (1.d0/12.d0*fzm2 -8.d0/12.d0*fzm1 +8.d0/12.d0*fzp1 -1.d0/12.d0*fzp2)/dr
        ddf = (-1.d0/12.d0*fzm2 +16.d0/12.d0*fzm1 -30.d0/12.d0*fz +16.d0/12.d0*fzp1 -1.d0/12.d0*fzp2)/dr**2.d0
        dz = -fz/(df -ddf/2.d0/df*fz)
        print *, abs(fz), abs(df), abs(ddf)

        z2 = z1 +dz
        if(real(1.d0 +abs(fz), kind = sp) == real(1.d0, kind = sp) .and. z2 == z1) then
            energy = z2
            return
        end if
        z1 = z2
        if(j == maxloop) then
            energy = z2
            print *, "Warning boundary #769: Do-loop arrive at 'maxloop'. "
        end if
    end do
end subroutine SUB_floquet_bound
! end sub-calculation ------------------------------










! ==================================================
! PROCESS
! ==================================================
! boundary matrix ----------------------------------
subroutine PROC_matching(mi, li)
    use math_const, only: i => math_i, pi => math_pi
    use hamiltonian, only: Amp_grid
    character(30), parameter :: &
        form_bound = "(1I15, 3ES25.10)", &
        form_ch    = "(1I15, 2ES25.10)", &
        form_num   = "(1I4.4)", &
        form_out   = "(1A15, 1ES15.3, 1ES15.3)"
    character(4) :: num_mi, num_li
    complex(dp) :: x
    integer(i1), parameter :: file_ch = 101
    integer(i4) :: j
    integer(i4), intent(in) :: mi, li
    real(dp) :: tmp

    if(mi == 0) then
        open(file_bound, file = "output/bound.data")
        if(allocated(Scatt)) deallocate(Scatt)
        if(allocated(CH))    deallocate(CH)
        allocate(Scatt(0:M, 0:L))
        allocate(CH(-F:F, 1:2))
        call SUB_initial_bound(li, tmp)
!         todo: find multi bound state
        Scatt(mi, li) = tmp
        write(file_log, form_out) "Bound State: ", tmp
        write(file_bound, form_bound) li, Amp_grid(mi), tmp, 0.d0
        write(file_bound, form_bound)
        write(file_bound, form_bound)

    else if(mi > 0) then
        x = Scatt(mi -1, li)
        if(aimag(x) == 0.d0) x = x -i*abs(x)*tiny3

!         note: do-loop code
        do j = 1, maxloop
            call SUB_floquet_bound(mi, li, x)
            if(aimag(x) > 0.d0) then
                x = conjg(x)
            else if(aimag(x) == 0.d0) then
                x = x -i*abs(x)*tiny3
            else
                exit
            end if
            if(j == maxloop) print *, "Warning boundary #115: Do-loop arrive at 'maxloop'. "
        end do

        if(.not. aimag(x) < 0.d0) stop "Error boundary #046: Ionization rate is minus. There is somthing wrong. "
        Scatt(mi, li) = x
        write(file_log, form_out) "Bound State: ", x
        write(file_bound, form_bound) li, Amp_grid(mi), real(x), aimag(x)
        write(num_mi, form_num) mi
        write(num_li, form_num) li
        open(file_ch, file = "output/ch_"//num_li//num_mi//".data")
        do j = -F, F
            write(file_ch, form_ch) j, real(CH(j, 2)), aimag((CH(j, 2)))
        end do
        close(file_ch)
    end if
end subroutine PROC_matching
! out ----------------------------------------------
subroutine PROC_boundary_out
    if(allocated(Scatt)) deallocate(Scatt)
    if(allocated(CH))    deallocate(CH)
    close(file_bound)
end subroutine PROC_boundary_out
! end out ------------------------------------------
end module boundary
