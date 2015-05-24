module boundary
    use kind_type
    use global 
    implicit none
    real   (dp), protected, private :: P 
    real   (dp), allocatable, protected, private :: d(:)
    complex(dp), allocatable, protected, private :: CH(:, :), U1(:, :), U2(:, :), W(:, :)
    complex(dp), save, allocatable, protected :: &
        lngth_R(:, :, :), vlcty_R(:, :, :), propa_R(:, :, :), modif_R(:, :, :), &
        A(:, :), S(:, :, :)
    ! temperater 
    complex(dp), parameter, private :: gamma = -0.5d0*(0.d0, 1.d0)*1.d3
    real   (sp), parameter, private :: tiny1  = 1.e-16, tiny2 = 1.e-8, tiny3 = 1.e-6 
    integer(i4), parameter, private :: cal = 500 

contains 


! ==================================================
! FUNCTION 
! ==================================================
! function R (real energy) -------------------------
function real_R(f1, f2, energy)
    use basis, only: H, E 
    integer(i4), intent(in) :: f1, f2 
    real   (dp), intent(in) :: energy
    real   (dp) :: real_R 
    real   (qp) :: sum 
    integer(i4) :: i
    sum = 0.d0 
    do i = 1, (2*F +1)*N 
        sum = sum +H(N, f1, i)*H(N, f2, i)/(E(i) -energy)
    end do 
    real_R = sum/(2.d0*Mass*Ra)
end function real_R
! function R (complex energy) ----------------------
function cmplx_R(f1, f2, energy)
    use basis, only: H, E 
    integer(i4), intent(in) :: f1, f2 
    complex(dp), intent(in) :: energy
    complex(dp) :: cmplx_R 
    complex(qp) :: sum 
    integer(i4) :: i
    sum = 0.d0 
    do i = 1, (2*F +1)*N 
        sum = sum +H(N, f1, i)*H(N, f2, i)/(E(i) -energy)
    end do 
    cmplx_R = sum/(2.d0*Mass*Ra)
end function cmplx_R
! function C ---------------------------------------
function transfer_C(n, er)
    use fgsl, only: fgsl_sf_bessel_jcn
    integer(i4), intent(in) :: n 
    real   (dp), intent(in) :: er 
    real   (dp) :: in1, in2, tmp0, tmp1, tmp2, tmpm, tmpr, dsum, transfer_C
    real   (qp) :: sum 
    integer(i4) :: m, dm, m0, m1, count
    in1 = Amp**2.d0/Freq**3.d0/8.d0 
    in2 = Amp/Freq*er 
    if(mod(abs(n), 2_i4) == 0) then 
        m0 = min(0, n/2)
        m1 = max(0, n/2) 
    else if(mod(abs(n), 2_i4) == 1) then 
        if(n > 0) then 
            m0 = 0 
            m1 = (n +1)/2
        else if(n < 0) then 
            m0 = (n -1)/2
            m1 = 0 
        end if 
    end if 
    sum  = 0.d0 
    do m = m0, m1 
        if(m == m0) tmp1 = fgsl_sf_bessel_jcn(m, in1)*fgsl_sf_bessel_jcn(n -2_i4*m, in2)
        tmp2 = fgsl_sf_bessel_jcn(m, in1)*fgsl_sf_bessel_jcn(n -2_i4*m, in2)
        sum  = sum +tmp2
    end do 
    tmp0 = max(abs(tmp1), abs(tmp2))
    dsum = 1.d0 
    dm = 0 
    count = 0 
    do while(dsum > 1.d-16)
        count = count +1 
        dm = dm +1 
        m  = m0 -dm 
            tmp1 = fgsl_sf_bessel_jcn(m, in1)*fgsl_sf_bessel_jcn(n -2_i4*m, in2)
            sum  = sum +tmp1 
        m  = m1 +dm 
            tmp2 = fgsl_sf_bessel_jcn(m, in1)*fgsl_sf_bessel_jcn(n -2_i4*m, in2)
            sum  = sum +tmp2 
        tmpm = max(abs(tmp1), abs(tmp2))
        tmpr = tmpm/tmp0 
        if(tmpr < 1.d0 .and. tmpm < abs(sum)) then 
            dsum = tmpm/(1.d0 -tmpr)/abs(sum)
        end if 
        if(mod(count, 100_i4) == 0) then 
            write(*, *) "Warning #702: Roop cycle", count 
        end if 
        tmp0 = tmpm 
    end do 
    write(*, *) "Info #703: Do-while roop exits", count 
    transfer_C = sum 
end function transfer_C
! determinant matrix WD1 ---------------------------
function det_WD1(x)
    use linear, only: det_cmplx 
    complex(dp), intent(in)  :: x 
    complex(dp), allocatable :: work(:, :)
    complex(dp) :: det_WD1, tmp1, tmp2 
    integer(i4) :: i 
    if(allocated(work)) deallocate(work)
    allocate(work(-F:F, -F:F))
    work(:, :) = -W(:, :)
    tmp1 = x**2.d0 
    do i = -F, F 
        tmp2 = 2.d0*d(i)*x 
        work(i, i) = work(i, i) +tmp1 +tmp2 
    end do 
    det_WD1 = det_cmplx(work)
    if(allocated(work)) deallocate(work)
end function det_WD1
! diff det matrix WD1 ------------------------------
function ddet_WD1(x)
    use linear, only: det_cmplx 
    complex(dp), intent(in)  :: x 
    complex(dp), allocatable :: work(:, :)
    complex(dp) :: ddet_WD1, tmp1, tmp2, tmp3, tmp4 
    complex(qp) :: sum 
    integer(i4) :: i, j 
    if(allocated(work)) deallocate(work)
    allocate(work(-F:F, -F:F))
    tmp1 = x**2.d0 
    tmp3 = 2.d0*x 
    sum = 0.d0 
    do j = -F, F 
        work(:, :) = -W(:, :)
        tmp4 = 2.d0*d(j)
        do i = -F, F 
            tmp2 = 2.d0*d(i)*x 
            work(i, i) = work(i, i) +tmp1 +tmp2 
        end do 
        work(:, j) = 0.d0 
        work(j, j) = tmp3 +tmp4 
        tmp2 = det_cmplx(work)
        sum = sum +tmp2 
    end do 
    ddet_WD1 = sum 
    if(allocated(work)) deallocate(work)
end function ddet_WD1
! find p -------------------------------------------
function find_p(x)
    complex(dp), intent(in) :: x 
    real   (sp) :: sign
    real   (dp) :: delta
    complex(dp) :: find_p, x1, x2, fx, dfx, dx1, dx2 
    integer(i4) :: num2, num3, count
    x1    = x 
    dx1   = 0.d0 
    delta = 1.e0 
    num2  = 0 
    num3  = 0 
    count = 0 
    do while(delta > tiny1 .and. num2 < 5 .and. num3 < 20)
        count = count +1 
        fx    =  det_WD1(x1) 
        print *, x, abs(fx) 
        dfx   = ddet_WD1(x1) 
        dx2   = -fx/dfx
!         sign  = real(dx1*conjg(dx2))
!         if(sign < 0.e0) then 
!             dx2 = 0.5d0*dx2
!         end if 
        x2    = x1 +dx2 
        delta = abs(dx2/x1) 
        if(delta < tiny3) then 
            num3 = num3 +1 
            if(delta < tiny2) then 
                num2 = num2 +1
            end if 
        end if 
        if(mod(count, 100_i4) == 0) then 
            write(*, *) "Warning #190: Roop cycle", count 
        end if 
        x1    = x2 
        dx1   = dx2 
    end do 
    write(*, *) "Info #191: Do-while roop exits", count 
    find_p = x2 
end function find_p 
! outer region wave function -----------------------
function wave_G(r, j1, j2) 
    use math_const, only: i => math_i
    real   (dp), intent(in) :: r
    integer(i4), intent(in) :: j1, j2 
    complex(dp) :: wave_G, p 
    p = CH(j2, 2)
    wave_G = exp(i*p*r)*A(j1, j2)
end function wave_G
! outer region differcial wave function ------------
function wave_dG(r, j1, j2) 
    use math_const, only: i => math_i
    real   (dp), intent(in) :: r
    integer(i4), intent(in) :: j1, j2 
    complex(dp) :: wave_dG, p 
    p = CH(j2, 2)
    wave_dG = i*p*exp(i*p*r)*A(j1, j2)
end function wave_dG
! end function -------------------------------------










! ==================================================
! SUB-CALCULATION 
! ==================================================
! length gauge R -----------------------------------
subroutine SUB_lngth_R(l, energy)
    character(60), parameter :: form_out = '(1A15, 1ES15.3, 1ES15.3, "i")'
    integer(i4),  intent(in) :: l 
    complex(dp),  intent(in) :: energy
    integer(i4) :: i, j 

    lngth_R(l, :, :) = 0.d0 
    do i = -F, F 
        do j = -F, F 
            lngth_R(i, j, l) = cmplx_R(i, j, energy)
        end do 
    end do 
    write(file_log, form_out) "(length) R: ", lngth_R(0, 0, l)
end subroutine SUB_lngth_R
! velocity gauge R ---------------------------------
subroutine SUB_vlcty_R(l)
    use linear, only: inverse_cmplx
    character(60), parameter  :: form_out = '(1A15, 1ES15.3, 1ES15.3, "i")'
    integer  (i4), intent(in) :: l 
    complex(dp),  allocatable :: rR(:, :), Y(:, :)
    real   (dp),  allocatable :: C(:, :), dC(:, :)
    complex(qp) :: sum 
    integer(i4) :: i, j, k 

    vlcty_R(l, :, :) = 0.d0 
    if(allocated(rR)) deallocate(rR)
    if(allocated(Y))  deallocate(Y)
    if(allocated(C))  deallocate(C)
    if(allocated(dC)) deallocate(dC)
    allocate(rR(-F:F, -F:F))
    allocate( Y(-F:F, -F:F))
    allocate( C(-F:F, -F:F))
    allocate(dC(-F:F, -F:F))

    rR(:, :) = lngth_R(:, :, l)
    call inverse_cmplx(rR)
    do j = -F, F 
        do i = -F, F 
            C (i, j) = transfer_C(i -j, Ra)
            dC(i, j) = Amp/Freq/2.d0*(transfer_C(i -j -1_i4, Ra) -transfer_C(i -j +1_i4, Ra))
        end do 
    end do 
    do j = -F, F 
        do i = -F, F 
            sum = 0.d0 
            do k = -F, F 
                sum = sum +C(i, k)*rR(k, j)/Ra 
            end do 
            Y(i, j) = sum 
        end do 
    end do 
    do i = -F, F 
        Y(i, :) = Y(i, :) +dC(i, :)
    end do 
    call inverse_cmplx(Y)
    do i = -F, F 
        do j = -F, F 
            sum = 0.d0 
            do k = -F, F 
                sum = sum +C(i, k)*Y(k, j)/Ra 
            end do 
            vlcty_R(i, j, l) = sum
        end do 
    end do 

    if(allocated(rR)) deallocate(rR)
    if(allocated(Y))  deallocate(Y)
    if(allocated(C))  deallocate(C)
    if(allocated(dC)) deallocate(dC)
    write(file_log, form_out) "(velocity) R: ", vlcty_R(0, 0, l)
end subroutine SUB_vlcty_R
! channel p ----------------------------------------
subroutine SUB_channel(energy)
    use math_const, only: i => math_i 
    use linear,     only: diag_her 
    integer  (i1), parameter :: file_assumed = 101, file_channel = 102
    character(45), parameter :: & 
        form_out   = '(1A15, 1ES15.3, 1ES15.3, "i")', & 
        form_label = '(15X, 1A15, 1I15.4)', & 
        form_gen   = '(1I25, 2ES25.10)'
    complex(dp), intent(in)  :: energy
    complex(dp), allocatable :: ksq(:)
    complex(dp) :: tmp1, tmp2, x1, x2
    complex(qp) :: sum 
    real   (dp) :: dAmp, mAmp 
    integer(i4) :: j1, j2, k, n
    integer(i4), save :: count = 0 
    character(4) :: num 

    count = count +1 
    if(allocated(ksq)) deallocate(ksq)
    allocate(ksq(-F:F)) 
    dAmp = Amp/dble(cal)

    tmp1 = energy
    do j1 = -F, F 
        tmp2 = 2.d0*Mass*(tmp1 +dble(j1)*Freq)
        CH(-j1, 1) = -tmp2**0.5d0 
        CH( j1, 2) =  tmp2**0.5d0 
    end do 

    do n = 1, cal 
        mAmp = dAmp*dble(n)
        P = mAmp/Freq 
        U2(:, :) = 0.d0 
        do j2 = -F +1, F
            U2(j2 -1, j2) = P/2.d0/i 
        end do 
        call diag_her(U2, d) 

        tmp1 = energy -mAmp**2.d0/Freq**2.d0/4.d0
        do j1 = -F, F 
            tmp2 = 2.d0*Mass*(tmp1 +dble(j1)*Freq)
            ksq(j1) = tmp2 
        end do 
        do j2 = -F, F 
            do j1 = -F, F 
                sum = 0.d0 
                do k = -F, F
                    sum = sum +conjg(U2(k, j1))*ksq(k)*U2(k, j2) 
                end do 
                W(j1, j2) = sum 
            end do 
        end do 

        do j2 = 1, 2 
            do j1 = -F, F 
                x1 = CH(j1, j2) 
                x2 = find_p(x1)
                CH(j1, j2) = x2 
            end do 
        end do 
    end do 

    write(num, '(1I4.4)') count 
    open(file_channel, file = 'output/channel_'//num//'.d') 
    do j2 = 1, 2 
        do j1 = -F, F 
            write(file_channel, form_gen) j1, real(CH(j1, j2)), aimag(CH(j1, j2))
        end do 
    end do 
    close(file_channel)

    if(allocated(ksq)) deallocate(ksq)
    write(file_log, form_out)   "channel p: ", CH(0, 1)
    write(file_log, form_out)   "",            CH(0, 2)
    write(file_log, form_label) "file label ", count 
end subroutine SUB_channel
! propa R ------------------------------------------
! subroutine SUB_propa_R(l)
!     use math_const, only: i => math_i
!     use linear, only: inverse_cmplx
!     integer(i4), intent(in)  :: l 
!     complex(dp), allocatable :: Y1(:, :), Y2(:, :)
!     complex(dp) :: tmp1, tmp2, tmp3, tmp4 
!     complex(qp) :: sum 
!     integer(i4) :: j1, j2, k1, k2
!     if(allocated(Y1)) deallocate(Y1)
!     if(allocated(Y2)) deallocate(Y2)
!     allocate(Y1(-F:F, -F:F))
!     allocate(Y2(-F:F, -F:F))
!     do j2 = -F, F
!     tmp4 = exp(-i*d(j2)*Ra)
!         do j1 = -F, F
!             tmp1 = exp(i*d(j1)*Ra)
!             sum  = 0.d0 
!                 do k2 = -F, F 
!                 tmp3 = U(k2, j2)
!                     do k1 = -F, F 
!                         tmp2 = conjg(U(k1, j1))*vlcty_R(l, k1, k2)
!                         sum  = sum +tmp1*tmp2*tmp3*tmp4 
!                     end do 
!                 end do 
!             Y2(j1, j2) = sum 
!         end do 
!     end do 
!     do j2 = -F, F 
!         Y1(:, j2) = i*Ra*Y2(:, j2)*d(j2)
!     end do 
!     do j2 = -F, F 
!         Y1(j2, j2) = Y1(j2, j2) +1.d0 
!     end do 
!     call inverse_cmplx(Y1)
!     do j2 = -F, F 
!         do j1 = -F, F 
!             sum = 0.d0 
!             do k1 = -F, F 
!                 sum  = sum +Y1(j1, k1)*Y2(k1, j2)
!             end do 
!             modif_R(j1, j2, l) = sum 
!         end do 
!     end do 
!     if(allocated(U))  deallocate(U)
!     if(allocated(d))  deallocate(d)
!     if(allocated(Y1)) deallocate(Y1)
!     if(allocated(Y2)) deallocate(Y2)
! !     if(.not. allocated(dr)) allocate(dr(-F:F, -F:F))
! !     dr(:, :) = 0.d0 
! !     do j2 = -F, F 
! !         do j1 = -F, F 
! !             dr(j1, j2) = abs(modif_R(j1, j2, l))**2.d0 
! !         end do 
! !     end do 
! !     tmp = maxval(dr(:, :))
! !     tmp = minval(dr(:, :))
! !     dr(:, :) = 0.d0 
! !     do j2 = -F, F 
! !         do j1 = -F, F 
! !             dr(j1, j2) = abs(modif_R(j1, j2, l) -conjg(modif_R(j2, j1, l)))**2.d0/abs(modif_R(j1, j2, l))**2.d0 
! !         end do 
! !     end do 
! !     open(100, file = 'output/modif_r.d')
! !     call plot_mat(100, dr)
! !     if(allocated(dr)) deallocate(dr)
! end subroutine SUB_propa_R
! modif R ------------------------------------------
subroutine SUB_modif_R(l)
    use math_const, only: i => math_i
    use linear, only: diag_her, inverse_cmplx
    character(60), parameter  :: form_out = '(1A15, 1ES15.3, 1ES15.3, "i")'
    integer(i4),   intent(in) :: l 
    complex(dp),  allocatable :: R(:, :)
    complex(dp) :: tmp1, tmp2 
    complex(qp) :: sum 
    integer(i4) :: j1, j2, k1, k2 

    modif_R(:, :, l) = 0.d0 
    if(allocated(R)) deallocate(R)
    allocate(R(-F:F, -F:F))
    R(:, :) = vlcty_R(:, :, l)
    call inverse_cmplx(R)
    do j1 = -F +1, F 
        R(j1 -1, j1) = R(j1 -1, j1) +0.5d0*Ra*P 
        R(j1, j1 -1) = R(j1, j1 -1) -0.5d0*Ra*P 
    end do 
    call inverse_cmplx(R)
    do j2 = -F, F 
        do j1 = -F, F 
            sum = 0.d0 
            do k1 = -F, F 
                tmp1 = conjg(U2(k1, j1))
                do k2 = -F, F 
                    tmp2 = R(k1, k2)*U2(k2, j2) 
                    sum = sum +tmp1*tmp2
                end do 
            end do 
            modif_R(j1, j2, l) = sum 
        end do 
    end do 
    if(allocated(R)) deallocate(R)
    write(file_log, form_out) "(modified) R: ", modif_R(0, 0, l)
end subroutine SUB_modif_R
! matrix A -----------------------------------------
subroutine SUB_A
    use linear, only: solve_cmplx_vec 
    use plot,   only: plot_mat
    complex(dp), allocatable :: work(:, :), y(:)
    complex(dp) :: x, tmp1, tmp2 
    real   (qp) :: sum 
    integer(i4) :: i, j 

    if(allocated(work)) deallocate(work)
    if(allocated(y))    deallocate(y)
    allocate(work(-F:F, -F:F))
    allocate(y(-F:F))
    do j = -F, F 
        work(:, :) = -W(:, :)
        x = CH(j, 2)
        tmp1 = x**2.d0 
        do i = -F, F 
            tmp2 = 2.d0*d(i)*x 
            work(i, i) = work(i, i) +tmp1 +tmp2 
        end do 
        y(:) = -work(:, -F)
        y(-F) = 1.d0 
        call solve_cmplx_vec(work(-(F -1):F, -(F -1):F), y(-(F -1):F))
        sum = 0.d0 
        do i = -F, F 
            sum = sum +abs(y(i))**2.d0 
        end do 
!         tmp2 = y(0)
!         tmp2 = tmp2/abs(tmp2)
!         y(:) = y(:)/tmp2 
        y(:) = y(:)/sum**0.5d0
        A(:, j) = y(:)
    end do 
!     open(100, file = 'output/mat_a.d')
!     call plot_mat(100, abs(A(:, :))**2.d0)
! !     call plot_mat(100, aimag(A(:, :)))
!     close(100)
    if(allocated(work)) deallocate(work)
    if(allocated(y))    deallocate(y)
end subroutine SUB_A
! matching boundary conditon -----------------------
subroutine SUB_matching(l)
    use math_const, only: i => math_i
    use linear,     only: det_cmplx, solve_cmplx_vec
    character(60), parameter :: form_out = '(1A15, 1ES15.3, 1ES15.3, "i")'
    integer(i4),  intent(in) :: l 
    complex(dp), allocatable :: Gap1(:, :), Gap2(:, :), work(:, :), y(:)
    complex(dp) :: tmp 
    complex(qp) :: sum1
    real   (qp) :: sum2 
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
            Gap1(j1, j2) = wave_G(Rap, j1, j2)
            Gap2(j1, j2) = wave_dG(Rap, j1, j2)
        end do 
    end do 
    do j1 = -F, F 
        Gap2(j1, :) = Gap2(j1, :) +i*d(j1)*Gap1(j1, :)
    end do 
    do j2 = -F, F 
        do j1 = -F, F 
            sum1 = 0.d0 
            do k = -F, F 
                sum1 = sum1 +modif_R(j1, k, l)*Rap*Gap2(k, j2)
            end do 
            work(j1, j2) = sum1 
        end do 
    end do 
    work(:, :) = Gap1(:, :) -work(:, :)

    y(:) = -work(:, -F)
    y(-F) = 1.d0 
    call solve_cmplx_vec(work(-(F -1):F, -(F -1):F), y(-(F -1):F))
    sum2 = 0.d0 
    do j1 = -F, F 
        sum2 = sum2 +abs(y(j1))**2.d0 
    end do 
    y(:) = y(:)/sum2**0.5d0 
    sum1 = 0.d0 
    do j1 = -F, F 
        sum1 = sum1 +work(-F, j1)*y(j1)
    end do 

    write(*, *) sum1 
    tmp = det_cmplx(work)

    if(allocated(Gap1)) deallocate(Gap1)
    if(allocated(Gap2)) deallocate(Gap2)
    if(allocated(work)) deallocate(work)
    if(allocated(y))    deallocate(y)
    write(file_log, form_out) "matching: ", tmp 
end subroutine SUB_matching
! end sub-calculation ------------------------------


! ==================================================
! TEST 
! ==================================================










! ==================================================
! PROCESS
! ==================================================
! boundary matrix ----------------------------------
subroutine PROC_boundary_mat(input_m, input_l)
    use math_const,  only: i => math_i
    use hamiltonian, only: coord_E
    use linear,      only: diag_her
    integer(i4), intent(in) :: input_m, input_l 
    complex(dp) :: energy

    energy = coord_E(input_m) +gamma

    if(allocated(lngth_R)) deallocate(lngth_R)
    if(allocated(vlcty_R)) deallocate(vlcty_R)
    allocate(lngth_R(-F:F, -F:F, 0:L))
    allocate(vlcty_R(-F:F, -F:F, 0:L))
!     if(.not. allocated(lngth_R)) then 
!         if(op_mat_r == 1) then 
!             allocate(lngth_R(-F:F, -F:F, 0:L))
!         else if(op_mat_r == 0) then 
!             allocate(lngth_R(-F:F, -F:F, input_l:input_l))
!         end if 
!     end if 
    call SUB_lngth_R(input_l, energy)
    call SUB_vlcty_R(input_l)

!     if(allocated(U1)) deallocate(U1)
    if(allocated(U2)) deallocate(U2)
    if(allocated(W))  deallocate(W)
    if(allocated(d))  deallocate(d)
    if(allocated(CH)) deallocate(CH)
!     allocate(U1(-F:F, -F:F))
    allocate(U2(-F:F, -F:F))
    allocate( W(-F:F, -F:F))
    allocate( d(-F:F))
    allocate(CH(-F:F, 1:2))
    call SUB_channel(energy)

!     if(allocated(propa_R)) deallocate(propa_R)
    if(allocated(modif_R)) deallocate(modif_R)
!     allocate(propa_R(-F:F, -F:F, 0:L))
    allocate(modif_R(-F:F, -F:F, 0:L))
!     call SUB_propa_R(input_l)
    call SUB_modif_R(input_l)

    if(allocated(A)) deallocate(A)
    allocate(A(-F:F, -F:F))
    call SUB_A
!     if(allocated(U1)) deallocate(U1)
    if(allocated(U2)) deallocate(U2)
    if(allocated(W))  deallocate(W)

    call SUB_matching(input_l)
    if(allocated(d))  deallocate(d)
    if(allocated(CH)) deallocate(CH)

!     if(.not. allocated(K)) then 
!         if(op_mat_k == 1) then 
!             allocate(K(0:L, -F:F, -F:F)) 
!         else if(op_mat_k == 0) then 
!             allocate(K(input_l:input_l, -F:F, -F:F))
!         end if 
!     end if 
!     if(F0 <= F) then 
!         call mat_K(input_l)
!     else if(F0 > F) then 
!         K(input_l, :, :) = 0.d0 
!     end if 
!     if(.not. allocated(S)) then 
!         if(op_mat_s == 1) then 
!             allocate(S(1:M, 0:L, -F:F)) 
!         else if(op_mat_s == 0) then 
!             allocate(S(input_m:input_m, input_l:input_l, -F:F))
!         end if 
!     end if 
!     if(F0 <= F) then 
!         call SUB_S(input_m, input_l)
!     else if(F0 > F) then 
!         S(input_m, input_l, :) = 0.d0 
!     end if 
!     if(op_mat_a == 1) then 
!         if(.not. allocated(A)) then 
!             if(op_mat_a == 1) then 
!                 allocate(A(0:L)) 
!             else if(op_mat_a == 0) then 
!                 allocate(A(input_l:input_l))
!             end if 
!         end if 
!         call mat_A(input_m, input_l) ! use only S matrix 
!     end if 

    if(op_mat_r == 0 .and. allocated(lngth_R)) deallocate(lngth_R)
    if(op_mat_r == 0 .and. allocated(vlcty_R)) deallocate(vlcty_R)
!     if(op_mat_r == 0 .and. allocated(propa_R)) deallocate(propa_R)
    if(op_mat_r == 0 .and. allocated(modif_R)) deallocate(modif_R)
!     if(op_mat_k == 0 .and. allocated(K)) deallocate(K)
!     if(op_mat_s == 0 .and. allocated(S)) deallocate(S)
!     if(op_mat_a == 0 .and. allocated(A)) deallocate(A)
end subroutine PROC_boundary_mat
! end boundary matrix ------------------------------
! bound states -------------------------------------
subroutine PROC_bound(l)
!     use fgsl, only: fgsl_sf_bessel_ksl_scaled
!     integer  (i1), parameter  :: file_bound = 101
!     character(30), parameter  :: form_bound = '(2ES25.10)', form_out = '(1A15, 1I15, 2ES15.3)'
    integer  (i4), intent(in) :: l 
!     real     (dp) :: dE, E0, Ei, tmp1, tmp2, d1, d2 
!     real     (dp) :: ka, sb_ka, diff_ka 
!     integer  (i4) :: i, j, m, num 
!     character (3) :: ch 

!     write(ch, '(I3.3)') l 
!     open(file_bound, file = "output/bound_"//ch//".d")
!     m   = 1000
!     num = 0 
!     do j = 1, N 
!         tmp1 = E(j)
!         tmp2 = E(j +1)
!         if(.not. tmp1 < 0.d0) exit 
!         if(.not. tmp2 < 0.d0) tmp2 = 0.d0 
!         E0  = tmp1
!         dE  = (tmp2 -tmp1)/dble(m)
!         d1  = 0.d0 
!         do i = 1, m -1 
!             Ei = E0 +dble(i)*dE 
!             ka = (2.d0*Mass*(-Ei))**0.5d0*Ra
!             sb_ka = fgsl_sf_bessel_ksl_scaled(l, ka)/exp(ka)*ka 
!             if(.not. l == 0) then 
!                 tmp1 = dble(l)*fgsl_sf_bessel_ksl_scaled(l, ka)/exp(ka)
!                 tmp2 = dble(l +1)*fgsl_sf_bessel_ksl_scaled(l +1_i4, ka)/exp(ka)
!             else if(l == 0) then 
!                 tmp1 = 0.d0 
!                 tmp2 = dble(l +1)*fgsl_sf_bessel_ksl_scaled(l +1_i4, ka)/exp(ka)
!             end if 
!             diff_ka = -ka**2.d0/dble(2*l +1)*(tmp1 +tmp2)
!             d2 = sb_ka -func_R(Ei)*(sb_ka +diff_ka)
!             write(file_bound, form_bound) Ei, d2
!             if(d1*d2 < 0.d0) then 
!                 num  = num +1 
!                 tmp1 = -dE*abs(d2)/(abs(d1) +abs(d2))
!                 tmp2 = Ei +tmp1 
!                 write(file_log, form_out) "Bound State: ", num, tmp2
!             end if 
!             d1 = d2 
!         end do 
!     end do 
!     close(file_bound)
end subroutine PROC_bound
! end bound states ---------------------------------
! out ----------------------------------------------
subroutine PROC_boundary_out
    if(allocated(lngth_R)) deallocate(lngth_R)
    if(allocated(vlcty_R)) deallocate(vlcty_R)
    if(allocated(A)) deallocate(A)
    if(allocated(S)) deallocate(S)
end subroutine PROC_boundary_out
! end out ------------------------------------------
end module boundary
