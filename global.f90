module global 
    use kind_type 
    implicit none
    integer (i1), parameter :: file_input = 5, file_log = 6 
    real    (dp), save :: Amp, Freq, Mass, Charge, Spin, Scatt, Bound, dr_p_drho
    integer (i4), save :: L, N, F, M, pr, ptheta 
    real    (dp), save, allocatable :: &
        coord_rho(:), coord_weight(:), coord_dshape(:, :), &
        H(:, :), E(:), R(:), K(:)
    complex (dp), save, allocatable :: S(:, :), A(:)
    character(1), save :: &
        op_ev, op_degree, op_aa, & 
        op_basis, op_bound, & 
        op_dcs, op_inner, op_outer, &
        op_tcs, op_ps
end module
