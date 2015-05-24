module global 
    use kind_type 
    implicit none
    integer (i1), parameter :: file_input = 5, file_log = 6 
    real    (dp), save :: Amp, Freq, Mass, Charge, Spin, Scatt, Ra, Rap
    integer (i4), save :: L, N, F, M, pr, ptheta 
    integer (i1), save :: & 
        op_mat_f = 0, op_mat_h = 0, op_mat_e = 0, &
        op_mat_r = 0, op_mat_k = 0, op_mat_s = 0, op_mat_a = 0, & 
        op_ev, op_degree, op_aa, & 
        op_basis, op_bound, & 
        op_dcs, op_inner, op_outer, &
        op_tcs, op_ps
end module
