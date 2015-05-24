module global
    use kind_const
    implicit none
    ! NOTE: global > file: file_input, file_log
    integer, parameter :: file_input = 5, file_log = 6
    ! NOTE: global > calculation: loopmax, r_N, amp_N, floq_N
    integer, parameter :: loop_max = 1000
    integer, target :: r_N, amp_N, floq_N
end module global
