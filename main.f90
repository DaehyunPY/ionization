program main
    use global
    use hamiltonian, only: process_coord
    use inner, only: process_basis_H1, process_basis_HF
    use outer, only: process_find_ground, process_find_bound
    implicit none
    integer :: m
    integer, pointer :: N

    N => amp_N
    ! note: 1.1 read input-file
    call process_input
    print *, "Reading input file has finished. "
    ! note: 1.2 inform data
    ! note: 1.3 make coordinates & hamiltonian
    call process_coord
    print *, "Building coordinates has finished. "

    ! note: 2.1 no-laser calculation
    ! note: 2.1.1 inner-region basis
    call process_basis_H1
    print *, "Calculating basis of H1 has finished. "
    ! note: 2.1.2 find first bound-state
    call process_find_ground
    print *, "Finding ground state has finished. "
    ! note: 2.1.3 plot first bound-state function

    ! note: 2.2 laser calculation
    do m = 0, N
        print *, m, N
        ! note: 2.2.1 inner-region basis
        call process_basis_HF(m)
        print *, "Calculating basis of HF has finished. "
        ! note: 2.1.2 find complex-energy
        call process_find_bound
        print *, "Finding bound state has finished. "
    end do
    print *, "Program has finished. "

contains









! ==============================================================================
! PROCESS: process_input
! ==============================================================================
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    subroutine process_input ! :::::::::::::::::::::::::::::::::::::::::::::::::
        ! It reads a input file and check if the values of the file are valid.
        use hamiltonian, only: freq, amp_max, amp_min, poten_type, r_inner
        character(*), parameter :: check = "Check your input value: "
        character(1) :: dummy
        logical :: op_input, op_log

        inquire(file="input.data", exist = op_input)
        if(op_input == .true.) then
            open(file_input, file = "input.data")
        end if

        read(file_input, *) ! ===
        read(file_input, *) ! laser
        read(file_input, *) ! ===
        read(file_input, *) ! ---
        read(file_input, *) dummy, freq
        read(file_input, *) dummy, amp_max
        read(file_input, *) dummy, amp_min
        read(file_input, *) ! -
        read(file_input, *) ! -
        if(.not. freq > 0.d0) stop check//"FREQUENCY. "
        if(.not. amp_max > 0.d0) stop check//"MAXIMUM AMPLITUDE. "
        if(.not. (amp_min >= 0.d0 .and. amp_min < amp_max)) stop check//"MINIMUM AMPLITUDE. "

        read(file_input, *) ! ===
        read(file_input, *) ! potential
        read(file_input, *) ! ===
        read(file_input, *) ! ---
        read(file_input, *) dummy, poten_type
        read(file_input, *) ! -
        read(file_input, *) ! -
        if(.not. poten_type == 0) stop check//"POTENTIAL TYPE. "

        read(file_input, *) ! ===
        read(file_input, *) ! calculation
        read(file_input, *) ! ===
        read(file_input, *) ! ---
        read(file_input, *) dummy, r_inner
        read(file_input, *) dummy, r_N
        read(file_input, *) dummy, amp_N
        read(file_input, *) dummy, floq_N
        read(file_input, *) ! -
        read(file_input, *) ! -
        if(.not. r_inner > 0.d0) stop check//"INNER REGION SIZE. "
        if(.not. r_N > 0) stop check//"GRID NUMBER OF r COORDINATES. "
        if(.not. amp_N > 0) stop check//"GRID NUMBER OF Amplitude COORDINATES. "
        if(.not. floq_N >= 0) stop check//"NUMBER OF FLOQUET BLOCK. "

        read(file_input, *) ! ===
        read(file_input, *) ! option
        read(file_input, *) ! ===
        read(file_input, *) ! ---
        read(file_input, *) ! unit
        read(file_input, *) !
        read(file_input, *) ! output
        read(file_input, *) dummy, op_log
        read(file_input, *) ! -
        read(file_input, *) ! -

        if(op_input == .true.) then
            close(file_input)
        end if
        if(op_log == .true.) then
            open(file_log, file = "log.data")
        end if
    end subroutine process_input ! :::::::::::::::::::::::::::::::::::::::::::::
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! END PROCESS ------------------------------------------------------------------

end program main


















































