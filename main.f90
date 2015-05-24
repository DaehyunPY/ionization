program main
<<<<<<< HEAD
<<<<<<< HEAD
    use kind_type
    use global
    use hamiltonian, only: PROC_input, PROC_inform, PORC_coord, PROC_prog_out 
    use basis,       only: PROC_H, PROC_basis_plot, PROC_bound
    use boundary,    only: PROC_boundary_mat, PROC_mat_out
    use inner,       only: PROC_inner_achive, PROC_inner_plot
    use outer,       only: PROC_outer_plot
    use asymptote,   only: PROC_CS_plot, PROC_E_vs_CS_plot, PROC_E_vs_PS_plot
    implicit none
    character(30), parameter :: & 
        form_step1 = '(1A20, 1I5)', &
        form_step2 = '(1A40, 1I8, 1A4, 1I8)', &
        form_time1 = '(1A25, 5X, 1ES15.3)', &
        form_time2 = '(1A15, 1ES15.3)'
    real     (sp) :: tt, t0, t1, t2 
    integer  (i4) :: i, j 

    call cpu_time(tt)
    call cpu_time(t1)
    write(*, *) "Reading input file..."
    call PROC_input
    write(file_log, *) "[PROCESS INPUT]"
    write(file_log, *)
        call PROC_inform
        call PORC_coord
        write(file_log, *)
    call cpu_time(t2)
    write(file_log, form_time1) "PROCESS RUNNING TIME: ", t2 -t1 
    write(file_log, *) 
    write(file_log, *) 
    write(file_log, *) 
    write(file_log, *) 
    write(file_log, *) 
    write(*, *) "Reading over."

    call cpu_time(t0)
    write(*, *) "Calculating..."
    write(file_log, *) "[PROCESS CALCULATE]"
    write(file_log, *) 
    do i = 0, L 
        call cpu_time(t1)
        write(file_log, form_step1) "ANGULAR MOMANTUM ", i
        call PROC_H(i) 
        if(op_basis == "Y") then 
            call PROC_basis_plot(i)
        end if 
        if(op_bound == "Y") then 
            call PROC_bound(i)
        end if 
        call cpu_time(t2)
        write(file_log, form_time2) "RUNNING TIME: ", t2 -t1
        write(file_log, *)
        do j = 1, M 
            call cpu_time(t1)
            write(file_log, form_step2) "STEP ", j, "of", M 
            if(op_dcs == "Y" &
                    .or. (op_inner == "Y" .or. op_outer == "Y") &
                    .or. (op_tcs   == "Y" .or. op_ps    == "Y")) then 
                call PROC_boundary_mat(j, i) 
            end if 
            if(op_inner == "Y" .and. M == 1) then 
                call PROC_inner_achive(i)
            end if 
            call cpu_time(t2) 
            write(file_log, form_time2) "RUNNING TIME: ", t2 -t1
            write(file_log, *) 
        end do 
    end do 
    call PROC_mat_out
    call cpu_time(t2) 
    write(file_log, form_time1) "PROCESS RUNNING TIME: ", t2 -t0 
    write(file_log, *) 
    write(file_log, *) 
    write(file_log, *) 
    write(file_log, *) 
    write(file_log, *) 
    write(*, *) "Calculating over."

    call cpu_time(t1)
    write(*, *) "Ploting..."
    write(file_log, *) "[PROCESS PLOT]"
    write(file_log, *) 
        if(op_inner == "Y") call PROC_inner_plot
        if(op_outer == "Y") call PROC_outer_plot 
        if(op_dcs   == "Y") call PROC_CS_plot 
        if(op_tcs   == "Y") call PROC_E_vs_CS_plot
        if(op_ps    == "Y") call PROC_E_vs_PS_plot
        write(file_log, *) 
    call cpu_time(t2)
    write(file_log, form_time1) "PROCESS RUNNING TIME: ", t2 -t1 
    write(file_log, *)
    write(file_log, *)
    write(file_log, *)
    write(file_log, *)
    write(file_log, *)
    write(*, *) "Ploting over."

    write(file_log, *) "PROGRAM OVER"
    write(file_log, *)
    call cpu_time(t2)
    write(file_log, form_time1) "PROGRAM RUNNING TIME: ", t2 -tt
    write(file_log, *)
    call PROC_prog_out 
    write(*, *) "Program over."
=======
    use global
    use hamiltonian, only: process_coord
    use inner, only: process_basis_H1, process_basis_HF
    use outer, only: process_find_ground, process_find_bound
    implicit none
=======
    use global
    use hamiltonian, only: process_coord
    use inner, only: process_basis_H1, process_basis_HF
    use outer, only: process_find_ground, process_find_bound
    implicit none
>>>>>>> 9f36052... new code
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

<<<<<<< HEAD
>>>>>>> 9f36052... new code
=======
>>>>>>> 9f36052... new code
end program main


















































