program main
    use kind_type
    use global
    use hamiltonian, only: PROC_input, PROC_inform, PORC_coord, PROC_hamiltonian_out
    use basis,       only: PROC_H, PROC_basis_plot, PROC_basis_break, PROC_basis_out
    use boundary,    only: PROC_bound, PROC_boundary_mat, PROC_boundary_out
!     use inner,       only: PROC_inner_achive, PROC_inner_plot
!     use outer,       only: PROC_outer_plot
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
        if(op_basis == 1) then 
            call PROC_basis_plot(i)
        end if 
        if(op_bound == 1) then 
            call PROC_bound(i)
        end if 
        call cpu_time(t2)
        write(file_log, form_time2) "RUNNING TIME: ", t2 -t1
        write(file_log, *)
        do j = 1, M 
            call cpu_time(t1)
            write(file_log, form_step2) "STEP ", j, "of", M 
            if(op_dcs == 1 &
                    .or. (op_inner == 1 .or. op_outer == 1) &
                    .or. (op_tcs   == 1 .or. op_ps    == 1)) then 
                call PROC_boundary_mat(j, i) 
            end if 
            if(op_inner == 1 .and. M == 1) then 
!                 call PROC_inner_achive(i)
            end if 
            call cpu_time(t2) 
            write(file_log, form_time2) "RUNNING TIME: ", t2 -t1
            write(file_log, *) 
        end do 
    end do 
    call PROC_basis_break
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
!         if(op_inner == 1) call PROC_inner_plot
!         if(op_outer == 1) call PROC_outer_plot 
        if(op_dcs   == 1) call PROC_CS_plot 
        if(op_tcs   == 1) call PROC_E_vs_CS_plot
        if(op_ps    == 1) call PROC_E_vs_PS_plot
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
    call PROC_boundary_out
    call PROC_basis_out
    call PROC_hamiltonian_out
    write(*, *) "Program over."
end program main
