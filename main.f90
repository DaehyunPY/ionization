program main
    use kind_type
    use global
    use hamiltonian, only: coord_E
    use hamiltonian, only: PROC_input, PROC_inform, PORC_coord, PROC_Poten_plot, PROC_out 
    use basis,       only: PROC_H, PROC_basis_plot
    use boundary,    only: PROC_boundary_mat
    use inner,       only: PROC_inner_achive, PROC_inner_plot
    use outer,       only: PROC_CS_achive, PROC_CS_plot, PROC_outer_plot, PROC_E_vs_CS_plot
    implicit none
    character(30), parameter :: form_out = '(1A25, 5X, 1ES15.3)'
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
        if(op_poten == "Y") then
            call PROC_Poten_plot
            write(file_log, *) "Potential function is ploted."
        end if
        write(file_log, *)
    call cpu_time(t2)
    write(file_log, form_out) "PROCESS RUNNING TIME: ", t2 -t1 
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
    do j = 1, M 
        Scatt = coord_E(j)
        do i = 0, L 
            call cpu_time(t1)
            write(file_log, *) "ROUND", i, coord_E(j)
                call PROC_H(i) 
                if(op_basis == "Y") &
                call PROC_basis_plot(i)
                if((op_dcs == "Y" .or. op_inner == "Y" .or. op_outer == "Y") .or. &
                    (op_tcs == "Y" .or. op_phase == "Y" .or. op_lt == "Y")) & 
                call PROC_boundary_mat(i) 
                if(op_inner == "Y") & 
                call PROC_inner_achive(i)
            call cpu_time(t2) 
            write(file_log, form_out) "RUNNING TIME: ", t2 -t1
            write(file_log, *) 
        end do 
        if(op_basis == "Y") & 
        write(file_log, *) "Basis function is ploted."
        if(op_tcs == "Y" .or. op_phase == "Y" .or. op_lt == "Y") & 
        call PROC_CS_achive(j)
        write(file_log, *)
    end do 
    call cpu_time(t2) 
    write(file_log, form_out) "PROCESS RUNNING TIME: ", t2 -t0 
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
        if(op_dcs   == "Y") then 
            call PROC_CS_plot 
            write(file_log, *) "Cross Section function is ploted."
        end if 
        if(op_inner == "Y") then 
            call PROC_inner_plot
            write(file_log, *) "Inner region function is ploted."
        end if 
        if(op_outer == "Y") then 
            call PROC_outer_plot 
            write(file_log, *) "Outer region function is ploted."
        end if 
        if(op_tcs   == "Y") then 
            call PROC_E_vs_CS_plot
            write(file_log, *) "Energy vs Cross Section function is ploted."
        end if
!         if(op_phase == "Y") 
!         if(op_lt    == "Y")
        write(file_log, *)
    call cpu_time(t2)
    write(file_log, form_out) "PROCESS RUNNING TIME: ", t2 -t1 
    write(file_log, *)
    write(file_log, *)
    write(file_log, *)
    write(file_log, *)
    write(file_log, *)
    write(*, *) "Ploting over."

    write(file_log, *) "PROGRAM OVER"
    write(file_log, *)
    call cpu_time(t2)
    write(file_log, form_out) "PROGRAM RUNNING TIME: ", tt -t1 
    write(file_log, *)
    call PROC_out 
    write(*, *) "Program over."
end program main
