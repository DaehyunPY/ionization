program main
    use kind_const
    use global
    use hamiltonian, only: PROC_input, PROC_inform, PORC_coord, PROC_hamiltonian_out
    use basis,       only: PROC_H, PROC_basis_break, PROC_basis_out
    use boundary,    only: PROC_matching, PROC_boundary_out
!     use inner,       only: PROC_inner_achive, PROC_inner_plot
!     use outer,       only: PROC_outer_plot
!     use asymptote,   only: PROC_CS_plot, PROC_E_vs_CS_plot, PROC_E_vs_PS_plot
    implicit none
    character(30), parameter :: &
        form_step1 = "(1A20, 1I5)", &
        form_step2 = "(1A40, 1I8, 1A4, 1I8)"
    character(60), parameter :: &
        form_time1 = "(1A25, 5X, 1I2.2, ' hr ', 1I2.2, ' min ', 1I2.2, ' sec ')", &
        form_time2 = "(1A15,     1I2.2, ' hr ', 1I2.2, ' min ', 1I2.2, ' sec ')"
    real     (sp) :: tt, t0, t1, t2
    integer  (i4) :: i, j, time, hr, min, sec


! ==================================================
! INPUT
! ==================================================
    call cpu_time(tt)
    call cpu_time(t1)
    write(*, *) "Reading input file..."

    ! note: 1.1. input parameters
    call PROC_input
    write(file_log, *) "[PROCESS INPUT]"
    write(file_log, *)

        ! note: 1.2. inform input parameters
        call PROC_inform

        ! note: 1.3. make coordination system
        call PORC_coord
        write(file_log, *)
    call cpu_time(t2)
    time = t2 -t1
    hr = time/3600
    min = (time -hr*3600)/60
    sec = (time -hr*3600 -min*60)
    write(file_log, form_time1) "PROCESS RUNNING TIME: ", hr, min, sec
    write(file_log, *)
    write(file_log, *)
    write(file_log, *)
    write(file_log, *)
    write(file_log, *)
    write(*, *) "Reading over."


! ==================================================
! CALCULATE
! ==================================================
    call cpu_time(t0)
    write(*, *) "Calculating..."
    write(file_log, *) "[PROCESS CALCULATE]"
    write(file_log, *)

    ! note: 2.1. angular momantum loop
    do i = 0, L
        write(file_log, form_step1) "ANGULAR MOMANTUM ", i
        write(file_log, *)

        ! note: 2.2. laser amplitude loop
        do j = 0, M
            call cpu_time(t1)
            write(file_log, form_step2) "STEP ", j, "of", M

            ! note: 2.3. calculate inner region basis
            call PROC_H(i, j)

            ! note: 2.4. matching boundary conditions
            call PROC_matching(j, i)
            call cpu_time(t2)
            time = t2 -t1
            hr = time/3600
            min = (time -hr*3600)/60
            sec = (time -hr*3600 -min*60)
            write(file_log, form_time2) "RUNNING TIME: ", hr, min, sec
            write(file_log, *)
        end do
    end do
    call cpu_time(t2)
    time = t2 -t0
    hr = time/3600
    min = (time -hr*3600)/60
    sec = (time -hr*3600 -min*60)
    write(file_log, form_time1) "PROCESS RUNNING TIME: ", hr, min, sec
    write(file_log, *)
    write(file_log, *)
    write(file_log, *)
    write(file_log, *)
    write(file_log, *)
    write(*, *) "Calculating over."


! ==================================================
! PLOT
! ==================================================
    call cpu_time(t1)
    write(*, *) "Ploting..."
    write(file_log, *) "[PROCESS PLOT]"
    write(file_log, *)

        ! NOTE: 3.1. wave function

        ! NOTE: 3.2. cross section
        write(file_log, *)
    call cpu_time(t2)
    time = t2 -t1
    hr = time/3600
    min = (time -hr*3600)/60
    sec = (time -hr*3600 -min*60)
    write(file_log, form_time1) "PROCESS RUNNING TIME: ", hr, min, sec
    write(file_log, *)
    write(file_log, *)
    write(file_log, *)
    write(file_log, *)
    write(file_log, *)
    write(*, *) "Ploting over."


! ==================================================
! OVER
! ==================================================
    write(file_log, *) "PROGRAM OVER"
    write(file_log, *)
    call cpu_time(t2)
    time = t2 -tt
    hr = time/3600
    min = (time -hr*3600)/60
    sec = (time -hr*3600 -min*60)
    write(file_log, form_time1) "PROGRAM RUNNING TIME: ", hr, min, sec
    write(file_log, *)
    call PROC_boundary_out
    call PROC_basis_out
    call PROC_hamiltonian_out
    write(*, *) "Program over."
end program main
