

MODULE COULOMB_EXPLOSION
    USE MOD_GLOBAL  ! Global structures to be used throughout the program.
    USE COMMANDS    ! Commands to run (i.e. mkdir)
    USE MOD_IO      ! Reads from control.inp and molecule.inp
    USE MOLDYNAMICS ! Movement/dynamics/interactions of atoms/molecules
    USE BOLTZMANN   ! Random initial ion-velocity

    implicit none

CONTAINS

SUBROUTINE initialize
  CALL cpu_time(program_start_time)
  write(*,*) "-=================================-"
  write(*,*) "-=# CLASSICAL COULOMB EXPLOSION #=-"
  write(*,*) "-=================================-"
  write(*,*) "Run Started. Initializing..."

  call prepare_output
  call read_molecule_input_file('molecule.inp')
  call read_control_input_file('control.inp')
  call read_seeds_input_file('seeds.inp')
  call compute_atomic_masses

  write(log_file,*)""
  write(log_file,*)"Initialization complete. Beginning computation..."
  write(log_file,*)""

  write(*,*)"Initialization complete. Beginning computation..."

END SUBROUTINE initialize


SUBROUTINE prepare_output
  character(len=15) :: formatted_datetime
  logical :: dir_exists
  character(len=255) :: output_dir_with_date_time
  character(len=255) :: trajectory_directory

  call get_formatted_datetime(formatted_datetime)
  
  call check_and_create_directory(output_dir, dir_exists)
  if (.not. dir_exists) then
    write(*,'(A, A, A, A)') "Output directory created: ", "'", trim(output_dir), "'"
  else
      write(*,'(A, A, A, A)') "'", trim(output_dir), "'", " directory found"
  end if

  output_dir_with_date_time = trim(adjustl(output_dir))//"/"//formatted_datetime

  call check_and_create_directory(output_dir_with_date_time, dir_exists)
  write(*,'(A)') " Output folder format:              DDMMYYYY_HHMMSS"
  write(*,'(A, A, A, A)') " Run output folder created: ", "'", trim(output_dir_with_date_time), "'"


  trajectory_directory = trim(adjustl(output_dir_with_date_time))//"/"//"trajectory"

  call check_and_create_directory(trajectory_directory, dir_exists)
  write(*,'(A, A, A, A)') " Trajectory output folder created: ", "'", trim(trajectory_directory), "'"

  ! Adds the paths to each of the file_names
  call append_paths_to_filenames(formatted_datetime, trajectory_directory)
  
  ! Open output files
  open(log_file, file=trim(adjustl(log_output_filename)))
  open(all_variable_file, file=trim(adjustl(all_variable_filename)))
  open(atom_info_file, file=trim(adjustl(atom_info_filename)))
  ! open the trajectory file for each run-->open in the "calculate" subroutine

  write(log_file,*) "Run Started and output directory successfully created. Initializing..."
  write(log_file,'(A, A)') " Output folder created: ",  trim(output_dir_with_date_time)
  write(log_file,'(A, A)') " Output files created: ",  trim(output_dir_with_date_time)
  write(log_file,'(A, A, A, A)') " all variable file created: ", "'", trim(all_variable_filename), "'"
  write(log_File,'(A, A, A, A)') " atom info file created: ", "'", trim(atom_info_filename), "'"
  write(log_file,'(A, A, A, A)') " Trajectory output folder created: ", "'", trim(trajectory_directory), "'"

END SUBROUTINE prepare_output


! Calculates the Coulomb Force at each time step in the simulation
! and propagates it throughout time. Sets a different seed value for each simulation
SUBROUTINE calculate
  ! Initialize the atom's velocity via Boltzmann distribution
  integer :: i, seed, unit_num
  character(len=255) :: full_trajectory_filename, seed_string

  do i=1, N_simulations
    ! reset the seed and perform the computation again
    seed = seed_array(i)
    ion_velocity_init_seed = seed

    write(seed_string, '(G0)') seed
    full_trajectory_filename = trim(adjustl(trajectory_filename)) // trim(adjustl(seed_string)) // ".xyz"
    unit_num = trajectory_file+1 
    open(unit_num,file=trim(adjustl(full_trajectory_filename)))

    write(log_file,'(A, A, A)') "r=", trim(adjustl(seed_string)), " computation started"

    ! Initial state of each simulations
    atom_position = atom_initial_position
    call calculate_atomic_velocities
    call calculate_force

    do iter=0, N_time_steps
        time = iter * time_step
    
        ! Output
        if (mod(iter,trajectory_output_frequency) == 0) then
          call update_trajectory_file(unit_num)
        end if

        ! Update via Verlet Algorithm (force, position, velocity)
        call calculate_force
        call calculate_position
        call calculate_velocity
    end do

    write(log_file,*) "Successfully finished trajectory file: ", full_trajectory_filename
    call update_atom_info_file
    close(unit_num)
    write(log_file,'(A, A, A)') "r=", trim(adjustl(seed_string)), " computation completed"
    write(log_file,*)
    write(*,'(A, A, A)') "  r=", trim(adjustl(seed_string)), " computation completed"

  end do

END SUBROUTINE calculate


SUBROUTINE cleanup
  real*8 :: program_elapsed_time

  call cpu_time(program_end_time)
  program_elapsed_time = program_end_time - program_start_time

  write(log_file,*)"Run finished. Elapsed time (seconds): ", program_elapsed_time
  write(*,*)"Run finished. Elapsed time (seconds): ", program_elapsed_time

  ! close the output files
  close(log_file)
  close(trajectory_file)
  close(atom_info_file)
  ! deallocate all of the matrices
  deallocate(atom_initial_position)
  deallocate(atom_position)
  deallocate(atom_velocity)
  deallocate(atom_acceleration)
  deallocate(atom_force)
  ! deallocate all of the arrays
  deallocate(atom_atomic_number)
  deallocate(atom_charge)
  deallocate(atom_mass)
  deallocate(seed_array)

END SUBROUTINE cleanup


END MODULE COULOMB_EXPLOSION






