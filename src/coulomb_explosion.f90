! Masterfile that performs to Coulomb explosion computations

MODULE COULOMB_EXPLOSION
    USE MOD_GLOBAL  ! Global structures to be used throughout the program.
    USE COMMANDS    ! Commands to run (i.e. mkdir)
    USE MOD_IO      ! Reads from control.inp and molecule.inp
    USE MOLDYNAMICS ! Movement/dynamics/interactions of atoms/molecules
    USE BOLTZMANN   ! Random initial ion-velocity

    implicit none

CONTAINS

SUBROUTINE initialize
  character(len=3) :: run_type_str
  character(len=256) :: full_run_trajectory_input_file

  call date_and_time(date=start_date, time=start_time)
  call cpu_time(program_CPU_start_time)

  write(*,*) "-=================================-"
  write(*,*) "-=# CLASSICAL COULOMB EXPLOSION #=-"
  write(*,*) "-=================================-"
  write(*,*) "Run Started. Initializing..."

  call prepare_output
  call read_control_input_file('control.inp')
  write(run_type_str, '(I0)') run_type

  if (run_type == 1) then
    print*, "run_type = ", trim(adjustl(run_type_str)), " | mode set to: BOLTZMANN VELOCITY DISTRIBUTION FROM GROUNDSTATE"
    molecule_filename_full = trim(adjustl(molecule_input_path)) // 'molecule.inp'
    call read_molecule_input_file(molecule_filename_full)
    print*, "molecule file name = ", molecule_filename_full
    seeds_filename_full = trim(adjustl(seeds_input_path)) // 'seeds.inp'
    call read_seeds_input_file(seeds_filename_full)
    print*, "seeds file name = ", seeds_filename_full
  else if (run_type == 2) then
    print*, "run_type = ", trim(adjustl(run_type_str)), " | mode set to: CONTINUE FROM TDDFT WITHOUT QUANTUM EFFECTS"
    moleculeformations_filename_full = trim(adjustl(moleculeformations_input_path)) // 'moleculeFormations.csv'
    call read_molformations_input_file(moleculeformations_filename_full)
    print*, "Successfully read from: moleculeFormations file: ", trim(adjustl(moleculeformations_filename_full))
    ! READ THE NUMBER OF ATOMS FROM THE FIRST TRAJECTORY FILE. ALSO READ IN THE SPECIES. NEED THIS INFO TO COMPUTE MASS
    full_run_trajectory_input_file = trim(adjustl(full_runs_dir_input_path)) // &
                                     trim(adjustl(full_runs_array(1))) // "/trajectory.xyz"
    call read_trajectory_full_run(full_run_trajectory_input_file)
    print*,'Molecule initial info pulled from: ',trim(adjustl(full_run_trajectory_input_file)) 
    call read_molformations_charges(moleculeformations_filename_full)
  else if (run_type == 3) then
    print*, "run_type = ", trim(adjustl(run_type_str)), " | mode set to: CUSTOM FRAGMENT AND VELOCITY INPUT"
    fragment_filename_full = trim(adjustl(fragment_input_path)) // 'fragment.inp'
    call read_fragment_input_file(fragment_filename_full)
    print*, "fragment file name = ", fragment_filename_full
    velocity_filename_full = trim(adjustl(velocity_input_path)) // 'velocity.inp'
    call read_velocity_input_file(velocity_filename_full)
    print*, "velocity file name = ", velocity_filename_full
  else
    print*, "ERROR: invalid run_type"
    print*, "  run_type = ", run_type, " Please set it to 1, 2, or 3"
    stop
  endif

  ! call compute atomic masses only if run_typ /= 3
  if (run_type /= 3) call compute_atomic_masses
  call open_optional_output_files

  call program_checks
  
  write(log_file,*)""
  write(log_file,*)"Initialization complete. Beginning computation..."
  write(log_file,*)""

  write(*,*)"Initialization complete. Beginning computation..."

END SUBROUTINE initialize


SUBROUTINE program_checks
  !TO-DO FINISH THIS
  integer :: number_of_lines

  if (output_atom_info .and. parallelization) then
    print*, 'ERROR, output_atom_info AND parallelization both set to true. Stopping.'
    stop
  end if

  if (run_type==1) then
    call count_lines(seeds_filename_full, number_of_lines)
    if (number_of_lines < N_simulations) then
      print*, 'ERROR, N_simulations is greater than the number of seeds (lines) in seeds.inp. Stopping.'
      stop
    endif
  endif

  if (run_type==2) then
    print*,"run_type=2 checks"
  endif

  
END SUBROUTINE program_checks


SUBROUTINE count_lines(filename, num_lines)
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(out) :: num_lines
  integer :: ios, line_count
  character(len=1000) :: line

  line_count = 0

  open(unit=10, file=filename, status='old', action='read', iostat=ios)
  if (ios /= 0) then
      print *, "Error opening file: ", filename
      num_lines = -1
      return
  end if

  do
      read(10, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line_count = line_count + 1
  end do

  close(10)
  num_lines = line_count
end SUBROUTINE count_lines



SUBROUTINE prepare_output
  call get_formatted_datetime(formatted_datetime)
  
  call check_and_create_directory(output_dir, dir_exists)
  if (.not. dir_exists) then
    write(*,'(A, A, A, A)') " Output directory created:  ", "'", trim(output_dir), "'"
  else
      write(*,'(A, A, A, A)') "'", trim(output_dir), "'", " directory found"
  end if

  output_dir_with_date_time = trim(adjustl(output_dir))//"/"//formatted_datetime

  call check_and_create_directory(output_dir_with_date_time, dir_exists)
  write(*,'(A)') " Output folder format:              DDMMYYYY_HHMMSS"
  write(*,'(A, A, A, A)') " Run output folder created: ", "'", trim(output_dir_with_date_time), "'"

  trajectory_directory = trim(adjustl(output_dir_with_date_time))//"/"//"trajectory"
  
  ! Adds the paths to each of the file_names
  call append_paths_to_filenames(formatted_datetime, trajectory_directory)
  
  open(log_file, file=trim(adjustl(log_output_filename)))
  open(all_variable_file, file=trim(adjustl(all_variable_filename)))

  write(log_file,*) "Run Started and output directory successfully created. Initializing..."
  write(log_file,'(A, A)') " Output folder created: ",  trim(output_dir_with_date_time)
  write(log_file,'(A, A)') " Output files created: ",  trim(output_dir_with_date_time)
  write(log_file,'(A, A, A, A)') " all variable file created: ", "'", trim(all_variable_filename), "'"
 

END SUBROUTINE prepare_output

SUBROUTINE open_optional_output_files
  if (output_trajectory) then
    call check_and_create_directory(trajectory_directory, dir_exists)
    write(*,'(A, A, A, A)') " Trajectory output folder created: ", "'", trim(trajectory_directory), "'"
  end if

  ! Open output files
  if (output_atom_info) open(atom_info_file, file=trim(adjustl(atom_info_filename)))
  ! open the trajectory file for each run-->open in the "calculate" subroutine

  if (output_atom_info) write(log_File,'(A, A, A, A)') " Atom info file created: ", "'", trim(atom_info_filename), "'"
  if (output_trajectory) write(log_file,'(A, A, A, A)') " Trajectory output folder created: ", "'", trim(trajectory_directory), "'"
END SUBROUTINE


SUBROUTINE calculate_simulation(i)
  integer, intent(in) :: i
  integer :: seed, unit_num
  character(len=255) :: full_trajectory_filename, seed_string

  ! Reset the seed and perform the computation again
  seed = seed_array(i)
  ion_velocity_init_seed = seed

  write(seed_string, '(G0)') seed
  if (output_trajectory) then
    full_trajectory_filename = trim(adjustl(trajectory_filename)) // &
                               "_r" // trim(adjustl(seed_string)) // &
                               ".xyz"
    unit_num = trajectory_output_file + 1 
    open(unit_num, file=trim(adjustl(full_trajectory_filename)))
  end if

  write(log_file, '(A, A, A)') "r=", trim(adjustl(seed_string)), &
                                " computation started"

  ! Initial state of each simulation
  atom_position = atom_initial_position
  call calculate_atomic_velocities
  call calculate_force

  do iter = 0, N_time_steps
      time = iter * time_step
  
      ! Output
      if (output_trajectory .and. mod(iter, trajectory_output_frequency) == 0) then
        call update_trajectory_file(unit_num)
      end if

      ! Update via Verlet Algorithm (force, position, velocity)
      call calculate_force
      call calculate_position
      call calculate_velocity
  end do

  if (output_trajectory) then
    write(log_file, *) "Successfully finished trajectory file: ", &
                      full_trajectory_filename
    close(unit_num)
  end if

  if (output_atom_info) call update_atom_info_file

  write(log_file, '(A, A, A)') "r=", trim(adjustl(seed_string)), &
                                " computation completed"
  write(log_file, *)
  write(*, '(A, A, A)') "  r=", trim(adjustl(seed_string)), &
                          " computation completed"
END SUBROUTINE calculate_simulation


SUBROUTINE simulate_cont_from_tddft(input_filename,i)
  character(len=256), intent(in) :: input_filename
  integer, intent(in) :: i
  integer :: j, unit_num
  character(len=255) :: full_trajectory_filename, seed_string

  !sam fix this later
  write(seed_string, '(G0)') i
  if (output_trajectory) then
    full_trajectory_filename = trim(adjustl(trajectory_filename)) // &
                               "_" // trim(adjustl(full_runs_array(i))) // &
                               ".xyz"
    unit_num = trajectory_output_file+(2*i)
    open(unit_num, file=trim(adjustl(full_trajectory_filename)))
  end if

  write(log_file, '(A, A, A)') trim(adjustl(full_runs_array(i))), &
                                " computation started"

  ! Initial state of each simulation
  call find_atom_kinetics_at_t(input_filename,i)

  ! Charge of each atom
  atom_charge = atom_charges_every_tddft(:,i)
  do j=1, N_total_atom
    write(log_file,*)"atom", j, "charge=",atom_charge(j)
  enddo

  do iter = traj_time_step_to_initialize, N_time_steps
      time = iter * time_step
  
      ! Output
      if (output_trajectory .and. mod(iter, trajectory_output_frequency) == 0) then
        call update_trajectory_file(unit_num)
      end if

      ! Update via Verlet Algorithm (force, position, velocity)
      call calculate_force
      call calculate_position
      call calculate_velocity
  end do

  if (output_trajectory) then
    write(log_file, *) "Successfully finished trajectory file: ", &
                      full_trajectory_filename
    close(unit_num)
  end if

  if (output_atom_info) call update_atom_info_file(i)

  write(log_file, '(A, A, A)') trim(adjustl(full_runs_array(i))), &
                                " computation completed"
  write(log_file, *)
  write(*, '(A, A, A)') "  ", trim(adjustl(full_runs_array(i))), &
                          " computation completed"
END SUBROUTINE simulate_cont_from_tddft


SUBROUTINE run_single_simulation
  integer :: i, unit_num
  character(len=255) :: full_trajectory_filename
  real*8 :: atom_velocity_magnitude

  if (output_trajectory) then
    full_trajectory_filename = trim(adjustl(trajectory_filename)) // ".xyz"
    unit_num = trajectory_output_file + 1 
    open(unit_num, file=trim(adjustl(full_trajectory_filename)))
  end if

  print*, "atom velocity: ", atom_velocity
  write(log_file, '(A)') "computation started"
  ! Initial state of each simulation
  atom_position = atom_initial_position
  call calculate_force

  do iter = time_step_start, N_time_steps + time_step_start
    time = iter * time_step

    ! Output
    if (output_trajectory .and. mod(iter, trajectory_output_frequency) == 0) then
      call update_trajectory_file(unit_num)
    end if

    ! Update via Verlet Algorithm (force, position, velocity)
    call calculate_force
    call calculate_position
    call calculate_velocity

  
    if (mod(iter,500) == 0) then
      atom_velocity_magnitude = sqrt(atom_velocity(1,1)**2 + &
                                      atom_velocity(2,1)**2 + &
                                      atom_velocity(3,1)**2)
      write(*,*) "iter=",iter," | atom_speed=",atom_velocity_magnitude
      write(*,*) "iter=",iter," | atom_force=",atom_force(1,1), atom_force(2,1), atom_force(3,1)
    end if

  end do

  if (output_trajectory) then
    write(log_file, *) "Successfully finished trajectory file: ", &
                      full_trajectory_filename
    close(unit_num)
  end if

  if (output_atom_info) call update_atom_info_file
END SUBROUTINE run_single_simulation


SUBROUTINE calculate
  use omp_lib
  integer :: i
  character(len=256) :: full_run_trajectory_input_file

  if (parallelization) then
    !$omp parallel private(i)
    !$omp master
    ! Print the number of threads used
    print *, "Running with parallelization, number of threads: ", omp_get_num_threads()
    !$omp end master

    !$omp do
    do i = 1, N_simulations
      call calculate_simulation(i)
    end do
    !$omp end do
    !$omp end parallel
  else
    ! In non-parallel mode, manually print the number of threads
    ! Here, typically, the number of threads is 1
    print *, "Running with no parallelization, number of threads: 1"
    
    if (run_type==1) then  
      do i = 1, N_simulations
        call calculate_simulation(i)
      end do
    else if (run_type==2) then
      do i = 1, N_simulations
        full_run_trajectory_input_file = trim(adjustl(full_runs_dir_input_path)) // &
                                     trim(adjustl(full_runs_array(i))) // "/trajectory.xyz"
        call simulate_cont_from_tddft(full_run_trajectory_input_file,i)
      end do
    else if (run_type==3) then
      call run_single_simulation
    endif
    
  endif

END SUBROUTINE calculate


SUBROUTINE cleanup
  real*8 :: program_CPU_total_time
  integer :: program_start_time
  integer :: program_end_time
  integer :: program_elapsed_time


  call date_and_time(date=end_date, time=end_time)
  call cpu_time(program_CPU_end_time)
  call time_to_seconds(start_date, start_time, program_start_time)
  call time_to_seconds(end_date, end_time, program_end_time)

  program_elapsed_time = program_end_time - program_start_time
  program_CPU_total_time = program_CPU_end_time - program_CPU_start_time

  write(*,*)"Run finished."
  if (parallelization) then
    write(*,*)" CPU time (seconds): ", program_CPU_total_time
    write(*,*)" Elapsed time (seconds): ", program_elapsed_time
  else
    write(*,*)" Elapsed time (seconds): ", program_CPU_total_time
  endif

  ! close the output files
  close(log_file)
  if (output_trajectory) close(trajectory_output_file)
  if (output_atom_info) close(atom_info_file)
  ! deallocate all of the matrices
  deallocate(atom_initial_position)
  deallocate(atom_position)
  deallocate(atom_velocity)
  deallocate(atom_acceleration)
  deallocate(atom_force)
  if (run_type==2) deallocate(atom_charges_every_tddft)
  ! deallocate all of the arrays
  deallocate(atom_atomic_number)
  deallocate(atom_charge)
  deallocate(atom_mass)
  if (run_type==1) deallocate(seed_array)
  if (run_type==2) deallocate(full_runs_array)


END SUBROUTINE cleanup


END MODULE COULOMB_EXPLOSION






