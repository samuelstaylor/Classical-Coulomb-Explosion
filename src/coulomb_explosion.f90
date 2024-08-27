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
  character(len=256) :: molecule_filename_full
  character(len=256) :: seeds_filename_full

  call date_and_time(date=start_date, time=start_time)
  call cpu_time(program_CPU_start_time)

  write(*,*) "-=================================-"
  write(*,*) "-=# CLASSICAL COULOMB EXPLOSION #=-"
  write(*,*) "-=================================-"
  write(*,*) "Run Started. Initializing..."

  call prepare_output
  call read_control_input_file('control.inp')

  molecule_filename_full = trim(adjustl(molecule_input_path)) // 'molecule.inp'
  call read_molecule_input_file(molecule_filename_full)

  seeds_filename_full = trim(adjustl(seeds_input_path)) // 'seeds.inp'
  call read_seeds_input_file(seeds_filename_full)

  print*,"molecule file name=",molecule_filename_full
  print*,"seeds file name=",seeds_filename_full

  call compute_atomic_masses
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

  call count_lines('seeds.inp', number_of_lines)
  if (number_of_lines < N_simulations) then
    print*, 'ERROR, N_simulations is greater than the number of seeds (lines) in seeds.inp. Stopping.'
    stop
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
                                trim(adjustl(seed_string)) // ".xyz"
    unit_num = trajectory_file + 1 
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


SUBROUTINE calculate
  use omp_lib
  integer :: i

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
    
    do i = 1, N_simulations
      call calculate_simulation(i)
    end do
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
  if (output_trajectory) close(trajectory_file)
  if (output_atom_info) close(atom_info_file)
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






