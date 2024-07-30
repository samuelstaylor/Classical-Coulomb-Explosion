

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

  call open_files
  call read_molecule_input_file('molecule.inp')
  call read_control_input_file('control.inp')


  ! global variables to not be set in control.inp
  ion_velocity_init_seed = 1
  ! OPTIONAL TO-DO: 
  ! maybe use an array of the initial seeds to generate from? or iterate it 0,1,2,..

  ! TO DO: Print to log file any info you want
  write(log_file,*)'All variables allocated and set, see: ', bare_all_variable_filename

  call compute_atomic_masses

END SUBROUTINE initialize


SUBROUTINE open_files
  character(len=15) :: formatted_datetime
  logical :: dir_exists
  character(len=255) :: output_dir_with_date_time


  call get_formatted_datetime(formatted_datetime)
  call append_paths_to_filenames(formatted_datetime)
  
  call check_and_create_directory(output_dir, dir_exists)
  if (.not. dir_exists) then
    write(*,'(A, A, A, A)') "'", trim(output_dir), "'", " directory does not exist. Creating directory..."
  else
      write(*,'(A, A, A, A)') "'", trim(output_dir), "'", " directory found."
  end if

  output_dir_with_date_time = trim(adjustl(output_dir))//"/"//formatted_datetime

  call check_and_create_directory(output_dir_with_date_time, dir_exists)
  write(*,'(A, A, A, A)') "output folder created: ", "'", trim(output_dir_with_date_time), "'"

  ! Open output files
  open(log_file,file=trim(adjustl(log_output_filename)))
  open(all_variable_file,file=trim(adjustl(all_variable_filename)))
  open(trajectory_file,file=trim(adjustl(trajectory_filename)))
  open(atom_info_file,file=trim(adjustl(atom_info_filename)))  !Boltzmann dist. to calculate the velocities

  write(log_file,*) "Run Started and output directory successfully created. Initializing..."


END SUBROUTINE open_files



SUBROUTINE get_formatted_datetime(output_string)
  implicit none
  character(len=15), intent(out) :: output_string
  character(len=8)  :: date
  character(len=8)  :: time
  character(len=5)  :: zone
  integer, dimension(8) :: values
  character(len=2) :: day_str, month_str, hour_str, minute_str, second_str
  character(len=4) :: year_str

  ! Obtain current date and time using keyword arguments
  call date_and_time(DATE=date, TIME=time, ZONE=zone, VALUES=values)

  ! Extract date and time components
  year_str = date(1:4)
  month_str = date(5:6)
  day_str = date(7:8)
  hour_str = time(1:2)
  minute_str = time(3:4)
  second_str = time(5:6)

  ! Create formatted date and time string
  output_string = trim(day_str) // trim(month_str) // trim(year_str) // '_' // &
                  trim(hour_str) // trim(minute_str) // trim(second_str)
END SUBROUTINE get_formatted_datetime


! Calculates the Coulomb Force at each time step in the simulation
! and propagates it throughout time 
SUBROUTINE calculate    
  ! Initialize the atom's velocity via Boltzmann distribution
  call calculate_atomic_velocities
  call calculate_force

  do iter=0, N_time_steps
      time = iter * time_step
  
      ! Output
      if (mod(iter,trajectory_output_frequency) == 0) then
        call update_trajectory_file(iter,time)
        call print_to_log_file
      end if

      ! Update via Verlet Algorithm (force, position, velocity)
      call calculate_force
      call calculate_position
      call calculate_velocity
  end do

  call update_atom_info_file

END SUBROUTINE calculate


subroutine cleanup
  real*8 :: program_elapsed_time

  ! close the output files
  close(log_file)
  close(trajectory_file)
  close(atom_info_file)
  ! deallocate all of the matrices
  deallocate(atom_position)
  deallocate(atom_velocity)
  deallocate(atom_acceleration)
  deallocate(atom_force)
  ! deallocate all of the arrays
  deallocate(atom_atomic_number)
  deallocate(atom_charge)
  deallocate(atom_mass)

  call cpu_time(program_end_time)
  program_elapsed_time = program_end_time - program_start_time
  write(*,*)"Run finished. Elapsed time (s): ", program_elapsed_time

end subroutine cleanup

END MODULE COULOMB_EXPLOSION






