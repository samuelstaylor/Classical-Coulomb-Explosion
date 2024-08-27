! Module that contains all of the input output operations used in the program

MODULE MOD_IO
  USE MOD_GLOBAL
  implicit none

CONTAINS


! INPUT SUBROUTINES
SUBROUTINE read_molecule_input_file(input_filename)
  character(len=*), intent(in) :: input_filename
  integer :: i, error_code, atom_counter
  character(len=256) :: line

  ! Open the file
  open(unit=molecule_file, file=trim(adjustl(input_filename)), status='old', action='read', iostat=error_code)
  if (error_code /= 0) then
    write(*,*) "***ERROR*** Error opening: ", input_filename
    write(*,*) "Ensure that ", input_filename, " is in your running directory"
    write(log_file,*) "***ERROR*** Error opening: ", input_filename
    write(log_file,*) "Ensure that ", input_filename, " is in your running directory"
    stop
  end if

  ! Read lines and skip those starting with '#'
  do
      read(molecule_file, '(A)', iostat=error_code) line
      if (error_code /= 0) exit

      ! Ignore anything on the line that appears after a '#'
      i = index(line, '#')
      if (i > 0) then
          line = line(1:i-1)
      end if

      ! Trim leading and trailing spaces
      line = adjustl(line)
      if (len_trim(line) == 0) cycle  ! Skip empty lines

      ! First non-comment line should be the number of atoms
      read(line, *, iostat=error_code) N_total_atom
      if (error_code == 0) exit
      if (error_code /= 0) then
      write(*,*) "***ERROR*** Error reading N_total_atom from ", input_filename
      write(log_file,*) "***ERROR*** Error reading N_total_atom from ", input_filename
      stop
      end if
  end do

  ! 1-dimensional arrays that will be allocated to the size N_total_atom
  allocate(atom_atomic_number(N_total_atom))
  allocate(atom_charge(N_total_atom))
  allocate(atom_mass(N_total_atom))

  ! 2-dimensional arrays that will be allocated to N_total_atoms columns
  allocate(atom_force(3,N_total_atom))
  allocate(atom_acceleration(3,N_total_atom))
  allocate(atom_velocity(3,N_total_atom))
  allocate(atom_position(3,N_total_atom))
  allocate(atom_initial_position(3,N_total_atom))
 
  ! Initialize all arrays to zero
  atom_initial_position = 0.0
  atom_position = 0.0
  atom_velocity = 0.0
  atom_acceleration = 0.0
  atom_force = 0.0
  atom_atomic_number = 0
  atom_mass = 0.0

  atom_counter = 0
  ! Read lines and skip those starting with '#'
  do while (atom_counter < N_total_atom)
    read(molecule_file, '(A)', end=110,iostat=error_code) line
    if (error_code /= 0) exit

    ! Ignore anything on the line that appears after a '#'
    i = index(line, '#')
    if (i > 0) then
        line = line(1:i-1)
    end if

    ! Trim leading and trailing spaces
    line = adjustl(line)
    if (len_trim(line) == 0) cycle  ! Skip empty lines

    ! First non-comment line should be the number of atoms
    ! Read the atom positions and atomic numbers
    atom_counter = atom_counter + 1        

    read(line, *, iostat=error_code) atom_initial_position(1, atom_counter), atom_initial_position(2, atom_counter), &
                atom_initial_position(3, atom_counter), atom_atomic_number(atom_counter), atom_charge(atom_counter)
    if (error_code /= 0) then
        write(*,*) "***ERROR*** Error reading atom data at atom: ", atom_counter
        write(log_file,*) "***ERROR*** Error reading atom data at atom: ", atom_counter
        stop
    end if
  end do

  ! Close the file
  110 close(molecule_file)  

  write(log_file,*) "Successfully read and initialized data from: ", input_filename
  call print_initial_atom_info_to_log

END SUBROUTINE read_molecule_input_file


SUBROUTINE read_control_input_file(input_filename)
  character(len=*), intent(in) :: input_filename
  character*256 :: file_line, the_key, value_string
  integer :: i, error_code

  ! Open the file
  open(unit=control_file, file=trim(adjustl(input_filename)), status='old', action='read', iostat=error_code)
  if (error_code /= 0) then
      write(log_file,*) "***ERROR*** Error opening: ", input_filename
      write(log_file,*) "***ERROR*** Ensure that ", input_filename, " is in your running directory"
      stop 
  end if

  do
    read(control_file,'(a)',end=110)file_line

    ! Ignore anything on the line that appears after a #
    i=index(file_line,'#')
    if(i>0) then
      file_line=file_line(1:i-1)
    endif

    i=index(file_line,'=')
    if(i==0) then
      cycle
    else
      read(file_line(1:i-1),'(a)')the_key
      read(file_line(i+1:),'(a)')value_string
      call convert_to_lowercase(the_key)
    
      select case(trim(adjustl(the_key)))
        case("molecule_input_path")
          molecule_input_path = trim(adjustl(value_string))
        case("seeds_input_path")
          seeds_input_path = trim(adjustl(value_string))       
        case("n_simulations")
          read(value_string,*)N_simulations
        case("n_time_steps")
          read(value_string,*)N_time_steps
        case("time_step")
          read(value_string,*)time_step 
        case("trajectory_output_frequency")
          read(value_string,*)trajectory_output_frequency
        case("temperature_ions")
          read(value_string,*)temperature_ions
        case("use_average_atomic_mass")
          read(value_string,*)use_average_atomic_mass 
        case("include_electron_mass")
          read(value_string,*)include_electron_mass
        case("output_trajectory")
          read(value_string,*)output_trajectory
        case("output_atom_info")
          read(value_string,*)output_atom_info 
        case("parallelization")
          read(value_string,*)parallelization 
        case default
          write(log_file,*)"ERROR: Invalid variable name: ",trim(adjustl(the_key))
          write(*,*)"ERROR: Invalid variable name: ",trim(adjustl(the_key))
          stop 		   
      end select
    endif
  enddo
  110 close(control_file)

  call check_and_fix_paths

  ! full list of variables
  write(all_variable_file,*) "  molecule_input_path=", molecule_input_path
  write(all_variable_file,*) "  seeds_input_path=", seeds_input_path
  write(all_variable_file,*) "  N_simulations=", N_simulations
  write(all_variable_file,*) "  N_time_steps=", N_time_steps
  write(all_variable_file,*) "  time_step=", time_step
  write(all_variable_file,*) "  trajectory_output_frequency=", trajectory_output_frequency
  write(all_variable_file,*) "  temperature_ions=", temperature_ions
  write(all_variable_file,*) "  use_average_atomic_mass=", use_average_atomic_mass
  write(all_variable_file,*) "  include_electron_mass=", include_electron_mass
  write(all_variable_file,*) "  output_trajectory=", output_trajectory
  write(all_variable_file,*) "  output_atom_info=", output_atom_info

  close(all_variable_file)

  write(log_file,*) "Successfully read and initialized data from: ", input_filename
  write(log_file,*)'All variables allocated and set, see: ', bare_all_variable_filename
  
END SUBROUTINE read_control_input_file


! Read in each seed from the 'seeds.inp' input file
SUBROUTINE read_seeds_input_file(input_filename)
  character(len=*), intent(in) :: input_filename
  integer :: i, error_code, seed_counter
  character(len=256) :: line

  ! Open the file
  open(unit=seeds_file, file=trim(adjustl(input_filename)), status='old', action='read', iostat=error_code)
  if (error_code /= 0) then
    write(*,*) "***ERROR*** Error opening: ", input_filename
    write(*,*) "Ensure that ", input_filename, " is in your running directory"
    write(log_file,*) "***ERROR*** Error opening: ", input_filename
    write(log_file,*) "Ensure that ", input_filename, " is in your running directory"
    stop
  end if

  ! 1-dimensional arrays that will be allocated to the size N_total_atom
  allocate(seed_array(N_simulations))
  
  seed_counter = 0
  seed_array = 0
  ! Read lines and skip those starting with '#'
  do while (seed_counter < N_simulations)
    read(seeds_file, '(A)', end=110,iostat=error_code) line
    if (error_code /= 0) exit

    ! Ignore anything on the line that appears after a '#'
    i = index(line, '#')
    if (i > 0) then
        line = line(1:i-1)
    end if

    ! Trim leading and trailing spaces
    line = adjustl(line)
    if (len_trim(line) == 0) cycle  ! Skip empty lines

    ! First non-comment line should be the number of atoms
    ! Read the atom positions and atomic numbers
    seed_counter = seed_counter + 1       

    read(line, *, iostat=error_code) seed_array(seed_counter)
    if (error_code /= 0) then
        write(*,*) "***ERROR*** Error reading seed data at seed line: ", seed_counter
        write(log_file,*) "***ERROR*** Error reading seed data at seed line: ", seed_counter
        stop
    end if
  end do

  ! Close the file
  110 close(seeds_file)  

  ! Error check
  if (N_simulations > seed_counter) then
    write(*,*) "***ERROR*** Error N_simulations must be <= the number of seeds in ", input_filename
    write(log_file,*) "***ERROR*** Error N_simulations must be <= the number of seeds in ", input_filename
    stop
  end if
  write(log_file,*) "Successfully read and initialized data from: ", input_filename

END SUBROUTINE read_seeds_input_file


SUBROUTINE convert_to_lowercase(string)
  integer          :: i
  character(len=*) :: string
  character        :: the_char

  do i=1,len(string)
    the_char=string(i:i)
    if((iachar(the_char)>=iachar("A")).AND.(iachar(the_char)<=iachar("Z"))) then
      the_char=achar(iachar(the_char)+(iachar('a')-iachar('A')))
    endif
    string(i:i)=the_char
  enddo
END SUBROUTINE convert_to_lowercase


SUBROUTINE check_and_fix_paths
  ! Ensure dft_input_path and velocity_output_path start and end with '/'
  if (len_trim(molecule_input_path) > 0 .and. molecule_input_path(1:1) /= ".") then
    if (molecule_input_path(1:1) /= '/') then
      molecule_input_path = '/' // trim(adjustl(molecule_input_path))
    endif
    if (molecule_input_path(len_trim(molecule_input_path):len_trim(molecule_input_path)) /= '/') then
      molecule_input_path = trim(adjustl(molecule_input_path)) // '/'
    endif
  endif
  if (len_trim(seeds_input_path) > 0 .and. seeds_input_path(1:1) /= ".") then
    if (seeds_input_path(1:1) /= '/') then
      seeds_input_path = '/' // trim(adjustl(seeds_input_path))
    endif
    if (seeds_input_path(len_trim(seeds_input_path):len_trim(seeds_input_path)) /= '/') then
      seeds_input_path = trim(adjustl(seeds_input_path)) // '/'
    endif
  endif

END SUBROUTINE check_and_fix_paths


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

SUBROUTINE append_paths_to_filenames(formatted_datetime, trajectory_directory)
  character(len=15), intent(in) :: formatted_datetime
  character(len=33), intent(in) :: trajectory_directory
  
  log_output_filename=trim(adjustl(output_dir))//"/"//formatted_datetime//"/"//bare_log_output_filename
  all_variable_filename=trim(adjustl(output_dir))//"/"//formatted_datetime//"/"//bare_all_variable_filename
  atom_info_filename=trim(adjustl(output_dir))//"/"//formatted_datetime//"/"//bare_atom_info_filename
  trajectory_filename=trim(adjustl(trajectory_directory))//"/"//bare_trajectory_filename

END SUBROUTINE append_paths_to_filenames


SUBROUTINE print_masses_to_log
  integer :: i

  write(log_file,*) 'Atom atomic number and mass (eV_fs^2/A^2):'
  do i=1, N_total_atom
      write(log_File,*) atom_atomic_number(i), atom_mass(i)
  end do

END SUBROUTINE print_masses_to_log


SUBROUTINE print_initial_atom_info_to_log
  integer :: i, elem

  write(log_file,*) 'Element, Atom atomic number, and charge:'
  do i=1, N_total_atom
    elem=atom_atomic_number(i)
    write(log_File,*) element_symbols(elem), elem, atom_mass(i)
  end do

  write(log_file,*) 'Atom positions (A) in x,y,z:'
  do i=1, N_total_atom
    write(log_file,*) atom_position(:,i)
  end do

END SUBROUTINE print_initial_atom_info_to_log


! OUTPUT SUBROUTINES
SUBROUTINE print_to_terminal
  ! Printing all dynamic information to the terminal. Useful for debugging
  integer :: i

  write(*,*) '#iter=', iter, "#time=", time

  write(*,*) 'Atom positions (A) in x,y,z:'
  do i=1, N_total_atom
      write(*,*) atom_position(:,i)
  end do

  write(*,*) 'Atom velocities (A/fs) in x,y,z:'
  do i=1, N_total_atom
      write(*,*) atom_velocity(:,i)
  end do

  write(*,*) 'Atom acceleration (A/fs^2) in x,y,z:'
  do i=1, N_total_atom
      write(*,*) atom_acceleration(:,i)
  end do

  write(*,*) 'Atom force (eV/A) in x,y,z:'
  do i=1, N_total_atom
      write(*,*) atom_force(:,i)
  end do
END SUBROUTINE print_to_terminal


SUBROUTINE print_to_log_file
  integer :: i, elem

  write(log_file,*) '#iter=', iter, '#time=', time

  write(log_file,*) 'Atom positions (A) in x,y,z:'
  do i=1, N_total_atom
    elem=atom_atomic_number(i)
    write(log_file,*) element_symbols(elem), atom_position(:,i)
  end do

  write(log_file,*) 'Atom velocities (A/fs) in x,y,z:'
  do i=1, N_total_atom
    elem=atom_atomic_number(i)
    write(log_file,*) element_symbols(elem), atom_velocity(:,i)
  end do

  write(log_file,*) 'Atom speed magnitude (A/fs):'
  do i=1, N_total_atom
    elem=atom_atomic_number(i)
    write(log_file,*) element_symbols(elem), SQRT(SUM(atom_velocity(:,i)**2))
  end do

  write(log_file,*) 'Atom acceleration (A/fs^2) in x,y,z:'
  do i=1, N_total_atom
    elem=atom_atomic_number(i)
    write(log_file,*) element_symbols(elem), atom_acceleration(:,i)
  end do

  write(log_file,*) 'Atom force (eV/A) in x,y,z:'
  do i=1, N_total_atom
    elem=atom_atomic_number(i)
    write(log_file,*) element_symbols(elem), atom_force(:,i)
  end do

  write(log_file,*) 'Atom force magnitude (A/fs^2):'
  do i=1, N_total_atom
    elem=atom_atomic_number(i)
    write(log_file,*) element_symbols(elem), SQRT(SUM(atom_force(:,i)**2))
  end do

  write(log_file,*) ""
END SUBROUTINE print_to_log_file

!for generating velocity.inp files
SUBROUTINE write_velocity_to_file(output_filename)
  character(len=*), intent(in) :: output_filename
  integer :: i, unit_num, error_code

  ! Open the file
  unit_num = 20
  open(unit=unit_num, file=trim(adjustl(output_filename)), status='replace', action='write', iostat=error_code)
  if (error_code /= 0) then
      write(*,*) "Error opening output file writing: ", output_filename
      stop
  end if

  ! write velocities to the file
  do i = 1, N_total_atom
      write(unit_num, *) atom_velocity(1, i), atom_velocity(2, i), atom_velocity(3, i)
  end do

  ! Close the file
  close(unit_num)

  write(log_file,*) "Successfully wrote velocities to file: ", output_filename
END SUBROUTINE write_velocity_to_file


subroutine update_trajectory_file(unit_num)
  !This subroutine updates xyz-file with atomic trajectories
  !The trajectory xyz-file is basically a big xyz-file containing
  !all consequtive snapshots (no blank lines between snapshots).
  !In the second line of each snapshot (second line of an xyz-file
  !can be anything and is ignored by most visualization software) we
  !put a comment line with the time [in femtoseconds] when the snapshot
  !was made
  !Input parameters:
  !  unit_num -- the unit number of the file to open
  integer, intent(in) :: unit_num 
  integer :: i,elem
  
  write(unit_num,'(1x,i5)') N_total_atom
  write(unit_num,'(1x,a,i6,a,f8.5)') '# iter =',iter,'  time[fs]=',time
  do i=1,N_total_atom
    elem=atom_atomic_number(i)
    write(unit_num,'(1x,a3,3(1x,e16.9))') &
        element_symbols(elem),atom_position(1,i),atom_position(2,i),atom_position(3,i)
  enddo
  flush(unit_num)
  end subroutine update_trajectory_file
  
SUBROUTINE update_atom_info_file
  !This subroutine updates csv-file with info about each of the atoms
  !from each of the simulations. 
  integer :: i, error_code, atom_zero_index
  character(len=255) :: atom_zero_index_string
  real*8 :: charge_sum

  ! write the atom info header
  WRITE(atom_info_file, '(A, I0, A)', ADVANCE='NO') "r_seed=", ion_velocity_init_seed, ", "
  do i = 1, N_total_atom
    atom_zero_index = i - 1
    write(atom_zero_index_string,*) atom_zero_index
    write(atom_info_file, '(A, A, A, A)', ADVANCE='NO') element_symbols(atom_atomic_number(i)), "[", &
                                                        trim(adjustl(atom_zero_index_string)), "], "
  end do
  write(atom_info_file,*) ""

  ! write charges
  write(atom_info_file, '(A)', ADVANCE='NO') "Charges, "
  do i = 1, N_total_atom
      write(atom_info_file, '(F15.10, A)', ADVANCE='NO') atom_charge(i), ", "
  end do
  write(atom_info_file,*) ""

  ! write Time
  write(atom_info_file, '(A)', ADVANCE='NO') "Time[fs], "
  do i = 1, N_total_atom
      write(atom_info_file, '(F10.6, A)', ADVANCE='NO') time, ", "
  end do
  write(atom_info_file,*) ""

  ! write charge Sum
  charge_sum = SUM(atom_charge)
  write(atom_info_file, '(A, F15.10, A)', ADVANCE='NO') "Charge Sum, ", charge_sum, ", "
  write(atom_info_file,*)

  ! write X Velocity
  write(atom_info_file, '(A)', ADVANCE='NO') "X Velocity[A/fs], "
  do i = 1, N_total_atom
      write(atom_info_file, '(F10.6, A)', ADVANCE='NO') atom_velocity(1, i), ", "
  END do
  write(atom_info_file,*) ""

  ! write Y Velocity
  write(atom_info_file, '(A)', ADVANCE='NO') "Y Velocity[A/fs], "
  do i = 1, N_total_atom
      write(atom_info_file, '(F10.6, A)', ADVANCE='NO') atom_velocity(2, i), ", "
  END do
  write(atom_info_file, '(A)') ""

  ! write Z Velocity
  write(atom_info_file, '(A)', ADVANCE='NO') "Z Velocity[A/fs], "
  do i = 1, N_total_atom
      write(atom_info_file, '(F10.6, A)', ADVANCE='NO') atom_velocity(3, i), ", "
  END do
  write(atom_info_file,*) ""

  ! write Speed
  write(atom_info_file, '(A)', ADVANCE='NO') "Speed[A/fs], "
  do i = 1, N_total_atom
      write(atom_info_file, '(F10.6, A)', ADVANCE='NO') SQRT(SUM(atom_velocity(:, i)**2)), ", "
  end do
  write(atom_info_file,*) ""
  write(atom_info_file,*) ""

  write(log_file,*) "Successfully updated atom info file: ", atom_info_filename


END SUBROUTINE update_atom_info_file

END MODULE MOD_IO
