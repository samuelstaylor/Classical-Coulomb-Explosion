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

  call allocate_and_initialize_arrays

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


SUBROUTINE allocate_and_initialize_arrays
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

  if (run_type==2) then
    allocate(atom_charges_every_tddft(N_total_atom,N_simulations))
    atom_charges_every_tddft=0.0
  endif
  
  if (finalize_info_at_cap) then
    allocate(times_at_cap(N_total_atom))
    allocate(atom_force_final(3,N_total_atom))
    allocate(atom_acceleration_final(3,N_total_atom))
    allocate(atom_velocity_final(3,N_total_atom))
    allocate(atom_position_final(3,N_total_atom))

    times_at_cap=0.0
    atom_force_final=0.0
    atom_acceleration_final=0.0
    atom_velocity_final=0.0
    atom_position_final=0.0
  endif

  ! Initialize all arrays to zero
  atom_initial_position = 0.0
  atom_position = 0.0
  atom_velocity = 0.0
  atom_acceleration = 0.0
  atom_force = 0.0
  atom_atomic_number = 0
  atom_mass = 0.0
END SUBROUTINE allocate_and_initialize_arrays


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
        case("run_type")
          read(value_string,*)run_type
        case("molecule_input_path")
          molecule_input_path = trim(adjustl(value_string))
        case("seeds_input_path")
          seeds_input_path = trim(adjustl(value_string))
        case("moleculeformations_input_path")
          moleculeformations_input_path = trim(adjustl(value_string)) 
        case("full_runs_dir_input_path")
          full_runs_dir_input_path = trim(adjustl(value_string))   
        case("traj_time_step_to_initialize")
          read(value_string,*)traj_time_step_to_initialize      
        case("fragment_input_path")
          fragment_input_path = trim(adjustl(value_string))   
        case("velocity_input_path")
          velocity_input_path = trim(adjustl(value_string))   
        case("time_step_start")
          read(value_string,*)time_step_start
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
        case("finalize_info_at_cap")
          read(value_string,*)finalize_info_at_cap
        case("cap_tolerance")
          read(value_string,*)cap_tolerance
        case("cap_left1")
          read(value_string,*)cap_left1
        case("cap_right1")
          read(value_string,*)cap_right1
        case("cap_left2")
          read(value_string,*)cap_left2
        case("cap_right2")
          read(value_string,*)cap_right2
        case("cap_left3")
          read(value_string,*)cap_left3
        case("cap_right3")
          read(value_string,*)cap_right3
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
  write(all_variable_file,*) "  run_type=", run_type
  write(all_variable_file,*) "  molecule_input_path=", molecule_input_path
  write(all_variable_file,*) "  seeds_input_path=", seeds_input_path
  write(all_variable_file,*) "  moleculeformations_input_path=", moleculeformations_input_path
  write(all_variable_file,*) "  full_runs_dir_input_path=", full_runs_dir_input_path
  write(all_variable_file,*) "  traj_time_step_to_initialize=", traj_time_step_to_initialize
  write(all_variable_file,*) "  fragment_input_path=", fragment_input_path
  write(all_variable_file,*) "  velocity_input_path=", velocity_input_path
  write(all_variable_file,*) "  time_step_start=", time_step_start
  write(all_variable_file,*) "  N_simulations=", N_simulations
  write(all_variable_file,*) "  N_time_steps=", N_time_steps
  write(all_variable_file,*) "  time_step=", time_step
  write(all_variable_file,*) "  trajectory_output_frequency=", trajectory_output_frequency
  write(all_variable_file,*) "  temperature_ions=", temperature_ions
  write(all_variable_file,*) "  use_average_atomic_mass=", use_average_atomic_mass
  write(all_variable_file,*) "  include_electron_mass=", include_electron_mass
  write(all_variable_file,*) "  output_trajectory=", output_trajectory
  write(all_variable_file,*) "  output_atom_info=", output_atom_info
  write(all_variable_file,*) "  finalize_info_at_cap=", finalize_info_at_cap
  write(all_variable_file,*) "  cap_tolerance=", cap_tolerance
  write(all_variable_file,*) "  cap_left1=", cap_left1
  write(all_variable_file,*) "  cap_right1=", cap_right1
  write(all_variable_file,*) "  cap_left2=", cap_left2
  write(all_variable_file,*) "  cap_right2=", cap_right2
  write(all_variable_file,*) "  cap_left3=", cap_left3
  write(all_variable_file,*) "  cap_right3=", cap_right3

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


SUBROUTINE read_molformations_input_file(input_filename)
  character(len=*), intent(in) :: input_filename
  integer :: i, error_code, run_counter
  character(len=512) :: line,lower_line
  logical :: match_found
  character(len=4) :: N_simulations_str

  ! Open the file to find number of simulations in the file
  open(unit=moleculeformations_file, file=trim(adjustl(input_filename)), status='old', action='read', iostat=error_code)
  if (error_code /= 0) then
    write(*,*) "***ERROR*** Error opening: ", input_filename
    write(*,*) "Ensure that ", input_filename, " is in your running directory"
    write(log_file,*) "***ERROR*** Error opening: ", input_filename
    write(log_file,*) "Ensure that ", input_filename, " is in your running directory"
    stop
  end if

  N_simulations=0
  do
    read(moleculeformations_file, '(a)', end=110) line

    ! Ignore anything on the line that appears after a ,
    i = index(line, ',')
    if (i > 0) then
        line = trim(adjustl(line(1:i-1)))
    endif

    ! Skip to the next iteration if the line is empty
    if (len(trim(line)) == 0) cycle

    call convert_to_lowercase(line)

    if (trim(adjustl(line)) == "speed[a/fs]") then
        N_simulations = N_simulations + 1
    endif
  enddo
  110 close(moleculeformations_file)

  write(N_simulations_str, '(I0)') N_simulations   ! Write integer N_simulations to string without leading spaces
  print*, 'Number of simulations is: ', trim(adjustl(N_simulations_str))


  ! 1-dimensional arrays that will be allocated to the size N_total_atom
  allocate(full_runs_array(N_simulations))

  ! Open the file to find fill the array of full run names
  open(unit=moleculeformations_file, file=trim(adjustl(input_filename)), status='old', action='read', iostat=error_code)
  if (error_code /= 0) then
    write(*,*) "***ERROR*** Error opening: ", input_filename
    write(*,*) "Ensure that ", input_filename, " is in your running directory"
    write(log_file,*) "***ERROR*** Error opening: ", input_filename
    write(log_file,*) "Ensure that ", input_filename, " is in your running directory"
    stop
  end if

  full_runs_array=""
  run_counter=1
  do
    read(moleculeformations_file, '(a)', end=111) line

    ! Ignore anything on the line after a comma
    i = index(line, ',')
    if (i > 0) then
        line = line(1:i-1)
    endif

    ! Trim and adjust the line once, reuse it
    line = trim(adjustl(line))

    lower_line=line
    call convert_to_lowercase(lower_line)

    ! Check if the line matches any of the excluded conditions
    match_found = (trim(adjustl(lower_line)) == "") .or. &
                  (trim(adjustl(lower_line)) == "densities") .or. &
                  (trim(adjustl(lower_line)) == "time[fs]") .or. &
                  (trim(adjustl(lower_line)) == "density sum") .or. &
                  (trim(adjustl(lower_line)) == "x velocity[a/fs]") .or. &
                  (trim(adjustl(lower_line)) == "y velocity[a/fs]") .or. &
                  (trim(adjustl(lower_line)) == "z velocity[a/fs]") .or. &
                  (trim(adjustl(lower_line)) == "speed[a/fs]")

    if (.not. match_found) then
        ! Check if we need to reallocate array size
        if (run_counter > size(full_runs_array)) then
          print*, "***ERROR*** More runs found than in size of full_runs_array"
          stop
        endif
        full_runs_array(run_counter) = line
        run_counter = run_counter + 1
    endif
  enddo
  111 close(moleculeformations_file)

  ! Error check
  if (N_simulations > run_counter) then
    write(*,*) "***ERROR*** Less runs found than are in full_runs_array from: ", input_filename
    write(log_file,*) "***ERROR*** Less runs found than are in full_runs_array from: ", input_filename
    stop
  end if
  write(log_file,*) "Successfully read and initialized data from: ", input_filename

END SUBROUTINE read_molformations_input_file


SUBROUTINE read_molformations_charges(input_filename)
  character(len=*), intent(in) :: input_filename
  integer :: i, j, error_code, run_counter
  character(len=512) :: line,lower_line,line_header,line_body
  logical :: match_found
  real :: line_array(N_total_atom)
  real :: temp_density

  
  ! Open the file to find fill the array of full run names
  open(unit=moleculeformations_file, file=trim(adjustl(input_filename)), status='old', action='read', iostat=error_code)
  if (error_code /= 0) then
    write(*,*) "***ERROR*** Error opening: ", input_filename
    write(*,*) "Ensure that ", input_filename, " is in your running directory"
    write(log_file,*) "***ERROR*** Error opening: ", input_filename
    write(log_file,*) "Ensure that ", input_filename, " is in your running directory"
    stop
  end if

  run_counter=1
  do
    read(moleculeformations_file, '(a)', end=110) line

    ! Ignore anything on the line after a comma
    i = index(line, ',')
    if (i > 0) then
        line_header = line(1:i-1)
        line_body= line(i+1:)
    endif

    ! Trim and adjust the line once, reuse it
    line_header = trim(adjustl(line_header))

    lower_line=line_header
    call convert_to_lowercase(lower_line)

    if (trim(adjustl(lower_line)) == "densities") then
      ! Split line_body by commas
      call split_line_to_real_array(line_body,line_array)
    
      ! Check if we need to reallocate atom_charge size based on N_simulations
      do i = 1, N_total_atom
        ! Convert the i-th element from line_array to a real and assign it to atom_charge
        temp_density=line_array(i)
        ! fill with densities
        atom_charges_every_tddft(i,run_counter) = temp_density
      end do
      run_counter = run_counter + 1
    endif
  enddo
  110 close(moleculeformations_file)

  !convert the densities to charges
  call convert_density_to_charges
END SUBROUTINE read_molformations_charges


SUBROUTINE read_fragment_input_file(input_filename)
  character(len=*), intent(in) :: input_filename
  integer :: i, error_code, atom_counter
  character(len=256) :: line

  ! Open the file
  open(unit=fragment_file, file=trim(adjustl(input_filename)), status='old', action='read', iostat=error_code)
  if (error_code /= 0) then
    write(*,*) "***ERROR*** Error opening: ", input_filename
    write(*,*) "Ensure that ", input_filename, " is in your running directory"
    write(log_file,*) "***ERROR*** Error opening: ", input_filename
    write(log_file,*) "Ensure that ", input_filename, " is in your running directory"
    stop
  end if

  ! Read lines and skip those starting with '#'
  do
    read(fragment_file, '(A)', iostat=error_code) line
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

  call allocate_and_initialize_arrays

  atom_counter = 0
  ! Read lines and skip those starting with '#'
  do while (atom_counter < N_total_atom)
    read(fragment_file, '(A)', end=110,iostat=error_code) line
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
                atom_initial_position(3, atom_counter), atom_mass(atom_counter), atom_charge(atom_counter)
    if (error_code /= 0) then
        write(*,*) "***ERROR*** Error reading atom data at atom: ", atom_counter
        write(log_file,*) "***ERROR*** Error reading atom data at atom: ", atom_counter
        stop
    end if
  end do

  ! Close the file
  110 close(fragment_file) 

  !initalize all elements in atom_atomic_number to 1 "hydrogen
  atom_atomic_number=1

  write(log_file,*) "Successfully read and initialized data from: ", input_filename
  call print_initial_atom_info_to_log

END SUBROUTINE read_fragment_input_file


SUBROUTINE read_velocity_input_file(input_filename)
  character(len=*), intent(in) :: input_filename
  integer :: i, error_code, atom_counter, num_atom
  character(len=256) :: line

  ! Open the file
  open(unit=velocity_file, file=trim(adjustl(input_filename)), status='old', action='read', iostat=error_code)
  if (error_code /= 0) then
    write(*,*) "***ERROR*** Error opening: ", input_filename
    write(*,*) "Ensure that ", input_filename, " is in your running directory"
    write(log_file,*) "***ERROR*** Error opening: ", input_filename
    write(log_file,*) "Ensure that ", input_filename, " is in your running directory"
    stop
  end if

  ! Read lines and skip those starting with '#'
  do
      read(velocity_file, '(A)', iostat=error_code) line
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
      read(line, *, iostat=error_code) num_atom
      if (error_code == 0) exit
      if (error_code /= 0) then
      write(*,*) "***ERROR*** Error reading num_atom from ", input_filename
      write(log_file,*) "***ERROR*** Error reading num_atom from ", input_filename
      stop
      end if
  end do

  ! CHECK THAT N_total_atom = num_atom
  if (num_atom /= N_total_atom) then
    write(*,*) "***ERROR*** N_total_atom must equal num_atom from ", input_filename
    write(log_file,*) "***ERROR*** N_total_atom must equal num_atom from ", input_filename
    stop
  end if

  atom_counter = 0
  ! Read lines and skip those starting with '#'
  do while (atom_counter < N_total_atom)
    read(velocity_file, '(A)', end=110,iostat=error_code) line
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

    read(line, *, iostat=error_code) atom_velocity(1, atom_counter), atom_velocity(2, atom_counter), &
                atom_velocity(3, atom_counter)
    if (error_code /= 0) then
        write(*,*) "***ERROR*** Error reading atom data at atom: ", atom_counter
        write(log_file,*) "***ERROR*** Error reading atom data at atom: ", atom_counter
        stop
    end if
  end do

  ! Close the file
  110 close(velocity_file)  

  write(log_file,*) "Successfully read and initialized data from: ", input_filename

END SUBROUTINE read_velocity_input_file


SUBROUTINE convert_density_to_charges
  integer :: i, elem

  ! Calculates assuming a neutral charge and most common isotope
  do i=1,N_total_atom
    elem=atom_atomic_number(i)
    if (elem/=0) then
      atom_charges_every_tddft(i,:)=element_valence_electrons(elem)-atom_charges_every_tddft(i,:)
    endif
  enddo
END SUBROUTINE convert_density_to_charges


SUBROUTINE split_line_to_real_array(line_body, line_array_out)
  implicit none
  ! Inputs
  character(len=*), intent(in) :: line_body
  ! Outputs
  real, intent(out) :: line_array_out(N_total_atom)
  
  ! Local variables
  character(len=32) :: line_array(N_total_atom)  ! Assuming each value is at most 20 characters long
  integer :: i, j, start, num_values
  real :: temp_real
  
  ! Initialize
  num_values = 0
  start = 1

  ! Loop over the line and split by commas
  do i = 1, len_trim(line_body)
      if (line_body(i:i) == ',') then
          num_values = num_values + 1
          line_array(num_values) = line_body(start:i-1)
          start = i + 1
      end if
  end do

  ! Capture the last value after the final comma
  if (start <= len_trim(line_body)) then
      num_values = num_values + 1
      line_array(num_values) = line_body(start:)
  end if

  ! Convert the strings in line_array to atom_charge
  do j = 1, num_values
      read(line_array(j), *, iostat=i) temp_real
      if (i == 0) then
          line_array_out(j) = temp_real
      else
          write(*,*) "Error converting string to real at index ", j
      end if
  end do
end subroutine split_line_to_real_array


SUBROUTINE read_trajectory_full_run(input_filename)
  character(len=*), intent(in) :: input_filename
  integer :: i, error_code, line_counter, atom_counter, atomic_number
  character(len=256) :: line
  character(len=3) :: atom_symbol

  ! Open the file to find number of simulations in the file
  open(unit=trajectory_input_file, file=trim(adjustl(input_filename)), status='old', action='read', iostat=error_code)
  if (error_code /= 0) then
    write(*,*) "***ERROR*** Error opening: ", input_filename
    write(*,*) "Ensure that ", input_filename, " is in your running directory"
    write(log_file,*) "***ERROR*** Error opening: ", input_filename
    write(log_file,*) "Ensure that ", input_filename, " is in your running directory"
    stop
  end if

  line_counter=0
  atom_counter=0
  do
    read(trajectory_input_file, '(a)', end=110) line
    line_counter=line_counter+1
    line=trim(adjustl(line))

    if (line_counter==1) then
      read(line,*)N_total_atom
      call allocate_and_initialize_arrays
    endif

    ! Ignore anything on the line that appears after a '#'
    i = index(line, '#')
    if (i > 0) then
        line = trim(adjustl(line(1:i-1)))
    endif

    ! Skip to the next iteration if the line is empty
    if (len(trim(line)) == 0) cycle

    ! Read in the atom types and coordinates SAM FIX THIS
    if (line_counter > 2) then
      atom_counter = atom_counter + 1
      if (atom_counter <= N_total_atom) then
        read(line, '(A1)', iostat=error_code) atom_symbol
        if (error_code == 0) then
          ! Search for the corresponding atomic number
          atomic_number = 0
          do i = 1, N_elements
            if (trim(adjustl(atom_symbol)) == trim(adjustl(element_symbols(i)))) then
                atomic_number = i
                exit
            endif
          enddo
          if (atomic_number > 0) then
            atom_atomic_number(atom_counter) = atomic_number
          else
            print *, "Unknown atom symbol: ", atom_symbol
          endif
        endif
      endif
    endif

    ! Reset atom counter at the start of each new block
    if (atom_counter == N_total_atom) then
        exit
    endif
  enddo
  110 close(trajectory_input_file)

  ! Error check
  if (size(atom_atomic_number) > atom_counter) then
    write(*,*) "***ERROR*** Less atoms found (atom_counter) than are in atom_atomic_number from: ", input_filename
    write(log_file,*) "***ERROR*** Less atoms found (atom_counter) than are in atom_atomic_number from: ", input_filename
    stop
  end if
  write(log_file,*) "Successfully read and initialized data from: ", input_filename

END SUBROUTINE read_trajectory_full_run


SUBROUTINE find_atom_kinetics_at_t(input_filename,j)
  character(len=*), intent(in) :: input_filename
  integer, intent(in) :: j
  integer :: file, i, error_code, line_counter, atom_counter,iter_num,next_iter,found_iter
  real*8 :: x,y,z
  character(len=256) :: line
  character(len=3) :: atom_symbol
  logical :: time_step_found,initial_positions_set
  real*8,dimension(3,N_total_atom) :: next_atom_position

  file=trajectory_input_file+(2*j)
  ! Open the file to find number of simulations in the file
  open(unit=file, file=trim(adjustl(input_filename)), status='old', action='read', iostat=error_code)
  if (error_code /= 0) then
    write(*,*) "***ERROR*** Error opening: ", input_filename
    write(*,*) "Ensure that ", input_filename, " is in your running directory"
    write(log_file,*) "***ERROR*** Error opening: ", input_filename
    write(log_file,*) "Ensure that ", input_filename, " is in your running directory"
    stop
  end if

  time_step_found=.FALSE.
  initial_positions_set=.FALSE.
  line_counter=0
  atom_counter=0
  iter_num=0
  found_iter=0
  next_iter=0

  do
    read(file, '(a)', end=110) line
    line_counter=line_counter+1
    line=trim(adjustl(line))

    !print*,"DEBUG: OG LINE IS:", line

    i = index(line, '#')
      if (i>0) then
        i = index(line, '=')
        if (i > 0) then
            line = trim(adjustl(line(i+1:)))
        endif

        i = index(line, " ")
        if (i > 0) then
          line = trim(adjustl(line(1:i-1)))
        endif

        read(line,*)iter_num

      endif

    ! Skip to the next iteration if the line is empty
    if (iter_num < traj_time_step_to_initialize) cycle


    if(iter_num >= traj_time_step_to_initialize) then
      if(time_step_found .and. iter_num>found_iter .and. next_iter<found_iter) then
        next_iter=iter_num !set the next iter number after the searched for iter number has been found
        cycle
      else if(.not.(time_step_found)) then
        time_step_found=.TRUE.
        found_iter=iter_num
        cycle
      endif
    endif    

    ! Read in the atom types and coordinates SAM FIX THIS
    if (time_step_found) then
      if (atom_counter <= N_total_atom .and. .not.(initial_positions_set)) then
        atom_counter = atom_counter + 1
        read(line, *, iostat=error_code) atom_symbol, x, y, z
        if (error_code == 0) then
          atom_position(1,atom_counter)=x
          atom_position(2,atom_counter)=y
          atom_position(3,atom_counter)=z
        endif
      else if (.not.(initial_positions_set)) then
        initial_positions_set=.TRUE.
      endif
      
      ! All of the checks to make sure it is at the next iteration after the found one
      if (next_iter > traj_time_step_to_initialize .and. &
      atom_counter > N_total_atom .and. &
      atom_counter <= (2*N_total_atom+1) .and. &
      initial_positions_set) then
        atom_counter = atom_counter + 1
        read(line, *, iostat=error_code) atom_symbol, x, y, z
        if (error_code == 0) then
          next_atom_position(1,atom_counter-N_total_atom-1)=x
          next_atom_position(2,atom_counter-N_total_atom-1)=y
          next_atom_position(3,atom_counter-N_total_atom-1)=z
        endif
      endif
    endif

    ! exit the loop after atom positions have been placed for set position and next position
    if (atom_counter >= 2*N_total_atom+1) then
        exit
    endif
  enddo
  110 close(file)

  ! do i=1, N_total_atom
  !   print*, "DEBUG: atom_position=",atom_position(:,i)
  ! end do
  ! do i=1, N_total_atom
  !   print*, "DEBUG: next_atom_position=",next_atom_position(:,i)
  ! end do
  do i=1, N_total_atom
    atom_velocity(:,i) = (next_atom_position(:,i) - atom_position(:,i)) / ((next_iter-found_iter)*time_step)
  end do
  ! do i=1, N_total_atom
  !   print*, "DEBUG: atom_velocity=",atom_velocity(:,i)
  ! end do
END SUBROUTINE find_atom_kinetics_at_t


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
  if (len_trim(moleculeformations_input_path) > 0 .and. moleculeformations_input_path(1:1) /= ".") then
    if (moleculeformations_input_path(1:1) /= '/') then
      moleculeformations_input_path = '/' // trim(adjustl(moleculeformations_input_path))
    endif
    if (moleculeformations_input_path(len_trim(moleculeformations_input_path):len_trim(moleculeformations_input_path)) /= '/') then
      moleculeformations_input_path = trim(adjustl(moleculeformations_input_path)) // '/'
    endif
  endif
  if (len_trim(full_runs_dir_input_path) > 0 .and. full_runs_dir_input_path(1:1) /= ".") then
    if (full_runs_dir_input_path(1:1) /= '/') then
      full_runs_dir_input_path = '/' // trim(adjustl(full_runs_dir_input_path))
    endif
    if (full_runs_dir_input_path(len_trim(full_runs_dir_input_path):len_trim(full_runs_dir_input_path)) /= '/') then
      full_runs_dir_input_path = trim(adjustl(full_runs_dir_input_path)) // '/'
    endif
  endif
  if (len_trim(fragment_input_path) > 0 .and. fragment_input_path(1:1) /= ".") then
    if (fragment_input_path(1:1) /= '/') then
      fragment_input_path = '/' // trim(adjustl(fragment_input_path))
    endif
    if (fragment_input_path(len_trim(fragment_input_path):len_trim(fragment_input_path)) /= '/') then
      fragment_input_path = trim(adjustl(fragment_input_path)) // '/'
    endif
  endif
  if (len_trim(velocity_input_path) > 0 .and. velocity_input_path(1:1) /= ".") then
    if (velocity_input_path(1:1) /= '/') then
      velocity_input_path = '/' // trim(adjustl(velocity_input_path))
    endif
    if (velocity_input_path(len_trim(velocity_input_path):len_trim(velocity_input_path)) /= '/') then
      velocity_input_path = trim(adjustl(velocity_input_path)) // '/'
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
  trajectory_filename=trim(adjustl(trajectory_directory))//"/"//trim(adjustl(bare_trajectory_filename))
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

  write(log_file,*) 'Element, Atom atomic number'
  do i=1, N_total_atom
    elem=atom_atomic_number(i)
    write(log_File,*) element_symbols(elem), elem
  end do

  write(log_file,*) 'Atom positions (A) in x,y,z:'
  do i=1, N_total_atom
    write(log_file,*) atom_initial_position(:,i)
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
  
SUBROUTINE update_atom_info_file(j)
  !This subroutine updates csv-file with info about each of the atoms
  !from each of the simulations. 
  integer, intent(in), optional :: j

  integer :: i, error_code, atom_zero_index
  character(len=255) :: atom_zero_index_string
  real*8 :: charge_sum

  ! write the atom info header
  if (run_type==1) then
    WRITE(atom_info_file, '(A, I0, A)', ADVANCE='NO') "r_seed=", ion_velocity_init_seed, ", "
  endif
  if (run_type==2) then
    WRITE(atom_info_file, '(A, A)', ADVANCE='NO') trim(adjustl(full_runs_array(j))), ", "
  endif
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
    if (finalize_info_at_cap) then
      write(atom_info_file, '(F10.6, A)', ADVANCE='NO') times_at_cap(i), ", "
    else
      write(atom_info_file, '(F10.6, A)', ADVANCE='NO') time, ", "
    endif
  end do
  write(atom_info_file,*) ""

  ! write charge Sum
  charge_sum = SUM(atom_charge)
  write(atom_info_file, '(A, F15.10, A)', ADVANCE='NO') "Charge Sum, ", charge_sum, ", "
  write(atom_info_file,*)

  ! write X Velocity
  write(atom_info_file, '(A)', ADVANCE='NO') "X Velocity[A/fs], "
  do i = 1, N_total_atom
    if (finalize_info_at_cap) then
      write(atom_info_file, '(F10.6, A)', ADVANCE='NO') atom_velocity_final(1, i), ", "
    else
      write(atom_info_file, '(F10.6, A)', ADVANCE='NO') atom_velocity(1, i), ", "
    endif
  END do
  write(atom_info_file,*) ""

  ! write Y Velocity
  write(atom_info_file, '(A)', ADVANCE='NO') "Y Velocity[A/fs], "
  do i = 1, N_total_atom
    if (finalize_info_at_cap) then
      write(atom_info_file, '(F10.6, A)', ADVANCE='NO') atom_velocity_final(2, i), ", "
    else
      write(atom_info_file, '(F10.6, A)', ADVANCE='NO') atom_velocity(2, i), ", "
    endif
  END do
  write(atom_info_file, '(A)') ""

  ! write Z Velocity
  write(atom_info_file, '(A)', ADVANCE='NO') "Z Velocity[A/fs], "
  do i = 1, N_total_atom
    if (finalize_info_at_cap) then
      write(atom_info_file, '(F10.6, A)', ADVANCE='NO') atom_velocity_final(3, i), ", "
    else
      write(atom_info_file, '(F10.6, A)', ADVANCE='NO') atom_velocity(3, i), ", "
    endif
  END do
  write(atom_info_file,*) ""

  ! write Speed
  write(atom_info_file, '(A)', ADVANCE='NO') "Speed[A/fs], "
  do i = 1, N_total_atom
    if (finalize_info_at_cap) then
      write(atom_info_file, '(F10.6, A)', ADVANCE='NO') SQRT(SUM(atom_velocity_final(:, i)**2)), ", "
    else
      write(atom_info_file, '(F10.6, A)', ADVANCE='NO') SQRT(SUM(atom_velocity(:, i)**2)), ", "
    endif
  end do
  write(atom_info_file,*) ""
  write(atom_info_file,*) ""

  write(log_file,*) "Successfully updated atom info file: ", atom_info_filename


END SUBROUTINE update_atom_info_file

END MODULE MOD_IO
