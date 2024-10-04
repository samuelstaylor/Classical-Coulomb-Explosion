! Declare the variables to be used in the program

MODULE MOD_GLOBAL
    USE PARAMETERS
    implicit none
    integer :: run_type=1

    ! Set in molecule.inp
    integer :: N_total_atom ! The number of atoms in the simultaion
   
    ! All of these can be set in control.inp
    character(len=256) :: molecule_input_path = "./"
    character(len=256) :: seeds_input_path = "./"
    character(len=256) :: moleculeformations_input_path = "./"
    character(len=256) :: full_runs_dir_input_path = "./"
    character(len=256) :: fragment_input_path = "./"
    character(len=256) :: velocity_input_path = "./"
    integer :: traj_time_step_to_initialize=0 ! Iteration number to read info from the trajectory file from 
    integer :: time_step_start=0 ! The iteration number to start the simulation at
    integer :: N_simulations=1 ! The number of simulations to run
    integer :: N_time_steps=120000 ! Number of timesteps (iterations) in the simulation 
    real*8  :: time_step=.001 ! Size of the time step to iterate program through
    integer :: trajectory_output_frequency=500 ! Frequency to output trajectory data
    real*8  :: temperature_ions=300 ! The temperature of the ions in K
    logical :: use_average_atomic_mass=.TRUE. ! Use mass of average number of nucleons for the mass calculation
    logical :: include_electron_mass=.FALSE. ! Add the mass of the electrons to the mass for the calculation
    logical :: output_trajectory=.TRUE.
    logical :: output_atom_info=.TRUE.
    logical :: parallelization=.TRUE.

    integer :: ion_velocity_init_seed=1 ! Values for this will be set in the seeds.inp file

    ! global variables to not be set in control.inp but used in the program
    real*8  :: time=0 ! The time (in fs) of the simulation at each step
    integer :: iter=0
    character(len=8) :: start_date, end_date
    character(len=8) :: start_time, end_time
    real*8  :: program_CPU_start_time
    real*8  :: program_CPU_end_time
    ! Used in the prepare_output subroutines
    character(len=15) :: formatted_datetime
    logical :: dir_exists
    character(len=255) :: output_dir_with_date_time
    character(len=255) :: trajectory_directory
    
    !full_filename
    character(len=256) :: molecule_filename_full = "molecule.inp"
    character(len=256) :: seeds_filename_full = "seeds.inp"
    character(len=256) :: moleculeformations_filename_full = "moleculeFormations.csv"
    character(len=256) :: trajectory_filename_full = "trajectory.xyz"
    character(len=256) :: fragment_filename_full = "fragment.inp"
    character(len=256) :: velocity_filename_full = "velocity.inp"


    ! file names with directories to be set
    character*255 :: log_output_filename=bare_log_output_filename
    character*255 :: all_variable_filename=bare_all_variable_filename
    character*255 :: trajectory_filename=bare_trajectory_filename
    character*255 :: atom_info_filename=bare_atom_info_filename

    ! 1-dimensional arrays that will be allocated to the size N_total_atom
    integer,dimension(:),allocatable :: atom_atomic_number
    real*8,dimension(:),allocatable :: atom_charge
    real*8,dimension(:),allocatable :: atom_mass

    ! 2-dimensional arrays that will be allocated to N_total_atoms columns
    real*8,dimension(:,:),allocatable :: atom_initial_position 
    real*8,dimension(:,:),allocatable :: atom_position 
    real*8,dimension(:,:),allocatable :: atom_velocity
    real*8,dimension(:,:),allocatable :: atom_force
    real*8,dimension(:,:),allocatable :: atom_acceleration
    real*8,dimension(:,:),allocatable :: atom_charges_every_tddft

    real*8,dimension(:),allocatable :: seed_array
    character(len=256),dimension(:),allocatable :: full_runs_array
    

    
END MODULE MOD_GLOBAL
