! Declare the variables to be used in the program

MODULE MOD_GLOBAL
    USE PARAMETERS
    implicit none

    ! Set in molecule.inp
    integer :: N_total_atom ! The number of atoms in the simultaion
   
    ! All of these can be set in control.inp
    integer :: run_type=1 ! The type of program you want to run
    integer :: N_simulations=1 ! The number of simulations to run
    integer :: N_time_steps=120000 ! Number of timesteps (iterations) in the simulation 
    real*8  :: time_step=.001 ! Size of the time step to iterate program through
    integer :: trajectory_output_frequency=500 ! Frequency to output trajectory data
    real*8  :: temperature_ions=300 ! The temperature of the ions in K
    logical :: use_average_atomic_mass=.TRUE. ! Use mass of average number of nucleons for the mass calculation
    logical :: include_electron_mass=.FALSE. ! Add the mass of the electrons to the mass for the calculation


    ! global variables to not be set in control.inp but used in the program
    integer :: ion_velocity_init_seed
    real*8  :: time ! The time (in fs) of the simulation at each step
    integer :: iter 
    real*8  :: program_start_time
    real*8  :: program_end_time

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
    real*8,dimension(:,:),allocatable :: atom_position 
    real*8,dimension(:,:),allocatable :: atom_velocity
    real*8,dimension(:,:),allocatable :: atom_force
    real*8,dimension(:,:),allocatable :: atom_acceleration

    
END MODULE MOD_GLOBAL
