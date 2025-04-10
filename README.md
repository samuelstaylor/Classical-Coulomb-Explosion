# Classical-Coulomb-Explosion
Program to simulate Coulomb explosion of molecules and generate output files.

**Description:** Simulates coulomb explosion of molecules given their atoms' element, initial position, and charge. 
Initializes the atoms' velocities to the Boltzmann distribution. Can perform an number of simulations in a given run
with varying seeds for the initialized velocities. Neglects quantum effects and only takes into account the Coulomb
force for the molecular dynamics. Used for Varga Groups TDDFT simulation to compare Coulomb explosion with and without
quantum effects.

**Program Input:** The program expects three different input files:
        
        control.inp file: Specifies the run conditions for the program. Key words to set in this file:
                run_type - The type of calculation perform: 1 (classical) or 2 (semiclassical)
                N_simulations - The number of Coulomb Explosion simulations to run and generate output for.
                N_time_steps - The number of time steps in each of the simulations.
                time_step - The time step of each simulation (in fs).
                trajectory_output_frequency - The frequency to output position data to the trajectory.xyz file. 
                use_average_atomic_mass - Use the average number of nucleons for the atom mass computation
                include_electron_mass - Logical: Include the electron difference in the mass computation.
                output_trajectory - Logical: Generate trajectory output files.
                output_atom_info - Logical: Generate atom info file for the run.

For Run-type 1 (classical):

        molecule.inp file: Must be in the same directory that you are in when you run the program. Needed in order to access the number of atoms, species, and positions of atoms in a molecule, and charge of each atom. 
        
        seeds.inp file: A list of seeds for each Coulomb explosion simulation.

        temperature_ions - Temperature of the molecule in kelvin


For Run-type 2 (semi-classical):

        moleculeformations_input_path - Path to pull molecule formations file from

        full_runs_dir_input_path - Path to all of the TDDFT runs

        traj_time_step_to_initialize= Iteration number to resume TDDFT simulations classically

**Program Output:** The program will output the following files:
        
        monitor.out file: The log file that tracks the status of the program and any errors

        all_variable_file.txt file: All the variables used in that particular run of the program.

        atom_info.csv file: Useful information about each of the atoms.

        trajectory_r<rval>.xyz file: Positions of the atoms at different time increments.

**Running the Program:** The program must be compiled in order to create an executable. There are two methods of doing this: either via the Makefile or manually compiling and then linking each file. See the .txt files included in the 'instruction' directory. Ensure that all of the input files are in your current directory when you are running the program. Run with the following command:

        ./<program_executable>


**Acknowledgements:** Kalman Varga, Cody Covington
