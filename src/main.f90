! Program:      Classical Coulomb Explosion
!
! Description:  Simulates Coulomb explosion of molecules given the following:
!                 1.  molecule.inp -- Number of atoms and atom initial positions 
!                 2.  control.inp  -- Simulation parameters (N_simulations, N_time_steps, temperature, etc.)
!                 3.  seeds.inp    -- Initial seeds for the Boltzmann distribution based atom initial velocities
!               Outputs the following:
!                 1.  *terminal*            -- Program live status
!                 2.  monitor.out           -- Log of the program and any errors.
!                 3.  all_variable_file.txt -- All the variables used in that particular run of the program
!                 3.  trajectory.xyz        -- Positions of the atoms at different time increments
!                 4.  atom_info.csv         -- Useful information about each of the atoms.
! Author:       Samuel S. Taylor, 2024

PROGRAM MAIN
    USE COULOMB_EXPLOSION
    implicit none

    call initialize
    call calculate
    call cleanup
    
END PROGRAM MAIN