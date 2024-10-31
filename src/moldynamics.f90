! Module that contains all of the subroutines for molecule dynamics computations

MODULE MOLDYNAMICS
  USE MOD_GLOBAL
  USE MOD_IO
  implicit none
  
CONTAINS

SUBROUTINE compute_atomic_masses
  !This subroutine computes the atomic masses of all atoms present in
  !the system
  integer :: i, elem

  ! Calculates assuming a neutral charge and most common isotope
  do i=1,N_total_atom
    elem=atom_atomic_number(i)
    if (elem/=0) then
      if (use_average_atomic_mass) then
        atom_mass(i)=element_num_nucleons(elem)*mass_convfactor
      else
        atom_mass(i)=element_average_atomic_mass(elem)*mass_convfactor
      endif
    endif
  enddo

  ! Removes (or adds) the electrons that are in each ion from the neutral mass
  if (include_electron_mass) then
    do i=1, N_total_atom
      elem=atom_atomic_number(i)
      if (elem/=0) then
        atom_mass(i)=atom_mass(i)-((elem-atom_charge(i))*ElectronMass)
      end if
    end do
  end if

  write(log_file,*) "All atomic masses successfully computed:"
  call print_masses_to_log

END SUBROUTINE compute_atomic_masses


SUBROUTINE calculate_force
  ! This subroutine calculates the force between ions and fills the atom_force matrix 
  ! In SI units Newtons (N)
    integer :: i, j
    real*8 :: r(3), f_vec(3), r_mag, w, pulse_force(3)
        
    ! Initialize atom_force array to zero
    atom_force=0.0

    j=1
    do i=1, N_total_atom
      do while (j < i)
        ! Calculate the distance vector and its magnitude
        r = atom_position(:, j) - atom_position(:, i)
        r_mag = SQRT(SUM(r**2))
      
        ! Coulomb's law in next following two lines of code:
        w=atom_charge(i)*atom_charge(i)/r_mag**3*E2

        ! Calculate the force vector
        f_vec = w * r 

        ! Add the force contribution to the i-th atom
        ! Force experienced on atom j from atom i
        atom_force(:, j) = atom_force(:, j) + f_vec
        
        ! Equal and opposite reaction: F on i from j
        atom_force(:, i) = atom_force(:, i) - f_vec
        
        j=j+1
      end do
      j=1
      ! Include the external electric field pulse via the dipole approximation
      if (include_pulse) then
        pulse_force=atom_charge(i)*pulse_array(iter)
        atom_force(:, i) = atom_force(:, i) + pulse_force

      endif
    end do  

    call calculate_acceleration
END SUBROUTINE calculate_force


SUBROUTINE calculate_acceleration
    integer :: i

    do i=1, N_total_atom
      atom_acceleration(:,i) = atom_force(:,i) / atom_mass(i)
    end do
END SUBROUTINE


SUBROUTINE calculate_position
    integer :: i

    do i=1, N_total_atom
      atom_position(:,i) = atom_position(:,i) + atom_velocity(:,i) * time_step + 0.5d0 * atom_acceleration(:,i) * time_step**2
    end do
END SUBROUTINE calculate_position


SUBROUTINE calculate_velocity
  ! This subroutine calcultes the velocity of each of the atoms and fills the atom_velocity matrix
  ! In units A/fs
    integer :: i
    
    do i=1, N_total_atom
      atom_velocity(:,i) = atom_velocity(:,i) + atom_acceleration(:,i) * time_step
    end do
END SUBROUTINE calculate_velocity


END MODULE MOLDYNAMICS