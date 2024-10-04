import numpy as np
import matplotlib.pyplot as plt

mass_C = 1243.7123176930008  # Mass of C atom in nano units (eV_fs^2/A^2)
mass_H = 103.64269314108340  # Mass of H atom in nano units (eV_fs^2/A^2)

class AngularDistribution:
    def __init__(self, file_path, element, atom_number, data_type):
        self.file_path = file_path
        self.thetas = []
        self.kinetic_energies = []
        if element.strip().lower().startswith('c'):
            self.element = "C"
            self.mass = mass_C
            self.color = 'r'
        elif element.strip().lower().startswith('h'):
            self.element = "H"
            self.mass = mass_H
            self.color = 'b'
        self.atom_number = atom_number
        if data_type.strip().lower().startswith('c'):  # classical
            self.data_type = 'Classical'
        if data_type.strip().lower().startswith('s'):  # semi-classical
            self.data_type = 'Semi-classical'
        if data_type.strip().lower().startswith('q'):  # quantum
            self.data_type = 'Quantum'

    def calculate_theta(self, x_vel, y_vel, z_vel):
        speed = np.sqrt(x_vel ** 2 + y_vel ** 2 + z_vel ** 2)
        theta = np.arccos(x_vel / speed)  # Angle in radians
        theta = np.degrees(theta)  # Convert to degrees

        if x_vel < 0:
            # If y velocity component is negative, adjust theta
            if y_vel < 0:
                theta = 180 + (180 - theta)
        else:
            # If y velocity component is negative, adjust theta
            if y_vel < 0:
                theta = - theta

        return theta

    def parse_data(self):
        with open(self.file_path, 'r') as file:
            lines = file.readlines()

        for i, line in enumerate(lines):
            if line.startswith('r_seed') or line.startswith('C2H2'):
                x_vel = float(lines[i + 4].split(',')[self.atom_number + 1].strip())
                y_vel = float(lines[i + 5].split(',')[self.atom_number + 1].strip())
                z_vel = float(lines[i + 6].split(',')[self.atom_number + 1].strip())
                speed = float(lines[i + 7].split(',')[self.atom_number + 1].strip())

                theta = self.calculate_theta(x_vel, y_vel, z_vel)
                kinetic_energy = 0.5 * self.mass * speed ** 2

                self.thetas.append(theta)
                self.kinetic_energies.append(kinetic_energy)

    def plot(self, ax):
        # Ensure that LaTeX-style formatting is used correctly for subscripts in labels
        ax.scatter(self.thetas, self.kinetic_energies, marker='o', color=self.color,
                   label=f'{self.element}$_{{{self.atom_number}}}$')
        # Set axis labels with bold font
        ax.set_xlabel(r'θ°', fontsize=12, fontweight='bold')
        ax.set_ylabel('Kinetic Energy (eV)', fontsize=12, fontweight='bold')
        ax.grid(True)
        ax.legend()
        if self.thetas[1] > -30 and self.thetas[1] < 30:
            ax.set_xlim(-30, 30)
        elif self.thetas[1] > 150 and self.thetas[1] < 210:
            ax.set_xlim(150, 210)
        else:
            ax.set_xlim(-180, 180)

# Create an instance of the AngularDistribution class for all atoms
def main():
    atom_types = ['C', 'C', 'H', 'H']
    '''
    file_path = 'data/c2h2_aces_cont_from_tddft_output/atom_info.csv'
    data_type = "Semi-classical"
    '''
    file_path = 'scripts/C2H2_boltzmann/moleculeFormations_14.csv'
    data_type = 'Quantum'
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))  # Create a 2x2 grid for subplots
    fig.suptitle(r'Angular Distribution C$_{2}$H$_{2}$ ('+data_type+')', fontsize=16, fontweight='bold')

    for atom_number, atom_type in enumerate(atom_types):
        angular_distribution = AngularDistribution(file_path=file_path,
                                                   element=atom_type,
                                                   atom_number=atom_number,
                                                   data_type=data_type)
        angular_distribution.parse_data()
        
        row = atom_number // 2  # Determine the row (0 or 1)
        col = atom_number % 2   # Determine the column (0 or 1)
        angular_distribution.plot(axes[row, col])  # Plot on the corresponding subplot

    plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout to make room for the title
    plt.savefig('images/angular_distribution.png')
    plt.show()


if __name__ == '__main__':
    main()
