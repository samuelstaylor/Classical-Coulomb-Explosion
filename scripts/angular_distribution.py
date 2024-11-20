import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Set global font to Times New Roman
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['axes.labelweight'] = 'bold'
mpl.rcParams['axes.labelsize'] = 24  # Increase label font size
mpl.rcParams['xtick.labelsize'] = 16  # Increase x tick font size
mpl.rcParams['ytick.labelsize'] = 16  # Increase y tick font size
mpl.rcParams['legend.fontsize'] = 16  # Legend font size

mass_C = 1243.7123176930008  # Mass of C atom in nano units (eV_fs^2/A^2)
mass_H = 103.64269314108340  # Mass of H atom in nano units (eV_fs^2/A^2)

class AngularDistribution:
    def __init__(self, file_path, element, atom_number, data_type, alpha):
        self.file_path = file_path
        self.thetas = []
        self.kinetic_energies = []
        self.alpha = alpha
        if element.strip().lower().startswith('c'):
            self.element = "C"
            self.mass = mass_C
            self.color = 'b'
        elif element.strip().lower().startswith('h'):
            self.element = "H"
            self.mass = mass_H
            self.color = 'r'
        self.atom_number = atom_number
        self.atom_number_sub = atom_number
        self.atom_number_to_subscript()
        if data_type.strip().lower().startswith('c'):  # classical
            self.data_type = 'Classical'
        if data_type.strip().lower().startswith('s'):  # semi-classical
            self.data_type = 'Semi-classical'
        if data_type.strip().lower().startswith('q'):  # quantum
            self.data_type = 'Quantum'
        self.SHOW_LEGENDS = False

    def atom_number_to_subscript(self):
        subscripts = ['₁', '₂', '₃', '₄']
        self.atom_number_sub = subscripts[self.atom_number]

    def calculate_theta(self, x_vel, y_vel, z_vel):
        speed = np.sqrt(x_vel ** 2 + y_vel ** 2 + z_vel ** 2)
        theta = np.degrees(np.arccos(x_vel / speed))  # Convert to degrees
        if x_vel < 0 and y_vel < 0:
            theta = 180 + (180 - theta)
        elif y_vel < 0:
            theta = -theta
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
        ax.scatter(self.thetas, self.kinetic_energies, marker='o', color=self.color,
                   label=f'{self.element}{self.atom_number_sub}', alpha=self.alpha)
        ax.grid(True)
        if self.SHOW_LEGENDS:
            ax.legend()
        if self.thetas[1] > -30 and self.thetas[1] < 30:
            ax.set_xlim(-15, 15)
        elif self.thetas[1] > 150 and self.thetas[1] < 210:
            ax.set_xlim(165, 195)
        else:
            ax.set_xlim(-180, 180)
        if self.element == "C":
            ax.set_ylim(10, 28)
        if self.element == "H":
            ax.set_ylim(18, 36)

# Create an instance of the AngularDistribution class for all atoms
def main():
    atom_types = ['C', 'C', 'H', 'H']
    
    file_path = 'data\\c2h2_quantum\\moleculeFormations_14.csv'
    file_path = 'data\\c2h2_classical\\atom_info.csv'
    file_path = 'data\\c2h2_semi_classical_pulse\\atom_info.csv'
    file_path = 'data\\c2h2_semi_classical\\atom_info.csv'


    data_type = 'q'
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    alpha = 0.5

    for atom_number, atom_type in enumerate(atom_types):
        angular_distribution = AngularDistribution(file_path=file_path,
                                                   element=atom_type,
                                                   atom_number=atom_number,
                                                   data_type=data_type,
                                                   alpha=alpha)
        angular_distribution.parse_data()
        
        row, col = divmod(atom_number, 2)
        ax = axes[row, col]
        angular_distribution.plot(ax)

        # Remove x-axis ticks and labels for the top row
        if row == 0:
            ax.set_xticklabels([])
            ax.set_xlabel('')  # Remove x-axis label

        # Remove y-axis ticks and labels for the right column
        if col == 1:
            ax.set_yticklabels([])
            ax.set_ylabel('')  # Remove y-axis label

    # Add global x and y labels
    # Add global x and y labels with adjusted positions
    fig.text(0.54, 0.04, r'θ°', ha='center', va='center', fontweight='bold', fontsize=24)  # Move x-axis label slightly to the right
    fig.text(0.04, 0.54, 'Kinetic Energy (eV)', ha='center', va='center', rotation='vertical', fontweight='bold', fontsize=24)  # Move y-axis label slightly up

    plt.tight_layout(rect=[0.05, 0.05, 1, 1])  # Adjust layout to make room for global labels
    plt.savefig('images/angular_distribution.png')
    plt.show()

if __name__ == '__main__':
    main()
