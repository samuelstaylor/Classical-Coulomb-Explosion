import csv
import matplotlib.pyplot as plt
import matplotlib as mpl

# Function to read speeds from a CSV file
def read_speeds(file_path):
    speeds = {'C[0]': [], 'C[1]': [], 'H[2]': [], 'H[3]': []}
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if 'Speed[A/fs]' in row:
                speeds['C[0]'].append(float(row[1]))
                speeds['C[1]'].append(float(row[2]))
                speeds['H[2]'].append(float(row[3]))
                speeds['H[3]'].append(float(row[4]))
    return speeds

# Set the font to DejaVu Sans (consistent font for text, subscripts, and numbers)
mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['axes.labelweight'] = 'bold'
mpl.rcParams['font.size'] = 12

# Read speeds from both quantum and classical CSV files

# run_type=1 settings
'''
quantum_speeds = read_speeds('scripts/C2H2_run_type1/moleculeFormations_14.csv')
classical_speeds = read_speeds('scripts/C2H2_run_type1/atom_info.csv')
'''

# run_type=2 settings

quantum_speeds = read_speeds('scripts/C2H2_run_type1/moleculeFormations_14.csv')
classical_speeds = read_speeds('aces_output/atom_info.csv')


# Function to format labels with subscripts
def format_label(atom_label):
    if '[' in atom_label and ']' in atom_label:
        base, subscript = atom_label.split('[')
        subscript = subscript.replace(']', '')
        return f'{base}$_{subscript}$'
    return atom_label

# Plot the distributions
atoms = ['C[0]', 'C[1]', 'H[2]', 'H[3]']
plt.figure(figsize=(12, 8))

for i, atom in enumerate(atoms, 1):
    plt.subplot(2, 2, i)
    plt.hist(quantum_speeds[atom], bins=30, alpha=0.5, label='Quantum', color='blue')
    plt.hist(classical_speeds[atom], bins=30, alpha=0.5, label='Classical', color='orange')
    formatted_label = format_label(atom)
    plt.title(f'Speed Distribution for {formatted_label}', weight='bold')
    plt.xlabel('Speed [Å/fs]', weight='bold')
    plt.ylabel('Frequency', weight='bold')
    plt.legend()

plt.tight_layout()
plt.savefig(f'images/speed_histogram.png')
plt.show()
