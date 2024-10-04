import csv
import matplotlib.pyplot as plt
import matplotlib as mpl

# Function to generate atom labels based on user input
def generate_atom_labels(num_carbons, num_hydrogens):
    labels = {}
    for i in range(num_carbons):
        labels[f'C[{i}]'] = []
    for i in range(num_hydrogens):
        labels[f'H[{i + num_carbons}]'] = []
    return labels

# Function to read speeds from a CSV file
def read_speeds(file_path, atom_labels):
    speeds = {label: [] for label in atom_labels}
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if 'Speed[A/fs]' in row:
                for i, label in enumerate(atom_labels):
                    speeds[label].append(float(row[i + 1]))  # Start from the second column for speeds
    return speeds

# Set the font to DejaVu Sans (consistent font for text, subscripts, and numbers)
mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['axes.labelweight'] = 'bold'
mpl.rcParams['font.size'] = 12

# User input for number of atoms
num_carbons = 4  # Adjust this for different molecules
num_hydrogens = 10  # Adjust this for different molecules

# Generate atom labels based on the number of atoms
atom_labels = generate_atom_labels(num_carbons, num_hydrogens)

# Read speeds from both quantum and classical CSV files
quantum_speeds = read_speeds('scripts/C4H10_boltzmann/moleculeFormations_14.csv', atom_labels)
classical_speeds = read_speeds('c4h10_aces_cont_from_tddft_output/atom_info.csv', atom_labels)

# Function to format labels with subscripts
def format_label(atom_label):
    if '[' in atom_label and ']' in atom_label:
        base, subscript = atom_label.split('[')
        subscript = subscript.replace(']', '')
        return f'{base}$_{subscript}$'
    return atom_label

# Plot the distributions
plt.figure(figsize=(14, 5))  # Adjust the figure size to match a 2x7 grid

num_rows = 2
num_cols = 7

for i, atom in enumerate(atom_labels, 1):
    plt.subplot(num_rows, num_cols, i)  # Create a 2x7 grid for the subplots
    plt.hist(quantum_speeds[atom], bins=30, alpha=0.5, label='Quantum', color='blue')
    plt.hist(classical_speeds[atom], bins=30, alpha=0.5, label='Classical', color='orange')
    formatted_label = format_label(atom)
    plt.title(f'{formatted_label}', weight='bold')
    plt.xlabel('Speed [Å/fs]', weight='bold')
    plt.ylabel('Frequency', weight='bold')
    plt.legend()

plt.tight_layout()
plt.savefig(f'images/speed_histogram.png')
plt.show()