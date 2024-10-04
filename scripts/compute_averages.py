import csv

N_carbon   = 4
N_hydrogen = 10
input_csv_file = 'data/c4h10_quantum/28/moleculeFormations_28.csv'

# Initialize sums and counters for densities and speeds for C₄H₁₀
density_sum = {f'C[{i}]': 0 for i in range(N_carbon)}
density_sum.update({f'H[{i}]': 0 for i in range(N_carbon, N_carbon+N_hydrogen)})
speed_sum = {f'C[{i}]': 0 for i in range(N_carbon)}
speed_sum.update({f'H[{i}]': 0 for i in range(N_carbon, N_carbon+N_hydrogen)})
count = 0

# Read the CSV file
with open(input_csv_file, 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        if 'Densities' in row:
            # Extract densities for C[0] to C[3] and H[4] to H[13]
            for i in range(N_carbon):
                density_sum[f'C[{i}]'] += float(row[i + 1])
            for i in range(N_carbon, N_carbon+N_hydrogen):
                density_sum[f'H[{i}]'] += float(row[i + 1])
                count += 1
        elif 'Speed[A/fs]' in row:
            # Extract speeds for C[0] to C[3] and H[4] to H[13]
            for i in range(N_carbon):
                speed_sum[f'C[{i}]'] += float(row[i + 1])
            for i in range(N_carbon, N_carbon+N_hydrogen):
                speed_sum[f'H[{i}]'] += float(row[i + 1])

# Calculate averages
average_density = {atom: density_sum[atom] / count for atom in density_sum}
average_speed = {atom: speed_sum[atom] / count for atom in speed_sum}

# Print results
print("Pulling from", count, "simulations...\n")
print("Average Densities:")
for atom, density in average_density.items():
    print(f"{atom}: {density}")

print("\nAverage Speeds:")
for atom, speed in average_speed.items():
    print(f"{atom}: {speed}")

# Example charges (adjust based on actual molecular charges if needed)
print()
for i in range(N_carbon):
    print(f"Charge of C[{i}]=", 4 - average_density[f'C[{i}]'])
for i in range(N_carbon, N_carbon+N_hydrogen):
    print(f"Charge of H[{i}]=", 1 - average_density[f'H[{i}]'])
