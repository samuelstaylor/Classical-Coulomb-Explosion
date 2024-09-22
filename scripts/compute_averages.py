import csv

# Initialize sums and counters for densities and speeds
density_sum = {'C[0]': 0, 'C[1]': 0, 'H[2]': 0, 'H[3]': 0}
speed_sum = {'C[0]': 0, 'C[1]': 0, 'H[2]': 0, 'H[3]': 0}
count = 0

# Read the CSV file
with open('scripts\\tddft_output\\\moleculeFormations_14.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        if 'Densities' in row:
            # Extract densities for C[0], C[1], H[2], H[3]
            density_sum['C[0]'] += float(row[1])
            density_sum['C[1]'] += float(row[2])
            density_sum['H[2]'] += float(row[3])
            density_sum['H[3]'] += float(row[4])
            count += 1
        elif 'Speed[A/fs]' in row:
            # Extract speeds for C[0], C[1], H[2], H[3]
            speed_sum['C[0]'] += float(row[1])
            speed_sum['C[1]'] += float(row[2])
            speed_sum['H[2]'] += float(row[3])
            speed_sum['H[3]'] += float(row[4])

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

'''
8/27 program ran... Here are the results:

Pulling from 200 simulations...

Average Densities:
C[0]: 1.9346947610792389
C[1]: 1.7448427843593408
H[2]: 0.0008513266054706691
H[3]: 7.47648574243741e-05

Average Speeds:
C[0]: 0.15821626000000008
C[1]: 0.15869406
H[2]: 0.8019361199999993
H[3]: 0.8087323350000003
'''
print()
print("Charge of C[0]=", 4-1.9346947610792389)
print("Charge of C[1]=", 4-1.7448427843593408)
print("Charge of H[2]=", 1-0.0008513266054706691)
print("Charge of H[3]=", 1-7.47648574243741e-05)

