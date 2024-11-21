import numpy as np

def read_xyz(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    num_atoms = int(lines[0].strip())
    num_timesteps = len(lines) // (num_atoms + 2)
    
    positions = []
    for t in range(num_timesteps):
        timestep_positions = []
        for i in range(2, num_atoms + 2):
            parts = lines[t * (num_atoms + 2) + i].split()
            timestep_positions.append([float(parts[1]), float(parts[2]), float(parts[3])])
        positions.append(timestep_positions)
    
    return np.array(positions), num_atoms, num_timesteps

def calculate_velocities(positions, dt):
    velocities = np.zeros_like(positions)
    # Backward difference for all points
    velocities[1:] = (positions[1:] - positions[:-1]) / dt
    # For the first point, we can assume the velocity is the same as the second point
    velocities[0] = velocities[1]
    # For the last point, we can assume the velocity is the same as the second to last point

    return velocities

def write_xyz(filename, velocities, num_atoms, num_timesteps,dt):
    with open(filename, 'w') as file:
        for t in range(num_timesteps - 1):
            file.write(f"{num_atoms}\n")
            file.write(f"# time step = {(t*dt)} fs\n")
            for i in range(num_atoms):
                file.write(f"Atom {velocities[t, i, 0]:.6f} {velocities[t, i, 1]:.6f} {velocities[t, i, 2]:.6f}\n")

def main():
    input_file = 'data\\c2h2_quantum\\trajectory\\trajectory_r1.xyz'
    output_file = 'velocity_r1.xyz'
    dt = 0.5  # Time step duration, adjust as needed
    
    positions, num_atoms, num_timesteps = read_xyz(input_file)
    velocities = calculate_velocities(positions, dt)
    write_xyz(output_file, velocities, num_atoms, num_timesteps,dt)

if __name__ == "__main__":
    main()