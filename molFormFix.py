import numpy as np
import pandas as pd

def read_trajectory_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    trajectory_data = {}
    num_atoms = int(lines[0].strip())
    for i in range(0, len(lines), num_atoms + 2):
        time_line = lines[i + 1].strip()
        time_fs = float(time_line.split('=')[1].split()[0]) / 1000
        positions = []
        for j in range(2, num_atoms + 2):
            parts = lines[i + j].split()
            positions.append([float(parts[1]), float(parts[2]), float(parts[3])])
        trajectory_data[time_fs] = np.array(positions)
    
    return trajectory_data, num_atoms

def calculate_velocity(positions_a,positions_b, dt):
    velocities = (positions_a - positions_b) / dt
    return velocities

def main():
    input_path = 'data/isoxazole_quantum/'
    input_csv = input_path + 'moleculeFormations_14.csv'
    output_csv = input_path + 'moleculeFormations_14_new.csv'
    
    df = pd.read_csv(input_csv, header=None)
    new_data = []
    
    for i in range(0, len(df), 8):
        sim_name = df.iloc[i, 0].split(',')[0]
        trajectory_file = input_path + f'trajectory/trajectory_r{sim_name.split("_")[1].split("r")[1]}.xyz'
        trajectory_data, num_atoms = read_trajectory_file(trajectory_file)
        
        dt = 0.5  # Assuming a time step of 0.5 fs, adjust as needed
        times = df.iloc[i + 2, 1:].values
        times = times[times != ' '].astype(float)
        velocities = []
        if i==0:
            print("sim_name",sim_name)
            print("times:", times)
            print(f"trajectory data at times[0]={times[0]}:\n",trajectory_data[times[0]])
        
        atom_count=0
        for t in times:
            if t in trajectory_data:
                positions_a = trajectory_data[t][atom_count]
                positions_b = trajectory_data[t - dt][atom_count]
                vel = calculate_velocity(positions_a,positions_b, dt)
                velocities.append(vel)
                if i==0:
                    print("vel",vel)
            else:
                velocities.append(np.zeros((num_atoms, 3)))
            atom_count+=1
        
        x_velocities = [vel[0] for vel in velocities]
        y_velocities = [vel[1] for vel in velocities]
        z_velocities = [vel[2] for vel in velocities]
        speeds = [np.linalg.norm(vel) for vel in velocities]
        print("speeds",speeds)
        
        new_data.extend(df.iloc[i:i+4].values.tolist())
        new_data.append(['X Velocity[A/fs]'] + [f'{vel:.6f}' for vel in x_velocities])
        new_data.append(['Y Velocity[A/fs]'] + [f'{vel:.6f}' for vel in y_velocities])
        new_data.append(['Z Velocity[A/fs]'] + [f'{vel:.6f}' for vel in z_velocities])
        new_data.append(['Speed[A/fs]'] + [f'{vel:.6f}' for vel in speeds])
        new_data.append([])
    
    new_df = pd.DataFrame(new_data)
    new_df.to_csv(output_csv, header=False, index=False)

if __name__ == "__main__":
    main()