import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib.ticker import FuncFormatter
import os
import matplotlib


# Set global font to Times New Roman
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['axes.labelweight'] = 'bold'
mpl.rcParams['axes.labelsize'] = 16  # increase label font size
mpl.rcParams['xtick.labelsize'] = 14  # increase x tick font size
mpl.rcParams['ytick.labelsize'] = 14  # increase y tick font size
matplotlib.rc('font', family='Times New Roman')  # set default font to Times New Roman

def read_pre_trajectory(file):
    with open(file, 'r') as f:
        lines = f.readlines()

    iterations = []
    positions1 = []
    positions2 = []
    distances = []
    # step is 4 because there are two molecules
    for i in range(0, len(lines), 4):
        iterations.append(int(lines[i+1].split('=')[1].split()[0]))
        atom1 = np.array(list(map(float, lines[i+2].split()[1:])))
        atom2 = np.array(list(map(float, lines[i+3].split()[1:])))
        positions1.append(atom1)
        positions2.append(atom2)
        distance = np.linalg.norm(atom1-atom2)
        distances.append(distance)

    return iterations, positions1, positions2, distances

def read_trajectory(file,iterations=[],positions1=[],positions2=[],distances=[]):
    with open(file, 'r') as f:
        lines = f.readlines()
        
    if iterations==[]:
        iterations = []
    if positions1==[]:
        positions1 = []
    if positions2==[]:
        positions2 = []
    if distances==[]:
        distances = []

    # step is 4 because there are two molecules
    for i in range(0, len(lines), 4):
        iterations.append(int(lines[i+1].split('=')[1].split()[0]))
        atom1 = np.array(list(map(float, lines[i+2].split()[1:])))
        atom2 = np.array(list(map(float, lines[i+3].split()[1:])))
        positions1.append(atom1)
        positions2.append(atom2)
        distance = np.linalg.norm(atom1-atom2)
        distances.append(distance)

    return iterations, positions1, positions2, distances

def calculate_velocity(positions, time):
    velocities = np.diff(positions, axis=0) / np.diff(time)[:, None]
    speeds = np.linalg.norm(velocities, axis=1)  # calculate speed from velocity
    return speeds

def calculate_acceleration(speeds, time):
    accelerations = np.diff(speeds) / np.diff(time[1:])  # calculate acceleration from speed
    return accelerations
    
def plot_distance(time, distances, directory):
    plt.plot(time, distances, color='black')
    plt.xlabel('Time (fs)')
    plt.ylabel('Distance (Å)')
    plt.grid(True)
    plt.tight_layout()  # Adjust layout to prevent overlapping content
    plt.savefig(os.path.join(directory, 'trajectory.png'))
    plt.show()
    plt.close()

def plot_positions(time, position1, position2, directory, label1, label2):
    distance1 = np.linalg.norm(position1, axis=1)  # calculate distance from position vector
    distance2 = np.linalg.norm(position2, axis=1)  # calculate distance from position vector
    plt.plot(time, distance1, color='blue', label='$\mathrm{'+label1+'}$')
    plt.plot(time, distance2, color='green', label='$\mathrm{'+label2+'}$')
    plt.xlabel('Time (fs)')
    plt.ylabel('Distance from origin (0,0,0) (Å)')
    plt.grid(True)
    plt.tight_layout()  # Adjust layout to prevent overlapping content
    plt.legend()
    plt.savefig(os.path.join(directory, 'position.png'))
    plt.show()
    plt.close()

def plot_velocity(time, speeds1, speeds2, directory, label1, label2):
    plt.plot(time[1:], speeds1, color='blue', label='$\mathrm{'+label1+'}$')
    plt.plot(time[1:], speeds2, color='green', label='$\mathrm{'+label2+'}$')
    plt.xlabel('Time (fs)')
    plt.ylabel('Speed (Å/fs)')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()  # Adjust layout to prevent overlapping content
    plt.savefig(os.path.join(directory, 'velocity.png'))
    plt.show()
    plt.close()

def plot_acceleration(time, accelerations1, accelerations2, directory, label1, label2):
    plt.plot(time[2:], accelerations1, color='blue', label='$\mathrm{'+label1+'}$')
    plt.plot(time[2:], accelerations2, color='green', label='$\mathrm{'+label2+'}$')
    plt.xlabel('Time (fs)')
    plt.ylabel('Acceleration (Å/fs²)')
    plt.grid(True)
    plt.ticklabel_format(axis='y')  # format y-axis in scientific notation
    plt.legend()
    plt.tight_layout()  # Adjust layout to prevent overlapping content
    plt.savefig(os.path.join(directory, 'acceleration.png'))
    plt.show()
    plt.close()
    
def plot_force(time, accelerations1, accelerations2, mass1, mass2, directory, label1, label2):
    force1 = mass1 * accelerations1  # calculate force from acceleration
    force2 = mass2 * accelerations2  # calculate force from acceleration
    plt.plot(time[2:], force1, color='blue', label='$\mathrm{'+label1+'}$')
    plt.plot(time[2:], force2, color='green', label='$\mathrm{'+label2+'}$')
    plt.xlabel('Time (fs)')
    plt.ylabel('Force (eV²/Å)')  # change y-label to 'Force'
    plt.grid(True)
    plt.ticklabel_format(axis='y')  # format y-axis in scientific notation
    plt.legend()
    plt.tight_layout()  # Adjust layout to prevent overlapping content
    plt.savefig(os.path.join(directory, 'force.png'))  # change filename to 'force.png'
    plt.show()
    plt.close()

trajectory_file_1='scripts\\C2H6_r13_trajectory\\remove_h_traj.xyz'
trajectory_file_2 = 'scripts\\C2H6_fragment\\trajectory.xyz'
label1="CH_2"
label2="CH_3"
#label1="C_2H"
#label2="CH_2"
mass1=1450.9977039751675
mass2=1554.6403971162511 
iterations = []
positions1 = []
positions2 = []
distances = []
#iterations, positions1, positions2, distances = read_pre_trajectory(trajectory_file_1)
iterations, positions1, positions2, distances = read_trajectory(trajectory_file_2,iterations, positions1, positions2, distances)
directory = os.path.dirname(trajectory_file_2)
time = np.array(iterations) / 1000  # convert iterations to time in fs
speeds1 = calculate_velocity(positions1, time)
speeds2 = calculate_velocity(positions2, time)
accelerations1 = calculate_acceleration(speeds1, time)
accelerations2 = calculate_acceleration(speeds2, time)
plot_distance(time, distances, directory)
plot_positions(time, positions1, positions2, directory,label1,label2)
plot_velocity(time, speeds1, speeds2, directory,label1,label2)
plot_acceleration(time, accelerations1, accelerations2, directory,label1,label2)
plot_force(time, accelerations1, accelerations2, mass1,mass2, directory,label1,label2)