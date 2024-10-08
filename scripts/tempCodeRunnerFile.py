def plot_distance(time, distances, directory, mode='', show=False):
    plt.plot(time, distances, color='black')
    if mode.lower().startswith('s'):
        plt.axvline(x=25, color='#ee87ee', linestyle='--')  # Add vertical dashed line at time=25
    plt.xlabel('Time (fs)')
    plt.ylabel('Distance (Å)')
    plt.ylim(0, 60)  # Set y-axis limits from 0 to 60
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(directory, 'distance.png'))
    if show:
        plt.show()
    plt.close()