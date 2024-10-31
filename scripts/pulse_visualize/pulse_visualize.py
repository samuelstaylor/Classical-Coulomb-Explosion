import matplotlib.pyplot as plt

# Set font properties
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 12

# Read data from pulse.txt
time = []
e_field = []

with open('scripts\\pulse_visualize\\pulse.txt', 'r') as file:
    for line in file:
        columns = line.split()
        time.append(float(columns[0]))
        e_field.append(float(columns[1]))

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(time, e_field, color='tab:blue', label='E-field (V/A)')
plt.xlim(0, 50)

# Set the color for the vertical line
vertical_line_color = "#D3B3E5"
plt.axvline(x=25, color=vertical_line_color, linestyle='--')

# Add label underneath the vertical line
plt.text(25, min(e_field) - 0.072 * (max(e_field) - min(e_field)), 't*', 
         horizontalalignment='center', verticalalignment='top', 
         fontsize=16, color='black', fontweight='bold', fontstyle='italic')

# Make axes labels bold and increase font size
plt.xlabel('Time (fs)', fontweight='bold', fontsize=20)
plt.ylabel('E-field (V/Å)', fontweight='bold', fontsize=20)

# Increase the font size of the tick mark labels
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

# Remove the title
# plt.title('Pulse Visualization')

plt.grid(True)

plt.tight_layout()

# Save the figure
plt.savefig('scripts/pulse_visualize/pulse_vis.png')

plt.show()