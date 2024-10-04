# FILEPATH: /c:/Users/Samuel/Research/varga/fortran/Classical-Coulomb-Explosion/scripts/C2H6_r13_trajectory/remove_h_from_traj.py

input_file = "scripts/C2H6_r13_trajectory/trajectory.xyz"
output_file = "scripts/C2H6_r13_trajectory/remove_h_traj.xyz"

with open(input_file, "r") as f_in:
    lines = f_in.readlines()

with open(output_file, "w") as f_out:
    for line in lines:
        temp_line = line.strip().lower()
        if not temp_line.startswith("h"):
            f_out.write(line.strip() + "\n")
