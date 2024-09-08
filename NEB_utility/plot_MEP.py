import os
import numpy as np
import matplotlib.pyplot as plt

N_replica = 64
cfgname = 'Cu30753_'

def extract_coordinates(file_name):
    coordinates = {}
    L = {}
    with open(file_name, 'r') as file:
        in_atoms_section = False
        for line in file:
            if "xlo" in line:
                parts = line.split()
                L[0] = float(parts[1]) - float(parts[0])
            if "ylo" in line:
                parts = line.split()
                L[1] = float(parts[1]) - float(parts[0])
            if "zlo" in line:
                parts = line.split()
                L[2] = float(parts[1]) - float(parts[0])

            if "Atoms " in line:
                in_atoms_section = True
                continue
            if in_atoms_section:
                parts = line.split()
                if(len(parts)<1):
                    continue
                atom_id = int(parts[0])
                coordinates[atom_id] = [float(parts[2]), float(parts[3]), float(parts[4])]
    return coordinates, L

def compute_disp(coord0, coord1, L, ld, dd, lim):
    # ld = layer direction, dd = direction of displacement
    disp = 0.0
    N_layer = 0
    for i in range(len(coord0)):
        atom_id = i + 1
        if(coord0[atom_id][ld] > lim):
            dx = coord1[atom_id][dd] - coord0[atom_id][dd]
            dx = dx - L[dd] * round(dx / L[dd])
            disp = disp + dx
            N_layer = N_layer + 1
    return disp/N_layer



# Function to read the last line of a file
def read_last_line(filename):
    with open(filename, 'rb') as file:
        file.seek(-2, 2)  # Jump to the second last byte
        while file.read(1) != b'\n':  # Until EOL is found
            file.seek(-2, 1)  # Keep jumping back by 2 bytes
        last_line = file.readline().decode()
    return last_line

# Function to reshape data after the 10th column into a 64x2 matrix
def reshape_data(last_line):
    # Split the line into columns (assuming whitespace delimiter)
    columns = last_line.split()
    
    # Extract the data after the 10th column
    data = columns[9:]
    
    # Convert the data into a numpy array of floats
    data_array = np.array(data, dtype=float)
    
    # Reshape the array into a 64x2 matrix
    reshaped_matrix = data_array.reshape(N_replica, 2)

    # Subtract the value in the first row, second column from all values in the second column
    reshaped_matrix[:, 1] -= reshaped_matrix[0, 1]    
    return reshaped_matrix

# Function to plot the data
def plot_data(matrix):
    x = matrix[:, 0]  # First column
    y = matrix[:, 1]  # Second column
    
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, marker='o', linestyle='-', color='b')
    plt.title('Minimum Energy Path')
    plt.xlabel('Reaction Coordinate')
    plt.ylabel('Energy (eV)')
    plt.show()

# Example usage
filename = 'log.lammps'
last_line = read_last_line(filename)
reshaped_matrix = reshape_data(last_line)

for i in range(N_replica):
    fn = cfgname+str(i)+'.cfg'
    os.system('rm tmp.lmp')
    os.system('atomsk '+fn+' tmp.lmp')
    if(i==0):
        coord0, L = extract_coordinates('tmp.lmp')
    else:
        coord1, _ = extract_coordinates('tmp.lmp')
        disp = compute_disp(coord0, coord1, L, 1, 2, 80.2)
        reshaped_matrix[i][0] = disp

    

with open('MEP', 'w') as file:
    for i in range(len(reshaped_matrix)):
        print(reshaped_matrix[i][0], reshaped_matrix[i][1])
        file.write("{:8.6f} {:6.3f}\n".format(reshaped_matrix[i][0], reshaped_matrix[i][1]))

