import os, math

Lmax = 200

def calc_dist(coord, xyz=[0,0,200], kxyz=[0,0,1]):
    dist = 0
    for i in range(3):
        dist += (coord[i]-xyz[i])**2*kxyz[i]
    return math.sqrt(dist)


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

# Key here !!! Parametrix bijective mapping from [0:1] to [0:1], where the parameter c determines the location of sharp peak
# Therefore the coordniates of atoms affect parameter c
def custom_map(x, c, a=20): # [0:1] -> [0:1]
    f  = 1.0/(1.0 + math.exp(a*(c-x)))
    f0 = 1.0/(1.0 + math.exp(a*c))
    f1 = 1.0/(1.0 + math.exp(a*(c-1)))
    g  = (f - f0)/(f1 - f0)
    return g

def substitute_coordinates(file_name, rank, num, coord0, coord1, L, neb=False):
    with open(file_name, 'r') as file:
        lines = file.readlines()

    with open(file_name, 'w') as file:
        in_atoms_section = False
        for line in lines:
            if "Atoms " in line:
                in_atoms_section = True
                if(neb==False):
                    file.write(line)
                continue
            if in_atoms_section:
                parts = line.split()
                if(len(parts)<1):
                    if(neb==False):
                        file.write("\n")
                    continue
                atom_id = int(parts[0])

                t = rank/(num-1)
                dist = calc_dist(coord0[atom_id]) 
                coeff = custom_map(t, dist/Lmax)

                for j in range(3):
                    dx = coord1[atom_id][j] - coord0[atom_id][j]
                    dx = dx - L[j] * round(dx / L[j])
                    x0 = coord0[atom_id][j]
                    x  = x0 + dx*coeff
                    parts[j+2] = "{:16.10f}".format(x)
                
                if neb:
                    parts[1] = ''
                file.write("\t" + "\t\t".join(parts) + "\n")
            else:
                if neb:
                    if " atoms" in line:
                        parts = line.split()
                        file.write(parts[0] + "\n")
                else:
                    file.write(line)

coord_init,  L  = extract_coordinates('init.lmp')
coord_final, _  = extract_coordinates('final.lmp')

N = 64
os.system('cp init.lmp 0.lmp')
os.system('cp final.lmp '+str(N-1)+'.lmp')

for i in range(N):
    os.system('cp final.lmp '+str(i)+'.lmp')
    substitute_coordinates(str(i)+'.lmp', i, N, coord_init, coord_final, L, neb=True)
