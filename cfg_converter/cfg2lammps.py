#!/usr/bin/env python3
import sys
print('Converting {0} to Lammps data input format'.format(sys.argv[1]))
filename = sys.argv[1]

coordinates = []
element = []
H = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]

with open(filename, 'r') as file:
    lines = file.readlines()

i_atom = 1
for line in lines:
    data = line.split()
    if line.startswith('#') or line.strip() == "":
        continue
    if line.startswith('N'):
        Natom = int(data[4])
        continue
    if line.startswith('H'):
        i = int(data[0][3])
        j = int(data[0][5])
        val = float(data[2])
        H[i][j] = val
        continue
    if len(data) >= 6:
        try:
            float(data[0])
        except:
            continue
        atom = {
            'mass': float(data[0]),
            'symbol': data[1],
            's1': float(data[2]),
            's2': float(data[3]),
            's3': float(data[4])
        }

        try:
            element.index(atom['symbol'])
        except:
            element.append(atom['symbol'])
        atom['type'] = element.index(atom['symbol']) + 1
        atom['x'] = atom['s1'] * H[1][1]  +  atom['s2'] * H[2][1]  +  atom['s3'] * H[3][1]
        atom['y'] = atom['s1'] * H[1][2]  +  atom['s2'] * H[2][2]  +  atom['s3'] * H[3][2]
        atom['z'] = atom['s1'] * H[1][3]  +  atom['s2'] * H[2][3]  +  atom['s3'] * H[3][3]
        atom['id'] = i_atom
        i_atom = i_atom + 1
        coordinates.append(atom)

file = open("restart.in", "w")

file.write('Lammps data input converted from CFG format\n\n')
file.write('{0:7d} atoms\n      0 bonds\n      0 angles\n      0 dihedrals\n      0 impropers\n\n'.format(Natom))
file.write('{0:7d} atom types\n      0 bond types\n      0 angle types\n      0 dihedral types\n      0 improper types\n\n'.format(len(element)))
file.write('      0 {0:8.4f} xlo xhi\n      0 {1:8.4f} ylo yhi\n      0 {2:8.4f} zlo zhi\n\nAtoms\n\n'.format(H[1][1],H[2][2],H[3][3]))

for atom in coordinates:
    file.write('      {0:6d}   {1:3d}   {2:20.15f} {3:20.15f} {4:20.15f}  0  0  0\n'.format(atom['id'], atom['type'], atom['x'], atom['y'], atom['z']))

file.close()
