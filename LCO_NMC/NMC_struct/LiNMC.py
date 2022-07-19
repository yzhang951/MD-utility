# Yin Zhang, 6/8/2020
# generate NMC strucuture for LAMMPS
# crystal structure from LiNiMnCoO2.cif, from http://www.crystallography.net/
import numpy as np
import math
import random
Nx = 1
Ny = 80
Nz = 5
ratio = 0.1
fp = open("LiNMC.data","w")
Natom = 0
Nmol = 0

a = 2.868
b = a
c = 14.213

va = np.array([a*math.cos(math.radians(30)),a*math.sin(math.radians(30)),0])
vb = np.array([b*math.cos(math.radians(150)),b*math.sin(math.radians(150)),0])
vc = np.array([0,0,c])

# new translation vectors
vx = va-vb
vy = va+vb
vz = vc

Lx = Nx*vx[0]
Ly = Ny*vy[1]
Lz = Nz*vz[2]

box = np.array([Lx,Ly,Lz])
Li = [] # Li site
M  = [] # Transition metal site
O  = [] # O site
Cell = [] # Unit cell
bonds = []
H = np.array([vx[0],vy[1],vz[2]])

Li.append(np.array([0,0,0]))
Li.append((va*2/3+vb*1/3+vc*1/3) % H)
Li.append((va*1/3+vb*2/3+vc*2/3) % H)
Li.append((Li[0]+va+np.array([0,0,-0.0001])) % H)
Li.append((Li[1]+va) % H)
Li.append((Li[2]+va) % H)


M.append(np.array([0,0,c/2]))
M.append((va*2/3+vb*1/3+vc*5/6) % H)
M.append((va*1/3+vb*2/3+vc*1/6) % H)
M.append((M[0]+va) % H)
M.append((M[1]+va) % H)
M.append((M[2]+va) % H)

O_c = 0.2411
O1 = np.array([0,0,c*O_c])
O2 = np.array([0,0,c*(1-O_c)])
O.append(O1)
O.append(O2)
O.append((va*2/3+vb*1/3+vc*1/3+O1) % H)
O.append((va*1/3+vb*2/3+vc*2/3+O1) % H)
O.append((va*2/3+vb*1/3+vc*1/3+O2) % H)
O.append((va*1/3+vb*2/3+vc*2/3+O2) % H)
O.append((O[0]+va) % H)
O.append((O[1]+va) % H)
O.append((O[2]+va) % H)
O.append((O[3]+va) % H)
O.append((O[4]+va) % H)
O.append((O[5]+va) % H)

q = np.array([0,1,2,4,3,-2])
Y = np.array([0,1,3.344,4.0,2.04,-2.96])
q_Y = q - Y

Cell = [Li,M,O]


fp.write("\n\n%d\tatoms\t# core and shell atoms\n" % (2*24*Nx*Ny*Nz))
fp.write("%d\tbonds\t# core and shell atoms\n" % (24*Nx*Ny*Nz))

fp.write("\n10\tatom types\t# 5 cores and 5 shells\n")
fp.write("5\tbond types\n")

fp.write("\n0.0 %.5f xlo xhi\n" % Lx)
fp.write("0.0 %.5f ylo yhi\n" % Ly)
fp.write("-100.0 %.5f zlo zhi\n" % (Lz+100))

fp.write("\nMasses\t\t# core/shell mass ratio = %f\n\n" % ratio)
fp.write("1  %.5f\t# Li core\n" % (6.94*(1-ratio)))
fp.write("2  %.5f\t# Ni core\n" % (58.693*(1-ratio)))
fp.write("3  %.5f\t# Mn core\n" % (54.938*(1-ratio)))
fp.write("4  %.5f\t# Co core\n" % (58.933*(1-ratio)))
fp.write("5  %.5f\t# O  core\n" % (15.999*(1-ratio)))
fp.write("6  %.5f\t# Li shell\n" % (6.94*ratio))
fp.write("7  %.5f\t# Ni shell\n" % (58.693*ratio))
fp.write("8  %.5f\t# Mn shell\n" % (54.938*ratio))
fp.write("9  %.5f\t# Co shell\n" % (58.933*ratio))
fp.write("10 %.5f\t# O  shell\n" % (15.999*ratio))


fp.write("\nAtoms\n\n")

perm6 = np.matrix([[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1],[1,0,0,0,0,0]])

for i in range(Nx):
    for j in range(Ny):
        for k in range(Nz):
            R = vx*i+vy*j+vz*k   # Translation vector
            for XYZ in Li:
                xyz = (XYZ+R) % box
                Natom = Natom + 1
                Nmol = Nmol + 1
                fp.write("%6d %6d %6d  %6.3f  %16.10f %16.10f %16.10f\n" %(Natom,Nmol,1,q_Y[1],xyz[0],xyz[1],xyz[2]))
                Natom = Natom + 1
                fp.write("%6d %6d %6d  %6.3f  %16.10f %16.10f %16.10f\n" %(Natom,Nmol,6,Y[1],xyz[0],xyz[1],xyz[2]))
                bonds.append([Nmol,1,Natom-1,Natom])
            p0  = np.array([2,4,3,3,4,2])
            rand = np.random.randint(0,5)
            perm = np.linalg.matrix_power(perm6,rand)
            perm234 = (p0*perm).tolist()[0]
            M0 = 0
            for XYZ in M:
                xyz = (XYZ+R) % box
                Natom = Natom + 1
                Nmol = Nmol + 1
                M_type = perm234[M0]
                M0 = M0 + 1
                fp.write("%6d %6d %6d  %6.3f  %16.10f %16.10f %16.10f\n" %(Natom,Nmol,M_type,q_Y[M_type],xyz[0],xyz[1],xyz[2]))
                Natom = Natom + 1
                fp.write("%6d %6d %6d  %6.3f  %16.10f %16.10f %16.10f\n" %(Natom,Nmol,M_type+5,Y[M_type],xyz[0],xyz[1],xyz[2]))
                bonds.append([Nmol,M_type,Natom-1,Natom])
            for XYZ in O:
                xyz = (XYZ+R) % box
                Natom = Natom + 1
                Nmol = Nmol + 1
                fp.write("%6d %6d %6d  %6.3f  %16.10f %16.10f %16.10f\n" %(Natom,Nmol,5,q_Y[5],xyz[0],xyz[1],xyz[2]))
                Natom = Natom + 1
                fp.write("%6d %6d %6d  %6.3f  %16.10f %16.10f %16.10f\n" %(Natom,Nmol,10,Y[5],xyz[0],xyz[1],xyz[2]))
                bonds.append([Nmol,3,Natom-1,Natom])
fp.write("\nBonds\n\n")

for bd in bonds:
    fp.write("%6d %6d %6d %6d\n" % (bd[0],bd[1],bd[2],bd[3]))

fp.close()

