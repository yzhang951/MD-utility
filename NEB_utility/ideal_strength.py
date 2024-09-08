import os
from ase import Atom, Atoms
from ase.build import bulk
from ase.cell import Cell
from ase.calculators.lammpsrun import LAMMPS
from ase.io import read
from ase.io.cfg import write_cfg
from ase.units import create_units
units = create_units('2014')
import numpy as np
import os
os.system('rm tmp/*')
os.system('export ASE_LAMMPSRUN_COMMAND=lmp_mc')

parameters = {'pair_style': 'eam/alloy',
              'pair_coeff': ['* * ../Cu_mishin1.eam.alloy Cu'],
              'atom_style':'atomic',
              'mass':['1 63.54'],
              }

files = ['../Cu_mishin1.eam.alloy']

a = 3.615000000000000

h = np.array([[a,0,0,],[0,a,0],[0,0,a]])
cell = Cell(h)

Cu1 = Atoms('Cu', positions=[cell.cartesian_positions([0,0,0])], cell=cell, pbc=[1,1,1])
Cu2 = Atoms('Cu', positions=[cell.cartesian_positions([0.5,0.5,0])], cell=cell, pbc=[1,1,1])
Cu3 = Atoms('Cu', positions=[cell.cartesian_positions([0.5,0,0.5])], cell=cell, pbc=[1,1,1])
Cu4 = Atoms('Cu', positions=[cell.cartesian_positions([0,0.5,0.5])], cell=cell, pbc=[1,1,1])


p  = Cu1 + Cu2 + Cu3 + Cu4

calc = LAMMPS(parameters=parameters, files=files, tmp_dir='tmp')
p.set_calculator(calc)

energy = p.get_potential_energy()
e0 = energy

# Deformation gradient, F=I+eps*(m*n), here m*n is the dyadic product of slip plane and slip direction
eps_max = 1.0
max_inc = 100
# m: slip plane normal, n: slip direction
m = np.array([ 1, 0, 0])
n = np.array([ 0, 2, 0])

for i in range(max_inc):
    eps = eps_max*(i/max_inc - 0.0)
    F = np.array([[1,0,0,],[0,1,0],[0,0,1]])+eps*np.tensordot(m,n,axes=0)
    h1 = np.matmul(F, h)
    #print(h1/a)
    p.set_cell(h1, scale_atoms=True)
    p.label=str(i)
    energy = p.get_potential_energy()
    s = p.get_stress() # Unit?
    stress = np.array([[s[0],s[5],s[4]],[s[5],s[1],s[3]],[s[4],s[3],s[2]]])/units['GPa']
    #print(stress)
    tmp = np.matmul(stress, n/np.linalg.norm(n))
    tau = np.matmul(m/np.linalg.norm(m), tmp)
    write_cfg('tmp/'+str(i)+'.cfg', p)
    print(eps, tau, energy-e0)
