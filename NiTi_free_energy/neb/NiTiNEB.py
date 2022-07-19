import os
from ase import Atom, Atoms
from ase.build import bulk
from ase.cell import Cell
#from ase.calculators.lammpsrun import LAMMPS
from tsase import neb
from tsase.calculators.lammps_ext import LAMMPS
from ase.io import read
import numpy as np
import os
os.system('rm trash/*')

parameters = {'pair_style': 'deepmd 000.pb',
              'pair_coeff': ['* *'],
              'atom_style':'charge',
              'mass':['1 58.69','2 47.867']}

files = ['000.pb']

a = 4.24199
b = 2.99954
c = 4.24199

h = np.array([[a,0,0,],[0,b,0],[0,0,c]])

cell110 = Cell(h)
Ni1 = Atoms('Ni', positions=[cell110.cartesian_positions([0.75,0,0.25])], cell=cell110, pbc=[1,1,1])
Ni2 = Atoms('Ni', positions=[cell110.cartesian_positions([0.25,0,0.75])], cell=cell110, pbc=[1,1,1])
Ti3 = Atoms('Ti', positions=[cell110.cartesian_positions([0.25,0.5,0.25])], cell=cell110, pbc=[1,1,1])
Ti4 = Atoms('Ti', positions=[cell110.cartesian_positions([0.75,0.5,0.75])], cell=cell110, pbc=[1,1,1])

p1 = Ni1 + Ni2 + Ti3 + Ti4


h = np.array([[4.63645,0,0,],[-0.312139,2.76295,0],[0,0,4.196]])

cell110 = Cell(h)
Ni1 = Atoms('Ni', positions=[cell110.cartesian_positions([0.675096,0.0217455,0.25])], cell=cell110, pbc=[1,1,1])
Ni2 = Atoms('Ni', positions=[cell110.cartesian_positions([0.324904,0.978254,0.75])], cell=cell110, pbc=[1,1,1])
Ti3 = Atoms('Ti', positions=[cell110.cartesian_positions([0.209995,0.421925,0.25])], cell=cell110, pbc=[1,1,1])
Ti4 = Atoms('Ti', positions=[cell110.cartesian_positions([0.790005,0.578075,0.75])], cell=cell110, pbc=[1,1,1])

p2 = Ni1 + Ni2 + Ti3 + Ti4

for p in p1:
    p.charge = 0.0
for p in p2:
    p.charge = 0.0 

calc = LAMMPS(parameters=parameters, files=files, tmp_dir='trash')

p1.set_calculator(calc)
p2.set_calculator(calc)

# external stress applied in the unit of GPa
stress=np.zeros((3,3))
#stress[2,2] = -1.0  #negative is tension

# initialize the band
nim = 16    # number of images, including end points
# no climbing image first
band = neb.ssneb(p1, p2, numImages = nim, express = stress)
#band = neb.ssneb(p1, p2, numImages = nim, method = 'ci', express = stress)
'''
#-------------- substitute for the above line if parallelize over image ------------------
band = neb.ssneb(p1, p2, numImages = nim, method = 'ci', express = stress, parallel = True)
'''

# to restart, uncomment the following lines which read the previous optimized images into the band
#for i in range(1,nim-1):
#    filename = str(i)+'.CON'
#    b = read(filename,format='vasp')
#    band.path[i].set_positions(b.get_positions())
#    band.path[i].set_cell(b.get_cell())


#opt = neb.qm_ssneb(band, maxmove = 0.10, dt = 0.05)
opt = neb.fire_ssneb(band, maxmove =0.1, dtmax = 0.1, dt=0.1)
opt.minimize(forceConverged=0.01, maxIterations = 300)


