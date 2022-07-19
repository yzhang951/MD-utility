import os
from ase import Atom, Atoms
from ase.build import bulk
from ase.cell import Cell
from ase.calculators.lammpsrun import LAMMPS
from ase.io import read, write
import numpy as np
import os

for i in range(11):
    if i<10:
        filename = 'linear00'+str(i)+'.cfg'
        p1 = read(filename)
    else:
        filename = 'linear0'+str(i)+'.cfg'
        p1 = read(filename)

    write('data.in', p1, format='lammps-data')
    os.system('mpiexec -np 8 lmp_deepmd -in lmp.in')
    os.system('cp log.lammps log.lammps.'+str(i))
