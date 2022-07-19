import os
from ase import Atom, Atoms
from ase.build import bulk
from ase.cell import Cell
from ase.calculators.lammpsrun import LAMMPS
from ase.io import read, write
import numpy as np
import os

for i in range(1,10):
    filename = str(i)+'.CON'
    p1 = read(filename, format='vasp')
    write('cfg/'+str(i)+'.cfg', p1, format='cfg')
