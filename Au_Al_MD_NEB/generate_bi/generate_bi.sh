export x='8'
export y='6'
export z='12'
xtal < fcc1.in
xtal < fcc2.in
mul 1.cfg $x $y $z 1.cfg
mul 2.cfg $x $y $z 2.cfg
vcut 1.cfg 0.5 0.5 0.5 1 0 0 1.cfg 1
vcut 2.cfg 0.5 0.5 0.5 -1 0 0 2.cfg 1

cfg2lammps 1.cfg 196.96654 Au
mv restart.in ../1.in
cfg2lammps 2.cfg 196.96654 Au
mv restart.in ../2.in
