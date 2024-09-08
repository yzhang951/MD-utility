export ASE_LAMMPSRUN_COMMAND=lmp_mc
ase_file="$(dirname $(python3 -c 'import ase; print(ase.__file__)'))/calculators/lammpsrun.py"
sed -i 's/line.startswith(_custom_thermo_mark)/line.strip\(\).startswith\("Step"\)/g' $ase_file

python3 -W "ignore" strength.py
