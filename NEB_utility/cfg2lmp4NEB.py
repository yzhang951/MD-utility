import os, math
last_dir = '../0.000'
N = 64
log_fp = last_dir+'/log.lammps'
output = os.popen('grep -B 1 Climb '+log_fp+' | head -n1 | cut -d " " -f1').read()
cfg_num = output.split()[0]

def lmp2NEB(file_name):
    with open(file_name, 'r') as file:
        lines = file.readlines()

    with open(file_name, 'w') as file:
        in_atoms_section = False
        for line in lines:
            if "Atoms " in line:
                in_atoms_section = True
                continue
            if in_atoms_section:
                parts = line.split()
                if(len(parts)<1):
                    continue
                atom_id = int(parts[0])
               
                parts[1] = ''
                file.write("\t" + "\t\t".join(parts) + "\n")
            else: 
                if " atoms" in line:     
                    parts = line.split()
                    file.write(parts[0] + "\n")


os.system('rm *.lmp')
for i in range(N):
    command = 'atomsk '+last_dir+'/Cu'+cfg_num+'_'+str(i)+'.cfg ' + str(i) + '.lmp'
    os.system(command)
    if i==0:
        os.system('cp 0.lmp init.lmp')
    lmp2NEB(str(i)+'.lmp')
