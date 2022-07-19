#!/usr/bin python
eV2J = 1.60218e-19          # electron
mu   = 28                   # GPa
b    = 2.8                  # Angstrom
Rmax = 1000                 # Angstrom
N    = 1000

T    = mu*1e9*b*b*1e-20/eV2J*1e-10  # Line tension in eV/A

f    = open("string.bond","w")
f.write("#Bond potential for dislocation line tension\n")
f.write("\nLT\n");
f.write("N %d\n\n" %N)

for i in range(1,N+1):
    f.write("%d %f %.6g %.6g\n" %(i, i*Rmax/N, T*i*Rmax/N, -T))
f.close()
