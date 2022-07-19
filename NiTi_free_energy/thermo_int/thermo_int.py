import numpy as np
from scipy.interpolate import interp1d
import os
import matplotlib.pyplot as plt


H = []
sig = []
for i in range(11):
    log = 'log.lammps.'+str(i) 
    os.system('grep Step -A 301 '+log+'>tmp')

    data = np.loadtxt('tmp', comments='Step')

    H6 = data[0,3:9]
    H.append(np.matrix([[H6[0],0,0],[H6[3],H6[1],0],[H6[4],H6[5],H6[2]]], dtype=np.double))
    s = data[10:,9:]
    s = np.mean(s,axis=0) 
    sig.append(np.matrix([[s[0],s[3],s[4]], [s[3],s[1],s[5]], [s[4],s[5],s[2]]], dtype=np.double)/-1E4) # Stress in GPa

F = []
F.append(0.0)
for i in range(1,11):
    dH = H[i] - H[i-1]
    det0 = np.linalg.det(H[i-1])
    det1 = np.linalg.det(H[i])
    func0 = det0*np.matmul(sig[i-1], np.linalg.inv(H[i-1]))
    func1 = det1*np.matmul(sig[i], np.linalg.inv(H[i]))
    dF = np.tensordot((func0+func1)/2.0, dH, axes=2)
    F.append(F[i-1]+dF)
print(F)
neb = np.loadtxt('fe.out', comments='Image')
y = neb[:,2]*1000/4

#x = np.linspace(0, 1, num=11, endpoint=True)
x = neb[:,1]/neb[10,1]
F = np.array(F)/160/864*1000

f = interp1d(x, F)
f2 = interp1d(x, F, kind='cubic')

f3 = interp1d(x, y, kind='cubic')

xnew = np.linspace(0, 1, num=41, endpoint=True)
plt.plot(x, F, 'bo', xnew, f2(xnew), 'b-', x, y, 'ro', xnew, f3(xnew),'r-')
plt.legend([None,'FE',None,'PE'])
plt.xlabel('Reaction coordinate')
plt.ylabel('Energy (meV/atom)')
plt.show()
