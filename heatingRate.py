# heatingRate.py
#
# create a line plot that shows possible heating rates and time scales
#
# LKS, November 2015
# constants
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.collections as collections

#
# heating rate eq is dW/dt =(1/2) (q^2 /m) Psi
q = 1.6*1e-19
m = 1.67*1e-27
t1=1000
t2=5000
t3=3600
Psi=np.logspace(-13,-10, 31) # logspaced bins


DW1=t1*(1/2.)*(q**2)*np.array(Psi)/(m*q)
DW2=t2*(1/2.)*(q**2)*np.array(Psi)/(m*q)
DW3=t3*(1/2.)*(q**2)*np.array(Psi)/(m*q)


#
# now plot this 
from numpy import ma
cmap=plt.cm.jet
# intitialize the figure
fig=plt.figure()
ax=fig.add_subplot(111)
ax.set_xlabel('$\Psi$ [V$^{2}$/m$^{2}$/Hz]', fontsize=30, fontweight='bold')
ax.set_ylabel('dW [eV]', fontsize=30, fontweight='bold')
plt.subplots_adjust(right=0.8, top=0.9, bottom=0.15, left=0.2)
ax.set_yscale('log')
ax.set_xscale('log')
#ax.set_ylim(10000, 30000)
#ax.set_xlim(2, 3)
ax.plot(Psi, DW1, lw=3, c='k')
ax.plot(Psi, DW2, lw=3, c='k')
ax.plot(Psi, DW3, lw=3, c='k', linestyle='--')
plt.fill_between(Psi, DW1, DW2, color='grey', alpha='0.5')
ax.tick_params(axis='both', which='major', labelsize=30)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
subdir_name='B_Freq_Map'
plt.draw()
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name, 0777) 
os.chdir(subdir_name)#
fig.set_size_inches(13,9)
plt.savefig('heatingRate.pdf')
plt.close(fig)
os.chdir('..')
