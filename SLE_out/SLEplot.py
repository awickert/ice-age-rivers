#! /usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt

ICE3G = np.genfromtxt('ICE3G_NA_SLE.txt')
ANU = np.genfromtxt('ANU_NA_SLE.txt')
ICE5G = np.genfromtxt('ICE5G_NA_SLE.txt')
ICE6G = np.genfromtxt('ICE6G_NA_SLE.txt')
G12 = np.genfromtxt('G12_NA_SLE.txt')
L2014global = np.genfromtxt('Lambeck_SL_curve.txt', skip_header=8)

def aa(ts):
  return (ts[:-1] + ts[1:]) / 2.

fig = plt.figure(figsize=(6,9))
ax1 = fig.add_subplot(211)

modern_mean_ice = ICE6G[-1,-1]

#plt.fill_between(L2014global[:,0], -L2014global[:,2]+L2014global[:,3]+modern_mean_ice, -L2014global[:,2]-L2014global[:,3]+modern_mean_ice, color='k', linewidth=.01, alpha=.6, label='L2014 global + mNA')
plt.plot(L2014global[:,0], -L2014global[:,2]+modern_mean_ice, color='k', linewidth=4, alpha=.6, label='L2014 Global + MNA')
plt.plot(ICE3G[:,0]/1000., ICE3G[:,2], color='green', linewidth=4, alpha=.6, label='ICE-3G')
plt.plot(ANU[:,0]/1000., ANU[:,2], color='purple', linewidth=4, alpha=.6, label='ANU')
plt.plot(ICE5G[:,0]/1000., ICE5G[:,2], color='blue', linewidth=4, alpha=.6, label='ICE-5G')
plt.plot(G12[:,0]/1000., G12[:,2], color='orange', linewidth=4, alpha=.6, label='G12')
plt.plot(ICE6G[:,0]/1000., ICE6G[:,2], color='red', linewidth=4, alpha=.6, label='ICE-6G')
plt.xlim(0, 20)
plt.legend(loc='upper left', fontsize=16)

plt.grid(b=True, which='both', axis='x', color='0.65',linestyle='-')

plt.ylabel('Ice sheet-volume [m SLE]', fontsize=16)

ax1 = fig.add_subplot(212)
plt.plot(aa(L2014global[:,0]), np.diff(-L2014global[:,2]*3.94765540E14)/np.diff(L2014global[:,0]*1000)/3.15E7/1E6, color='k', linewidth=4, alpha=.6)
plt.plot(aa(ICE3G[:,0]/1000.), np.diff(ICE3G[:,1])/np.diff(ICE3G[:,0])/3.15E7/1E6, color='green', linewidth=4, alpha=.6)
plt.plot(aa(ANU[:,0]/1000.), np.diff(ANU[:,1])/np.diff(ANU[:,0])/3.15E7/1E6, color='purple', linewidth=4, alpha=.6)
plt.plot(aa(ICE5G[:,0]/1000.), np.diff(ICE5G[:,1])/np.diff(ICE5G[:,0])/3.15E7/1E6, color='blue', linewidth=4, alpha=.6)
plt.plot(aa(G12[:,0]/1000.), np.diff(G12[:,1])/np.diff(G12[:,0])/3.15E7/1E6, color='orange', linewidth=4, alpha=.6)
plt.plot(aa(ICE6G[:,0]/1000.), np.diff(ICE6G[:,1])/np.diff(ICE6G[:,0])/3.15E7/1E6, color='red', linewidth=4, alpha=.6)
plt.xlim(0, 20)
plt.ylim(-0.1, 0.5)

plt.ylabel('Total meltwater discharge [Sv]', fontsize=16) #m$^3$ s$^-1$
plt.xlabel('Age [ka]', fontsize=16)

plt.grid(b=True, which='both', axis='x', color='0.65',linestyle='-')

plt.tight_layout()

plt.savefig('/home/awickert/Dropbox/Papers/InProgress/DrainageMethods/figures/inprogress/Vice_Qice_with_global.svg')
#plt.savefig('/home/awickert/Dropbox/Papers/InProgress/DrainageMethods/figures/Vice_Qice_with_global.pdf')

plt.show()

