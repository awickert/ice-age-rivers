#! /usr/bin/env python

# Goes in ~/models/jxm_model/SLcurves_River_mouths_5G

from glob import glob # to get a list of files of interest
import numpy as np
from matplotlib import pyplot as plt
import os

# Import files
topofiles5G = glob('/home/awickert/Dropbox/jxm_model/SLcurves_River_mouths_5G/topo*.txt') # Elevation w.r.t. sea level at that time
topofilesG12 = glob('/home/awickert/Dropbox/jxm_model/SLcurves_River_mouths_Gregoire/topo*.txt')
topofiles3G = glob('/home/awickert/Dropbox/jxm_model/SLcurves_River_mouths_3G/topo*.txt')
topofilesANU = glob('/home/awickert/Dropbox/jxm_model/SLcurves_River_mouths_ANU/topo*.txt')
topofiles6G = glob('/home/awickert/Dropbox/jxm_model/SLcurves_River_mouths_6G/topo*.txt')

lat_lon_name = np.loadtxt('/home/awickert/Dropbox/jxm_model/SLcurves_River_mouths_5G/site_list_lat_lon_name', delimiter=', ', dtype=str)
lat = lat_lon_name[:,0].astype(float)
lon = lat_lon_name[:,1].astype(float)
rivernames = lat_lon_name[:,2]

# zero-pad for ICE-5G
for i in range(len(topofiles5G)):
  topofile = topofiles5G[i]
  if len(topofile) == 9:
    print topofiles6G[i]
    topofiles5G[i] = topofile[:4] + '0' + topofile[-5:]
    os.rename(topofile, topofiles5G[i])

# and ICE-6G
for i in range(len(topofiles6G)):
  topofile = topofiles6G[i]
  if len(os.path.split(topofile)[1]) == 9:
    print topofiles6G[i]
    topofiles6G[i] = topofile[:-5] + '0' + topofile[-5:]
    os.rename(topofile, topofiles6G[i])

topofiles5G = sorted(topofiles5G) # need to get them into order
topofiles6G = sorted(topofiles6G) # need to get them into order
topofilesG12 = sorted(topofilesG12) # need to get them into order
topofilesANU = sorted(topofilesANU)[:-1] # need to get them into order, ends after present
topofiles3G = sorted(topofiles 3G)[:-1] # need to get them into order, ends after present (not in time stamps)

time5G = np.loadtxt('/home/awickert/Dropbox/jxm_model/SLcurves_River_mouths_5G/step_time')[:,-1]
timeG12 = np.linspace(21, 0, 211)[:-1] # Also have a file for this
timeANU = np.loadtxt('/home/awickert/Dropbox/jxm_model/SLcurves_River_mouths_ANU/ages_ANU')[:-1] # ends after present
time3G = np.loadtxt('/home/awickert/Dropbox/jxm_model/SLcurves_River_mouths_3G/ages_3G')
time6G = np.hstack(( np.array([25, 24, 23, 22, 21]), np.arange(20.5, -.1, -.5) ))

def SLcurves(topofiles, time):

  # Correct for modern topography (on same spherical harmonic expansion)
  topo_raw_now = np.loadtxt(topofiles[-1])
  if topo_raw_now.ndim == 1:
    topo_raw_now = np.array([topo_raw_now]) # add dimension
  inner_correction = -0-topo_raw_now[:,2]
  print inner_correction

  # Elevation matrix
  elev = np.zeros(( len(time), 2, topo_raw_now.shape[0] ))
  for i in range(topo_raw_now.shape[0]):
    elev[:,0,i] = time

  for i in range(len(topofiles)):
    elev[i,1,:] = np.loadtxt(topofiles[i])[:,-1] + inner_correction

  return elev

elev5G = SLcurves(topofiles5G, time5G)
elevG12 = SLcurves(topofilesG12, timeG12)
elevANU = SLcurves(topofilesANU, timeANU)
elev3G = SLcurves(topofiles3G, time3G)
elev6G = SLcurves(topofiles6G, time6G)

timeG12 = np.concatenate((timeG12, [0]))
elevG12 = np.concatenate((elevG12, [2*[9*[0]]]), axis=0)

time_LGM5G = list(time5G[55:])
time_LGM5G.pop(-2)

G12_global_extra = np.genfromtxt('/home/awickert/Dropbox/jxm_model/SLcurves_River_mouths_Gregoire/SL_corrections_21_0.5_0_ka.txt')
from scipy.interpolate import interp1d
f = interp1d(time_LGM5G[::-1], G12_global_extra)
G12_global_extra_all_times = f(timeG12[::-1])

for i in range(len(rivernames)):
  elevG12[:, 1, i] -= G12_global_extra_all_times

# Save
# np.savetxt('sea_level_curves.txt', elev, delimiter=' ')

# Plot
# LGM to present
# Y TICK LABEL CHANGES!!!
fig, axes = plt.subplots(elev5G.shape[-1], sharex=True, sharey=False, figsize=(6, 1.8*elev5G.shape[-1]))
for i in range(len(axes)):
  ax = axes[i]
  ax.set_title(rivernames[i], fontsize=12, fontweight='bold') # horizontalalignment='left', 
  if i == elev5G.shape[-1]/2: # use floor division if necessary, "+1" handled by 0-indexing
    ax.set_ylabel('Local relative sea level [m]', fontsize=20, fontweight='bold')
  #ax.plot(time5G, esl, 'k-')
  ax.plot(time3G, -elev3G[:,1,i], color='green', linewidth=4, alpha=.5)
  ax.plot(time5G, -elev5G[:,1,i], 'b-', linewidth=4, alpha=.5)
  ax.plot(time6G, -elev6G[:,1,i], color='red', linewidth=4, alpha=.5)
  ax.plot(timeANU, -elevANU[:,1,i], color='purple', linewidth=4, alpha=.5)
  ax.plot(timeG12, -elevG12[:,1,i], color='orange', linewidth=4, alpha=.5)
  ax.set_xlim(0,21)
  ax.locator_params(axis='y', nbins=6)
plt.xlabel('Age [ka]', fontsize=20, fontweight='bold')
fig.tight_layout()
fig.subplots_adjust(hspace=.5, bottom=.07)
plt.savefig('/home/awickert/Dropbox/Papers/InProgress/DrainageMethods/figures/inprogress/sea_level_at_mouths.svg')
plt.show()

# Location map
plt.figure(2)
m = Basemap(width=7000000,height=6500000,
            resolution='l',projection='laea',\
            lat_ts=52,lat_0=52,lon_0=-100.)
m.drawcoastlines()
m.drawcountries()
m.drawstates()
x, y = m(lon, lat)
#xoffset = 
m.plot(x, y, 'ko', markersize=12, marker='d')
for label, xpt, ypt in zip(rivername, x, y):
  plt.text(xpt+100000, ypt+50000, label)

plt.show()

