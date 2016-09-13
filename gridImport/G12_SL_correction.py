import os
import glob
import numpy as np
from Scientific.IO.NetCDF import NetCDFFile
from scipy.io import loadmat
from scipy.interpolate import interp1d
import time
import grass.script as grass
from grass.script import array as garray
from scipy.interpolate import interp1d
import multiprocessing
import re
import fnmatch

age_SLcorr = np.genfromtxt('rest_of_world_SL_ICE5G_for_G12.txt', delimiter='|')
age_SLcorr = age_SLcorr[::-1]
age = age_SLcorr[:,0]
SLcorr = age_SLcorr[:,-1]

f = interp1d(age, SLcorr)

topomaps = sorted( grass.parse_command('g.list', type='raster', pattern='topo_??????').keys() )[::-1]
topomaps = np.array(topomaps)
# on time-steps
ages = []
ages_numeric = []
for topomap in topomaps:
  # on time-step
  ages.append(re.findall('\d+', topomap)[0])
  # Use float instead of int to avoid undesired
  # floor division issues
  ages_numeric.append(float(ages[-1]))
ages = np.array(ages)
ages_numeric = np.array(ages_numeric)

SLcorr_ts = f(ages_numeric)

def worker(name, correction):
  starttime = time.time()
  outgrid = garray.array()
  print name, ': reading.'
  outgrid.read(name)
  print name, ': correcting.'
  outgrid[...] += correction
  print "    Writing numpy array to GRASS GIS."
  outgrid.write(name, overwrite=True)
  print '    Deleting outgrid object' # possibly unnecessary b/c in function, but playing safe.
  del outgrid
  print "    ", name, ':', time.time() - starttime, "seconds elapsed."

# Parallel: but limited by memory
# 8 Gb per process, so 6 is reasonable. 8, maybe pushing it but speeds things up.
# 8 also starts to run towards the limit of processors (I have 4 other operations going on and 12 cores)
n_file_start = 0
n_simultaneous_jobs = 8
p = None # Nothing until defined by a job
for i in range( int(np.ceil( len(SLcorr_ts) / float(n_simultaneous_jobs) )) ):
  topomaps_subset = topomaps[n_file_start:n_file_start+n_simultaneous_jobs]
  jobs = []
  for topomap in topomaps_subset:
    correction = float(SLcorr_ts[topomaps == topomap])
    print topomap, correction
    p = multiprocessing.Process(target=worker, args=(topomap, correction))
    jobs.append(p)
    p.start()
  if p:
    while p.is_alive():
      pass
  n_file_start += n_simultaneous_jobs



# Check which ones were messed up
# g.region -p n=12:00:30n s=12n w=58:10:30w e=58:10w
a = garray.array()
z = []
for topomap in topomaps:
  print topomap
  a.read(topomap)
  z.append(float(a))

plt.plot(ages_numeric, z); plt.show()


a = garray.array()
z = []
for age in ages:
  print age
  a.read('topo_'+age)
  z.append(float(a))




# Old way -- better.
# Have SLcorr defined backwards -- SL fall + (topo rise)
grass.run_command('g.region', rast='topo_000000', flags='p')
for i in range(len(ages)-1):
  print 'topo_'+ages[i], ages_numeric[i], SLcorr_ts[i]
  grass.run_command('g.rename', rast='topo_'+ages[i]+',tmp', overwrite=True)
  mcstr = 'topo_'+ages[i]+' = tmp + 2*'+str(SLcorr_ts[i]) # 2* to fix previous backwards attempt
  print mcstr
  grass.mapcalc(mcstr)
  
