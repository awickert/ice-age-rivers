#! /usr/bin/python

import os
import glob
import numpy as np
from Scientific.IO.NetCDF import NetCDFFile
from scipy.io import loadmat
import time
import grass.script as grass
from grass.script import array as garray
from scipy.interpolate import interp1d
import multiprocessing
import re

# Add important grids
inloc='IceGlobal30as_template'
grass.run_command('g.region', n=90, s=-90, w=-180, e=180, res='0:0:30')
try:
  grass.run_command('r.proj', loc=inloc, input='cellsize_meters2')
except:
  pass
try:
  grass.run_command('r.proj', loc=inloc, input='toponow')
except:
  pass
try:
  grass.run_command('v.proj', loc=inloc, input='river_mouth_regions')#, overwrite=True)
except:
  pass
try:
  grass.run_command('r.proj', loc=inloc, input='river_mouth_regions')#, overwrite=True)
except:
  pass
  
# Longitude boundaries
welon=0
eelon=360

def interpZonal30(topo,welon,eelon):
  """
  "topo" is the input file name/path
  "welon" is the western east-longitude
  "eelon" is the eastern east-longitude

  Interpolates along latitude bands to 30-arcsecond to match the resolution 
  produced by JXM's group's spherical harmonic code in the meridional orientation

  Returns topoInterp, an interpolated file
  
  This function uses longitudinal limits and incorporates rollover at the edges 
  of the map
  """
  
  import numpy as np
  from scipy.interpolate import interp1d

  topo = np.loadtxt(topo)

  # copy first column onto right-hand end of topo array to allow interpolation
  # to wrap around
  # Will probably need a better solution in the end, with rollover, but this
  # should work for stuff away from the edges well enough
  col1 = topo[:,0]
  col1=np.array([col1])
  col1 = col1.transpose()
  topo = np.hstack((topo,col1))

  # Possibly unnecessary except for the wraparound assist

  # east-longitudes: 513+1 columns from 0-360
  x = np.linspace(0,360,514)
  # Constrain to given range

  # Interpolate

  # Create desired x-resolution (east-longitudes)
  if welon > eelon:
    eelon = eelon + 360 # Untested
  xnew = np.arange(welon+1/240.,eelon,1/120.)

  # All at once

  # Create a function f to perform 1D interpolation on the rows
  # Only interpolating through the desired output region, not to edges of
  # input data (those are there to help the cubic interp)
  
  # NOTE -- I REMOVED THE PARTS THAT LET YOU SUBSET THE MAP!
  f = interp1d(x,topo,kind='cubic')

  # Interpolate
  topoInterp = f(xnew) # <-- really quick! just a couple minutes for 70
                       # degrees of latitude (1 arcmin)!
                       # Even less when lon is constrained!
                      
  return topoInterp
  

def interpolate(item):
  print "***", item, "***"
  print "    Interpolating grid for", item
  topogrid_tmp = interpZonal30(item,welon,eelon)
  topogrid = np.zeros((21600,43200), dtype='float32')
  print "    N-S resizing for GEBCO for", item
  for i in range(topogrid.shape[0]):
    topogrid[i,:] = np.average(topogrid_tmp[i:i+2,:], axis=0)
  print "    DEM genereation for", item, "complete."
  return topogrid


# Time step and file lists
files = np.array(sorted(glob.glob("topo*.txt")))

ages = []
for fname in files:
  ages += re.findall('\d+', fname)
ages = np.array(ages)
ages_numeric = ages.astype(int)

# Climate model run (TraCE-21K) goes back only to 22 ka
ages = ages[(ages_numeric > 22000) + (ages_numeric == 0)]
ages_numeric = ages_numeric[(ages_numeric > 22000) + (ages_numeric == 0)]


#############################
# STARTS TO USE MEMORY HERE #
#############################

###########################
# TOPOGRAPHY / BATHYMETRY #
###########################

# These maps run from 0 to 360
grass.run_command('g.region', n=90, s=-90, w=0, e=360, res='0:0:30')

# Start out by gridding at time = 0 (for corrections to modern topography)
print "Loading modern topography derived from GEBCO_08 (30 arcsecond) and"
print "  interpolated etopo1 (60 arcsecond) bedrock surface data"
toponow_temp = garray.array()
toponow_temp.read('toponow')
toponow = toponow_temp.astype('int')
del toponow_temp
print "Generating array to correct spherical harmonic bumpiness"
try:
  # Write "correction" to location in future versions?
  # Yes -- starting with G12
  correction = toponow - interpolate(files[ages_numeric == 0][0])
except:
  print "Projecting to present -- hard-coded for G12"
  t500 = interpolate(files[0])
  #t600 = interpolate(files[1])
  #diff = t500 - t600
  diff = (toponow - t500)/5.
  #del t600
  # NEW METHOD -- NO DIFFERENCE EXCEPT FOR FIXING MISTAKE (HOW?) FOR 
  # TOPO_000100
  """
  correction = toponow - (t500 + 5*diff)
  outarray = garray.array()
  outarray[...] = correction
  outarray.write('correction')
  outarray[...] = t500 + diff + correction
  outarray.write('topo_000400', overwrite=True)
  outarray[...] = t500 + 2*diff + correction
  outarray.write('topo_000300', overwrite=True)
  outarray[...] = t500 + 3*diff + correction
  outarray.write('topo_000200', overwrite=True)
  outarray[...] = t500 + 4*diff + correction
  outarray.write('topo_000100', overwrite=True)
  outarray[...] = toponow
  outarray.write('topo_000000', overwrite=True)
  """
  outarray = garray.array()
  outarray[...] = t500 + diff
  outarray.write('topo_000400', overwrite=False)
  outarray[...] = t500 + 2*diff
  outarray.write('topo_000300', overwrite=False)
  outarray[...] = t500 + 3*diff
  outarray.write('topo_000200', overwrite=False)
  outarray[...] = t500 + 4*diff
  outarray.write('topo_000100', overwrite=False)
  outarray[...] = toponow
  outarray.write('topo_000000', overwrite=True)
  del t500
  del diff
  del outarray

  for topomap in 'topo_000400', 'topo_000300', 'topo_000200', 'topo_000100':
    grass.run_command('r.colors', map=topomap, color='etopo2')

def worker(name, correction):
  starttime = time.time()
  outgrid = garray.array()
  print name, ': interpolating.'
  outgrid[...] = interpolate(item) + correction
  print "    Writing numpy array to GRASS GIS."
  outgrid.write(name, overwrite=True)
  print 'Deleting outgrid object' # possibly unnecessary b/c in function, but playing safe.
  del outgrid
  print "    ", name, ':', time.time() - starttime, "seconds elapsed."

# Parallel: but limited by memory
# 8x2 Gb per process, so 3 is reasonable. 4, maybe pushing it but speeds things up.
n_file_start = 0
n_simultaneous_jobs = 4
p = None # Nothing until defined by a job
for i in range( int(np.ceil( len(files) / float(n_simultaneous_jobs) )) ):
  files_subset = files[n_file_start:n_file_start+n_simultaneous_jobs]
  jobs = []
  for item in files_subset:
    age = re.findall('\d+', item)[0]
    name = 'topo_' + age
    # Do not overwrite preexisting rasters
    # Good if this memory-intensive process crashes
    if grass.parse_command('g.list', type='raster', pattern=name):
      print name, "already exists. Skipping."
    else:
      p = multiprocessing.Process(target=worker, args=(name, correction))
      jobs.append(p)
      p.start()
  if p:
    while p.is_alive():
      pass
  n_file_start += n_simultaneous_jobs

# Give a better colormap
for topomap in sorted(grass.parse_command('g.list', type='raster', pattern='topo_*')):
  grass.run_command('r.colors', map=topomap, color='etopo2')

# COMMENT THIS LATER -- FOR EXTRA IMPORT STEP WHILE ADDING TIME-STEPS
# M180_180 REQUIRED FOR PROJECTION TO ANTARCTICA
for age in ages[1:]:
  grass.mapcalc('topo_m180_180_'+age+' = topo_'+age)

for age in ages:
  grass.run_command('r.proj', location='ICE6G_global', input='ice_raw_m180_180_'+age, out='ice_raw_'+age, overwrite=True)

#for age in ['023000', '024000', '025000', '026000']:
#  grass.run_command('r.proj', location='globalNoAnt', input='topo_m180_180_'+age, out='topo_'+age, overwrite=True)

for age in ages:
  grass.run_command('r.proj', location='ICE6G_global', input='topo_m180_180_'+age, out='topo_'+age, overwrite=True)

for age in ages:
  grass.run_command('r.proj', location='ICE6G_global', input='ice_'+age, out='topo_'+age, overwrite=True)

for age in ages:
  grass.run_command('r.proj', location='ICE6G_global', input='ice_'+age, out='topo_'+age, overwrite=True)


"""
#band-aid
for age in ages[10:12]:
  grass.run_command('g.region', rast='topo_000000')
  print age
  grass.run_command('r.proj', location='G12_global', input='topo_'+age, overwrite=True)
"""

######################
# WATER BALANCE MAPS #
######################

# Set region -- set up for WB (water balance)
#grass.run_command('g.region', n=nlat, s=slat, w=wlon, e=elon, res=0.25)
# Need these boundaries to project something that isn't NULL! (-180, 180)
grass.run_command('g.region', n=90, s=-90, w=-180, e=180, res=0.25)

# Intelligently import these so they correspond to the ages of the time-steps
# Use weighted averaging.
#wb_ages = np.arange(22000, -1, -100)

for age in ages:
  print '***', age, '***'
  try:
    wb = 'wb_'+age
    grass.run_command('r.proj', loc=inloc, input=wb, quiet=True, overwrite=True)
  except:
    print "Linearly interpolating water balance approx. for", age, "years ago"
    lowage = '%06d' %int(np.floor(int(age)/100.)*100)
    highage = '%06d' %int(np.ceil(int(age)/100.)*100)
    grass.run_command('r.proj', location=inloc, input='wb_'+lowage, output='_wb_'+lowage, quiet=True, overwrite=True)
    grass.run_command('r.proj', location=inloc, input='wb_'+highage, output='_wb_'+highage, quiet=True, overwrite=True)
    weight_high = ( float(age) - float(lowage) ) / ( float(highage) - float(lowage) )
    weight_low = 1 - weight_high
    mcstr = 'wb_'+age+' = ('+str(weight_high)+' * _wb_'+highage+') + ('+str(weight_low)+' * _wb_'+lowage+')'
    grass.mapcalc(mcstr, overwrite=True)
    grass.run_command('g.remove', type='raster', pattern='_wb_*', flags='f')
    
