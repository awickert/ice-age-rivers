#! /usr/bin/python

# BEFORE RUNNING THIS, "toponow" NEEDS to be added (some other good too):
"""
r.proj loc=NA30as in=PminusETnow # Find codes from when I imported this!
r.proj loc=NA30as in=cellsize_km2
r.proj loc=NA30as in=cellsize_meters2
r.proj loc=NA30as in=zeros
r.proj loc=NA30as in=toponow # Not necessary in this version of code but maybe in future
v.proj loc=NA30as in=river_mouth_regions
r.proj loc=NA30as in=river_mouth_regions
"""

import os
import glob
import numpy as np
from Scientific.IO.NetCDF import NetCDFFile
from scipy.io import loadmat
import time
from grass.script import array as garray
from scipy.interpolate import interp1d
import multiprocessing

# Longitude boundaries
welon=190
eelon=320

def interpZonalAdvanced(topo,welon,eelon):
  """
  "topo" is the input file name/path
  "welon" is the western east-longitude
  "eelon" is the eastern east-longitude

  Interpolates along latitude bands to 1 arcminute to match the resolution 
  produced by JXM's group's spherical harmonic code in the meridional orientation

  Returns topoInterp, an interpolated file
  
  This function uses longitudinal limits and incorporates rollover at the edges 
  of the map
  """
  
  topo = np.loadtxt(topo)

  # copy first column onto right-hand end of topo array to allow interpolation
  # to wrap around
  # Will probably need a better solution in the end, with rollover, but this
  # should work for stuff away from the edges well enough
  col1 = topo[:,0]
  col1 = np.array([col1])
  col1 = col1.transpose()
  topo = np.hstack((topo,col1))

  # Possibly unnecessary except for the wraparound assist

  # east-longitudes: 513+1 columns from 0-360
  x = np.linspace(0,360,514)
  # Constrain to given range
  # Indices: [0] needed to get to 1D array
  westx_i = max((x<=welon).nonzero()[0])
  eastx_i = min((x>=eelon).nonzero()[0])
  # Give some padding for interpolation help
  westx_i = westx_i-3
  eastx_i = eastx_i+3
  # Calculate the latitudes that these respectively pertain to
  westx = x[westx_i]
  if eastx_i>=topo.shape[1]:
    eastx = x[eastx_i - t56.shape[1]]
  else:
    eastx = x[eastx_i]
  # Recalculate the x-value
  # +1 because both ends are counted
  xConstrained = np.linspace(westx,eastx,eastx_i-westx_i+1)
  
  # Start to make the constrained / wraparound-corrected topo file
  topoConstrained = np.array([topo[:,westx]]).transpose()
  
  # +1 west b/c we have already started to create this array.
  # +1 east so we don't short-change the east side
  for i in range(westx_i+1,eastx_i+1):
    # Rollover eastwards (westwards already good w/ negative indexing)
    if i>=topo.shape[1]:
      j = j -topo.shape[1]
    else:
      j = i
    # Concatenate onto the side until the array is built
    topoConstrained = np.hstack((topoConstrained, \
      np.array([topo[:,j]]).transpose()))

  # Interpolate

  # Create desired x-resolution (east-longitudes)
  if welon > eelon:
    eelon = eelon + 360 # Untested
  xnew = np.arange(welon+1/240.,eelon,1/120.)

  # All at once

  # Create a function f to perform 1D interpolation on the rows
  # Only interpolating through the desired output region, not to edges of
  # input data (those are there to help the cubic interp)
  f = interp1d(xConstrained,topoConstrained,kind='cubic')

  # Interpolate
  topoInterp = f(xnew) # <-- really quick! just a couple minutes for 70
                       # degrees of latitude (1 arcmin)!
                       # Even less when lon is constrained!
                      
  return topoInterp


def interpolate(item):
  print "***", item, "***"
  print "    Interpolating grid."
  topogrid_tmp = interpZonalAdvanced(item,welon,eelon)
  topogrid = np.zeros((8400,15600), dtype='float32')
  print "    N-S resizing for GEBCO."
  for i in range(8400):
    topogrid[i,:] = np.average(topogrid_tmp[i:i+2,:], axis=0)
  return topogrid


# Time step and file lists
files = np.array(sorted(glob.glob("topo*.txt")))
ages_numeric = np.loadtxt('ages')

# Start out by gridding at time = 0 (for corrections to modern topography)
#topo_spharm_now = interpolate(files[ages_numeric == 0])
print "Loading modern topography derived from GEBCO_08 (30 arcsecond) and"
print "interpolated etopo1 (60 arcsecond) bedrock surface data"
toponow_temp = garray.array()
toponow_temp.read('toponow')
toponow = toponow_temp.astype('int')
del toponow_temp
print "Generating array to correct spherical harmonic bumpiness"
correction = toponow - interpolate(files[ages_numeric == 0][0])

def worker(name, correction):
  starttime = time.time()
  outgrid = garray.array()
  print name, ': interpolating.'
  outgrid[...] = interpolate(item) + correction
  print "    Writing numpy array to GRASS GIS."
  outgrid.write(name)
  print 'Deleting outgrid object' # possibly unnecessary b/c in function, but playing safe.
  del outgrid
  print "    ", name, ':', time.time() - starttime, "seconds elapsed."

jobs = []
for item in files:
  name = os.path.splitext(files[0])[0]
  p = multiprocessing.Process(target=worker, args=(name, correction))
  jobs.append(p)
  p.start()


