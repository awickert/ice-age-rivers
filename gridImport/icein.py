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

import Fcn30as as fcn
import os
import glob
import numpy as np
from Scientific.IO.NetCDF import NetCDFFile
from scipy.io import loadmat
import time
from grass.script import array as garray

# Longitude boundaries
welon=190
eelon=320

def interpolate(item):
  print "***", item, "***"
  print "    Interpolating grid."
  topogrid_tmp = fcn.interpZonalAdvanced(item,welon,eelon)
  topogrid = np.zeros((8400,15600))
  print "    N-S resizing for GEBCO."
  for i in range(8400):
    topogrid[i,:] = np.average(topogrid_tmp[i:i+2,:], axis=0)
  return topogrid


# Time step and file lists
files = np.array(sorted(glob.glob("topo*.txt")))
ages_numeric = np.loadtxt('ages')

# Confine to 22ka (LGM-ish) to present
lgm_to_present = (ages_numeric < 22) * (ages_numeric >= 0)
files = files[lgm_to_present]
ages_numeric = ages_numeric[lgm_to_present]

# turn "files" back into a list
#files = list(files)

# Generate output names with ages
ages = []
for item in ages_numeric:
  strage = '%06.2f' %item
  strage = strage[:3] + '_' + strage[4:] + 'k' # Replace decimal with underscore
  ages.append(strage)
ages = np.array(ages)

# Start out by gridding at time = 0 (for corrections to modern topography)
#topo_spharm_now = interpolate(files[ages_numeric == 0])
print "Loading modern topography derived from GEBCO_08 (30 arcsecond) and"
print "interpolated etopo1 bedrock surface data"
toponow = loadmat('../ICE3G_custom_uniform_grid/gebco_plus_etopo/gebco_and_etopo.mat')['map_data'] # This is the one local idiosyncracy I think I have
print "Generating array to correct spherical harmonic bumpiness"
correction = toponow - interpolate(files[ages_numeric == 0][0])

outgrid = garray.array()
for item in files:
  starttime = time.time()
  name = 'topo_' + ages[files == item][0]
  print "***", item, ":", name, "***"
  outgrid[...] = interpolate(item) + correction
  print "    Writing numpy array to GRASS GIS."
  outgrid.write(name)
  print "    ", time.time() - starttime, "seconds elapsed."





