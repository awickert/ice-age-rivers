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
import fnmatch

inloc = 'NA30asANU'
#inloc = 'Gregoire2012ice_latlon'
#grass.run_command('g.region', n='85N', s='26:45:30.019213N', w='175W', e='6:14:52.615673E', rows=287, cols=261) # G12

grass.run_command('g.region', n=85, s=15, w=-170, e=-40, res='0:0:30')

allmaps = grass.parse_command('r.proj', loc=inloc, flags='l').keys()

icemaps = sorted( fnmatch.filter(allmaps, 'ice_[0-9][0-9]*') )

for icemap in icemaps:
  age = re.findall('\d+', icemap)
  age = ''.join(age)
  if inloc == 'NA30as': # ICE-5G
    age += '00'
  if inloc == 'NA30asANU':
    age += '0'
  print icemap, '-->', 'ice_'+age
  grass.run_command('r.proj', location=inloc, input=icemap, output='ice_'+age) # G12
  grass.run_command('r.proj', location=inloc, input=icemap, output='ice_raw_'+age) # G12

"""
# Fix error
for icemap in icemaps:
  age = re.findall('\d+', icemap)
  age = ''.join(age)
  if inloc == 'NA30asANU':
    age += '0'
  grass.run_command('g.rename', rast='ice_raw_'+age+',ice_'+age)
"""

# G12
for age in ['000400', '000300', '000200', '000100', '000000']:
  grass.run_command('g.copy', rast='ice_raw_000500,ice_raw_'+age)

# G12
# /home/awickert/RUNS_G12/outaw090CVM2VM2_gregoireout
icemapsNew = sorted(glob.glob('*.txt'))
icemapsNew.append('topo021000.txt')
agesNew = []
for fname in icemapsNew:
  agesNew.append(re.findall('\d+', fname)[0])

# /home/awickert/Dropbox/Gregoire2012_ice_model/latlon_grids
icemapsOld = sorted(glob.glob('*.txt'))
  
for i in range(1, len(icemapsNew)+1):
  print icemapsOld[-i], 'ice_global_'+agesNew[-i]+'.txt'
  os.rename(icemapsOld[-i], 'ice_global_'+agesNew[-i]+'.txt')

# Before interpolation step  
ice_raw_import = sorted( grass.parse_command('g.list', type='raster', pattern='ice_raw_??????').keys() )
ages = []
ice = []
for item in ice_raw_import:
  ages.append(re.findall('\d+', item)[-1])
  ice.append('ice_' + ages[-1])

"""
icemapsNew = sorted(glob.glob('*.txt'))
for fname in icemapsNew:
  agesNew.append(re.findall('\d+', fname)[0])

grass.run_command('g.region', flags='p', 
"""
# Not sure if I trust these maps. Back to the initial NetCDF and project.
# Using code from above. But at least I renumbered these files!


# For ICE-6G
# Get ICE-6G files
# First, navigate to ICE-6G folder.
# Then...
files = sorted(glob.glob('*.nc'))
from Scientific.IO.NetCDF import NetCDFFile

grass.run_command('g.region', n=90, s=-90, w=0, e=360, res=1)
outarray = garray.array()
for fname in files:
  ncfile = NetCDFFile(fname)
  H_ice = ncfile.variables['stgit'][:][::-1]
  outarray[...] = H_ice
  age = re.findall('\d+', fname)[-1]
  outname = 'ice_raw_'+age
  print outname
  outarray.write(outname, overwrite=True)
    
grass.run_command('g.region', n=85, s=15, w=-170, e=-40, res='0:0:30')
ice_raw_import = sorted( grass.parse_command('g.list', type='raster', pattern='ice_raw_??????').keys() )
ages = []
ice = []
for item in ice_raw_import:
  ages.append(re.findall('\d+', item)[-1])
  ice.append('ice_' + ages[-1])
ages = np.array(ages)
ice = np.array(ice)
ages_numeric = ages.astype(int)
ice = ice[ages_numeric <= 22000]
#def resample_ice_clever(self):
# Using for Gregoire et al. (2012) model and ICE-6G
# Iterative interp for earlier models (especially with the blocky ICE-5G, cylinders of ICE-3G, ANU just along for the ride)

# ICE-6G and G12
for i in range(len(ice)):
  grass.run_command('g.region', n=85, s=15, w=-170, e=-40, res='0:0:30')
  #grass.run_command('r.resamp.bspline', input=ice_raw_import[i], output=ice[i], method='bilinear', lambda_=0.8, memory=10000, ns_step=0.25, ew_step=0.25, overwrite=True)
  print ice[i]
  grass.run_command('r.resamp.interp', input=ice_raw_import[i], output=ice[i], method='bicubic', overwrite=True)
  grass.mapcalc(ice[i]+' = '+ice[i]+' * ('+ice[i]+' > 0)', overwrite=True)

