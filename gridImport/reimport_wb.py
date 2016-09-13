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

# Time step and file lists
files = np.array(sorted(glob.glob("topo*.txt")))

ages = []
for fname in sorted(grass.parse_command('g.list', type='raster', pattern='topo_0*').keys()):
  ages += re.findall('\d+', fname)
ages = np.array(ages)
ages_numeric = ages.astype(int)

# Climate model run (TraCE-21K) goes back only to 22 ka
ages = ages[ages_numeric <= 22000]
ages_numeric = ages_numeric[ages_numeric <= 22000]

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
    
