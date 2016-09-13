# For ICE-6G
# All top material from spherical_interpolation.py
from grass import script as grass
from grass.script import array as garray
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate as interp
from Scientific.IO.NetCDF import NetCDFFile
import re
import glob

def make_map(name, outname, res, outres=None, sdiv=30):
  grass.run_command('g.region', n=90, s=-90, w=-180, e=180, res=res)
  a = garray.array() # var of interest
  theta = garray.array() # colat
  phi = garray.array()   # e-lon
  a.read(name, null=np.nan)
  
  grass.mapcalc("lats = y()", overwrite=True)
  grass.mapcalc("lons = x()", overwrite=True)
  theta.read('lats')
  theta = 90 - theta # colat
  phi.read('lons')
  phi += 180 # e-lon

  theta *= np.pi/180.
  phi *= np.pi/180.

  theta1 = theta.ravel()
  phi1   = phi.ravel()
  a1     = a.ravel()

  theta1 = theta1[np.isnan(a1) == False]
  phi1 = phi1[np.isnan(a1) == False]
  a1 = a1[np.isnan(a1) == False]

  if outres:
    grass.run_command('g.region', n=90, s=-90, w=-180, e=180, res=outres)
    thetaout = garray.array() # colat
    phiout = garray.array()   # e-lon
    grass.mapcalc("outlats = y()", overwrite=True)
    grass.mapcalc("outlons = x()", overwrite=True)
    # Sloppy, but don't need these anymore
    thetaout.read('outlats')
    thetaout = 90 - thetaout # colat
    phiout.read('outlons')
    phiout += 180 # e-lon
    thetaout *= np.pi/180.
    phiout *= np.pi/180.
  else:
    thetaout=theta
    phiout=phi

  # If a1 too small, can produce bad results
  lut = interp.SmoothSphereBivariateSpline( theta1, phi1, a1, s=int(np.floor(len(theta1)/sdiv)), eps=1E-10)
  data_global = lut(thetaout[:,0], phiout[0,:])
  outarray = garray.array()
  outarray[...] = data_global
  outarray.write(outname, overwrite=True)
  grass.run_command('g.region', n=90, s=-90, w=-180, e=180, res=res)

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
    
lats = ncfile.variables['lat'][:]
lons = ncfile.variables['lon'][:]

for fname in files:
  age = re.findall('\d+', fname)[-1]
  name = 'ice_raw_'+age
  outname = 'ice_'+age
  res=1
  sdiv=1

make_map('ice_raw_021000', 'ice_021000', outres='1', sdiv=3)


  """
  outarray = garray.array()
  outarray[...] = data_global
  outarray.write(outname, overwrite=True)
  grass.run_command('g.region', n=90, s=-90, w=-180, e=180, res=res)

  
