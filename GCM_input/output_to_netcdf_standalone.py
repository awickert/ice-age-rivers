# 23 November 2016
# Explicitly start GRASS GIS before running this.
# Writes pre-calculated grids to NetCDF

from Scientific.IO.NetCDF import NetCDFFile
import numpy as np
import shutil
import os
import sys
import subprocess
import glob
import sys
from grass import script as grass

location = str(grass.parse_command('g.gisenv', get='LOCATION_NAME').keys()[0])

def exists(name, datatype):
  YN = len(grass.parse_command('g.list', type=datatype, pattern=name))
  return YN

ncfiles_for_mask__name = sorted(glob.glob('/home/awickert/Dropbox/Papers/InProgress/26to0kaHadCM3/GIS/ocean_masks/*.nc'))

# Define ages
ops = sorted(grass.parse_command('g.list', type='rast', patt='ocean_plus_shore_*').keys())
ages = []
#agesint = []
for opsi in ops:
  ages.append(opsi[-6:])
  #agesint.append(int(opsi[-6:]))


#################
# OUTPUT NETCDF #
#################

ii = 0
for ncfile_for_mask__name in ncfiles_for_mask__name:

  # VERY SPECIFIC TO THIS DATA SET -- IS HARD-CODED HERE BUT WILL NEED TO BE CHAGNED IN GENERAL
  ncfile_basename = os.path.basename(ncfile_for_mask__name).split('.')[0] 
  # Output file name
  outname = '/home/awickert/Dropbox/Papers/InProgress/26to0kaHadCM3/MeltwaterInputFiles/'+location+'_'+ncfile_basename+'_Q_m3_s.nc'
  # Input GIS ocean grid: I have to import this *and*
  sea_grid_points = 'sea_grid_points_'+ncfile_basename

  if os.path.isfile(outname):
    print ""
    print "*****"
    print 'Skipping ocean mask', ncfile_basename
    print 'Output file already exists.'
    print "*****"
    print ""
  else:
    print ""
    print "*****"
    print 'Using ocean mask', ncfile_basename
    print "*****"
    print ""

  Qout = []
  t = []

  for age in ages[:10]:
    print ""
    print "***"
    print age
    print "***"
    print ""

    ncfile_for_mask = NetCDFFile(ncfile_for_mask__name, 'r')

    testgrid = 1 - ncfile_for_mask.variables['lsm'].getValue()
    testgrid[testgrid == 0] = np.nan

    # lsm = land-sea mask. 0 over ocean, 1 over land.
    discharge_grid = 0 * ncfile_for_mask.variables['lsm'].getValue() # Just a grid of the right size
                                                                     # Masking done in GRASS GIS location
    lats = ncfile_for_mask.variables['latitude'].getValue()
    elons = ncfile_for_mask.variables['longitude'].getValue()
    #lons = elons.copy()
    #lons[lons>180] -= 360 -- could test for this, but know that this set has 0--360 in Qll

    # Any possible floating point precision issues can be rounded;
    Qll = np.array(grass.vector_db_select('discharge_to_coast_'+ncfile_basename+'_'+age).values()[0].values(), dtype=float)
    
    for row in Qll:
      # Summing these into gridded bins here
      discharge_grid[lats == round(row[-1],3), elons == round(row[-2],2)] += row[1]
      #discharge_grid[np.asarray(np.ix_(lats == round(row[-1],3), lons == round(row[-2],2))).squeeze()]

    t.append(age)
    Qout.append(discharge_grid.copy()) # "copy" important! otherwise points to zeroed grid.
    print "Total discharge = ", np.sum(Qout[-1]), 'm3/s'
    discharge_grid *= 0

  # Towards netcdf export
  #newnc = NetCDFFile('test2.nc', 'w')
  #shutil.copyfile('qrparm.waterfix.hadcm3_bbc15ka.nc', dst)
  newnc = NetCDFFile(outname, 'w')
  newnc.createDimension('t', len(t))
  newnc.createDimension('longitude', len(lons))
  newnc.createDimension('latitude', len(lats))
  newnc.createVariable('t', 'i', ('t',))
  newnc.createVariable('longitude', 'f', ('longitude',))
  newnc.createVariable('latitude', 'f', ('latitude',))
  newnc.createVariable('discharge', 'd', ('t', 'latitude', 'longitude'))
  newnc.variables['t'][:] = t
  newnc.variables['longitude'][:] = elons
  newnc.variables['latitude'][:] = lats
  newnc.variables['discharge'][:] = np.array(Qout)
  newnc.close()

  ii += 1

print ""
print "*****"
print ii/float(len(ncfiles_for_mask__name)) * 100, '%'
print "*****"
print ""

