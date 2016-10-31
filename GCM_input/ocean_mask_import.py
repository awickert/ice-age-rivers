# Land--ocean mask upload

#from grass import script as grass
from Scientific.IO.NetCDF import NetCDFFile
import numpy as np
import shutil
import os
import sys
import subprocess
import glob

# specify (existing) location and mapset
#location = "G12_NA"
gisdb    = "/data4/grassdata"
location = "ocean_model_grid"
mapset   = "PERMANENT"
#gisdb    = "/media/awickert/data4/grassdata"

# path to the GRASS GIS launch script
# MS Windows
grass7bin_win = r'C:\OSGeo4W\bin\grass71svn.bat'
# uncomment when using standalone WinGRASS installer
# grass7bin_win = r'C:\Program Files (x86)\GRASS GIS 7.0.0beta3\grass70.bat'
# Linux
grass7bin_lin = 'grass71'
# Mac OS X
# this is TODO
grass7bin_mac = '/Applications/GRASS/GRASS-7.1.app/'
 
# DATA
# define GRASS DATABASE
# add your path to grassdata (GRASS GIS database) directory
#gisdb = os.path.join(os.path.expanduser("~"), "grassdata")
# the following path is the default path on MS Windows
# gisdb = os.path.join(os.path.expanduser("~"), "Documents/grassdata") 
 
########### SOFTWARE
if sys.platform.startswith('linux'):
    # we assume that the GRASS GIS start script is available and in the PATH
    # query GRASS 7 itself for its GISBASE
    grass7bin = grass7bin_lin
elif sys.platform.startswith('win'):
    grass7bin = grass7bin_win
else:
    raise OSError('Platform not configured.')
 
# query GRASS 7 itself for its GISBASE
startcmd = [grass7bin, '--config', 'path']
 
p = subprocess.Popen(startcmd, shell=False,
                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()
if p.returncode != 0:
    print >>sys.stderr, "ERROR: Cannot find GRASS GIS 7 start script (%s)" % startcmd
    sys.exit(-1)
gisbase = out.strip('\n\r')
 
# Set GISBASE environment variable
os.environ['GISBASE'] = gisbase
# the following not needed with trunk
os.environ['PATH'] += os.pathsep + os.path.join(gisbase, 'extrabin')
# add path to GRASS addons
home = os.path.expanduser("~")
os.environ['PATH'] += os.pathsep + os.path.join(home, '.grass7', 'addons', 'scripts')
 
# define GRASS-Python environment
gpydir = os.path.join(gisbase, "etc", "python")
sys.path.append(gpydir)
 
########### DATA
# Set GISDBASE environment variable
os.environ['GISDBASE'] = gisdb
 
# import GRASS Python bindings (see also pygrass)
import grass.script as gscript
import grass.script.setup as gsetup

###########
# launch session
gsetup.init(gisbase,
            gisdb, location, mapset)
 
gscript.message('Current GRASS GIS 7 environment:')
print gscript.gisenv()


############################
# MAIN PROGRAM STARTS HERE #
############################

from grass import script as grass
from grass.script import array as garray

ncfiles_for_mask__name = sorted(glob.glob('/home/awickert/Dropbox/Papers/InProgress/26to0kaHadCM3/GIS/ocean_masks/*.nc'))

i = 0
for ncfile_for_mask__name in ncfiles_for_mask__name:
  # VERY SPECIFIC TO THIS DATA SET -- IS HARD-CODED HERE BUT WILL NEED TO BE CHAGNED IN GENERAL
  ncfile_basename = os.path.basename(ncfile_for_mask__name).split('.')[0] 
  # Import
  ncfile_for_mask = NetCDFFile(ncfile_for_mask__name, 'r')
  lsm = ncfile_for_mask.variables['lsm'].getValue()
  oceanmask = np.flipud(lsm)
  oceanmask[oceanmask == 1] = np.nan # remove land
  # Region
  lats = ncfile_for_mask.variables['latitude'][:]
  dlat = np.mean(np.diff(lats)) # bug possible here, but should work for regular grids
  s = lats[0] - dlat/2
  n = lats[-1] + dlat/2
  lons = ncfile_for_mask.variables['longitude'][:]
  dlon = np.mean(np.diff(lons)) # bug possible here, but should work for regular grids
  w = lons[0] - dlon/2
  e = lons[-1] + dlon/2
  grass.run_command('g.region', w=w, e=e, n=n, s=s, rows=oceanmask.shape[0], cols=oceanmask.shape[1])
  # ASSUMING NO CHANGE IN GRID SIZE!
  try:
    grass.mapcalc('lon = x()', overwrite=True)
  except:
    pass
  try:
    grass.mapcalc('lat = y()', overwrite=True)
  except:
    pass
  GRASSmask = garray.array()
  GRASSmask[...] = oceanmask
  GRASSmask.write(ncfile_basename, overwrite=True)
  grass.run_command('r.to.vect', input=ncfile_basename, output='sea_grid_points_'+ncfile_basename, type='point', overwrite=True)
  grass.run_command('v.db.dropcolumn', map='sea_grid_points_'+ncfile_basename, columns='value')
  grass.run_command('v.db.addcolumn', map='sea_grid_points_'+ncfile_basename, columns='lon double precision, lat double precision')
  grass.run_command('v.what.rast', map='sea_grid_points_'+ncfile_basename, raster='lon', column='lon')
  grass.run_command('v.what.rast', map='sea_grid_points_'+ncfile_basename, raster='lat', column='lat')
  
  i += 1
  
  print ""
  print "*****"
  print i/float(len(ncfiles_for_mask__name)) * 100, '%'
  print "*****"
  print ""

