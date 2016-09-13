#from grass import script as grass
from Scientific.IO.NetCDF import NetCDFFile
import numpy as np
import shutil
import os
import sys
import subprocess

# specify (existing) location and mapset
#location = "G12_NA"
gisdb    = "/media/awickert/data4/grassdata"
#location = "ICE6G_Eurasia"
location = "ICE6G_globalNoAnt"
mapset   = "PERMANENT"
#gisdb    = "/media/awickert/data4/grassdata"

# Output file name
#outname = 'Eurasia_Q_m3_s.nc'
outname = 'GlobalNoAnt_ICE6G_Q_m3_s.nc'

# Outside of GRASS

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

""" 
gscript.message('Available raster maps:')
for rast in gscript.list_strings(type = 'rast'):
    print rast
 
gscript.message('Available vector maps:')
for vect in gscript.list_strings(type = 'vect'):
    print vect
""" 


# Thinning prevents double-counting, but I also want to make sure that I don't miss cells that get to coast in between thinned areas.
# Will have to think about this again.

from grass import script as grass

os.chdir('/home/awickert/Dropbox/Papers/InProgress/MWP1a_Runoff_Circulation/program')

# Lats and lons for export
ncfile = NetCDFFile('qrparm.waterfix.hadcm3_bbc15ka.nc', 'r')

grass.run_command('g.region', rast='topo_000000')

discharge_grid = 0 * ncfile.variables['field672'].getValue()[0,0]
lats = ncfile.variables['latitude'].getValue()
elons = ncfile.variables['longitude'].getValue()
lons = elons.copy()
lons[lons>180] -= 360

ops = sorted(grass.parse_command('g.list', type='rast', patt='ocean_plus_shore_*').keys())
ages = []
agesint = []
for opsi in ops:
  ages.append(opsi[-6:])
  agesint.append(int(opsi[-6:]))

Qout = []
t = []

# sea_grid_points
"""
try:
  grass.mapcalc('lon = x()')
except:
  pass
try:
  grass.mapcalc('lat = y()')
except:
  pass
"""

isll = grass.parse_command('g.proj', flags='g')['epsg'] == '4326'

#grass.run_command('r.to.vect', input='lon', output='sea_grid_points', type='point')
#grass.run_command('v.db.addtable', map='sea_grid_points', columns='sea_grid_lon double precision, sea_grid_lat double precision')
try:
  grass.run_command('v.proj', location='ocean_model_grid', input='sea_grid_points')
except:
  pass

for age in ages:
  print ""
  print "***"
  print age
  print "***"
  print ""
  # Abs because negative just means that there is offmap drainage -- happens when it hits the coast, should flag 
  # to maintain this positive in fugure revisions.
  # Choose depending on whether I want meltwater or full runoff
  # grass.mapcalc('discharge_to_coast_'+age+' = abs('+'ocean_plus_shore_'+age+' + accumulation_'+age+')', overwrite=True)
  grass.run_command('g.region', rast='topo_000000')
  try:
    grass.mapcalc('discharge_to_coast_'+age+' = abs('+'ocean_plus_shore_'+age+' + accumulation_ice_'+age+')', overwrite=True)
    grass.run_command('r.null', map='discharge_to_coast_'+age, setnull=0) # speeds up vector map creation
  except:
    pass
  if isll:
    grass.run_command('g.region', w=-180, e=180)
  try:
    grass.run_command('r.to.vect', input='discharge_to_coast_'+age, output='discharge_to_coast_'+age, type='point', column='discharge_m3_s', overwrite=True)
  except:
    try:
      grass.run_command('g.region', rast='topo_000000')
      grass.mapcalc('discharge_to_coast_'+age+' = abs('+'ocean_plus_shore_'+age+' + accumulation_ice_'+age+')', overwrite=True)
      grass.run_command('r.null', map='discharge_to_coast_'+age, setnull=0) # speeds up     pass
      grass.run_command('g.region', w=-180, e=180)
      grass.run_command('r.to.vect', input='discharge_to_coast_'+age, output='discharge_to_coast_'+age, type='point', column='discharge_m3_s', overwrite=True)
    except:
      print age, 'ERROR'
  try:
    # Top commented because....
    # This works only if it is lat/lon; will fail for projected grids
    #grass.run_command('v.distance', from_='discharge_to_coast_'+age, to='sea_grid_points', upload='to_x,to_y', column='sea_grid_lon,sea_grid_lat')
    # Uses stored values to work for projected grids
    try:
      tmp = grass.vector_db_select('discharge_to_coast_'+age, columns='sea_grid_lat').values()[0].values()
      if tmp[0] == ['']:
        if isll:
          grass.run_command('v.distance', _from='discharge_to_coast_'+age, to='sea_grid_points', upload='to_x,to_y', column='sea_grid_lon,sea_grid_lat')
        else:
          grass.run_command('v.distance', from_='discharge_to_coast_'+age, to='sea_grid_points', upload='to_attr', to_column='lon', column='sea_grid_lon')
          grass.run_command('v.distance', from_='discharge_to_coast_'+age, to='sea_grid_points', upload='to_attr', to_column='lat', column='sea_grid_lat')
    except:
      grass.run_command('v.db.addcolumn', map='discharge_to_coast_'+age, columns='sea_grid_lon double precision, sea_grid_lat double precision')
      if isll:
        grass.run_command('v.distance', _from='discharge_to_coast_'+age, to='sea_grid_points', upload='to_x,to_y', column='sea_grid_lon,sea_grid_lat')
      else:
        grass.run_command('v.distance', from_='discharge_to_coast_'+age, to='sea_grid_points', upload='to_attr', to_column='lon', column='sea_grid_lon')
        grass.run_command('v.distance', from_='discharge_to_coast_'+age, to='sea_grid_points', upload='to_attr', to_column='lat', column='sea_grid_lat')
  except:
    pass # No discharge points!
print ""
print "***"
print "Output Step"
print "***"
print ""

for age in ages:
  print ""
  print "***"
  print age
  print "***"
  print ""
  Qll = np.array(grass.vector_db_select('discharge_to_coast_'+age).values()[0].values(), dtype=float) # Any possible point precision issues can be rounded;
  for row in Qll:
    # Summing these into gridded bins here
    discharge_grid[lats == round(row[-1],3), lons == round(row[-2],2)] += row[1]
  t.append(age)
  Qout.append(discharge_grid.copy()) # "copy" important! otherwise points to zeroed grid.
  print "Total discharge = ", np.sum(Qout[-1])
  discharge_grid *= 0

# Towards netcdf export
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

