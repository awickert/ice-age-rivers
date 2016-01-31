from grass import script as grass
from Scientific.IO.NetCDF import NetCDFFile
import numpy as np
import shutil

# Thinning prevents double-counting, but I also want to make sure that I don't miss cells that get to coast in between thinned areas.
# Will have to think about this again.

# Lats and lons for export
ncfile = NetCDFFile('qrparm.waterfix.hadcm3_bbc15ka.nc', 'r')

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

for age in ages:
  print ""
  print "***"
  print age
  print "***"
  print ""
  grass.run_command('g.region', rast='toponow')
  grass.mapcalc('tmp4 = ocean_plus_shore_'+age+' + '+'sl_binary_'+age+' * 0', overwrite=True)
  grass.run_command('r.grow', input='tmp4', output='tmp5', old=-1, new=1, overwrite=True)
  grass.run_command('r.colors', map='tmp5', color='rainbow')
  grass.run_command('r.null', map='tmp5', setnull=-1)
  grass.run_command('r.thin', input='tmp5', output='tmp6', overwrite=True)
  grass.mapcalc("tmp7 = tmp6 * accumulation_ice_"+age, overwrite=True)
  grass.run_command("discharge_to_coast_"+age+" = tmp7", overwrite=True) # anything negative reaching the coast just means no water, doesn't take away from rivers
  grass.mapcalc("discharge_to_coast_"+age+" = tmp7 * (tmp7 > 0)", overwrite=True) # anything negative reaching the coast just means no water, doesn't take away from rivers
  # Not setting these new zeros to null: 0 is where there is potential discharge but isn't any
  #Commenting out this part b/c no longer just providing this gridded output
  #grass.run_command('g.region', res='0:30')
  #grass.run_command('r.resamp.stats', method='sum', input='tmp8', output='tmp9', overwrite=True)
  grass.run_command('r.null', map='discharge_to_coast_'+age, setnull=0) # speeds up vector map creation
  grass.run_command('r.to.vect', input='discharge_to_coast_'+age, output='discharge_to_coast_'+age, type='point', column='discharge_m3_s', overwrite=True)
  grass.run_command('v.db.addcolumn', map='discharge_to_coast_'+age, columns='sea_grid_lon double precision, sea_grid_lat double precision')
  grass.run_command('v.distance', _from='discharge_to_coast_'+age, to='sea_grid_points', upload='to_x,to_y', column='sea_grid_lon,sea_grid_lat')

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
    discharge_grid[lats == round(row[-1],3), lons == round(row[-2],2)] += row[1]
  t.append(age)
  Qout.append(discharge_grid.copy()) # "copy" important! otherwise points to zeroed grid.
  print "Total discharge = ", np.sum(Qout[-1])
  discharge_grid *= 0

# Towards netcdf export
#shutil.copyfile('qrparm.waterfix.hadcm3_bbc15ka.nc', dst)
newnc = NetCDFFile('Q_m3_s.nc', 'w')
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

