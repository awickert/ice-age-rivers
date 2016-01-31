import os
import sys
from glob import glob
import numpy as np
from grass.script import core as grass
from grass.script import mapcalc
from grass.script import array as garray
from grass.script import db


# IMPORT FILES

wbdir = '/home/awickert/Documents/iceage_climate/water_balance_GRASS'
wbfiles = sorted(glob(wbdir + '/wb_*.txt'))

wb = []
for wbpath in wbfiles:
  wbfile = os.path.basename(wbpath)
  wb.append(wbfile[:-4])

for i in range(len(wbfiles)):
  grass.run_command('r.in.ascii', input=wbfiles[i], output=wb[i], overwrite=True)


# INTERPOLATE AT LOW RES TO GET ALL NEARSHORE NULL VALUES
age = []
for w in wb:
  age.append(w[3:])
  
for i in range(len(wb)):
  grass.run_command('g.rename', rast=wb[i] + ',wb_orig_' + age[i], overwrite=True)

reg = grass.region() # Get original region parameters
grass.run_command('g.region', flags='p', res=3.75)
for i in range(len(wb)):
#  os.system("r.fillnulls input=wb_orig_" + age[i] + " output=" + wb[i] + " --o")
  grass.run_command('r.fillnulls', input='wb_orig_' + age[i], output=wb[i], overwrite=True)



# NOW, RESAMPLE TO HIGH RES
# Actually, no resampling needed :) Just doing it nearest-neighbor style
# Actually, it is needed: areas of cells are important and need to be remembered 

for i in range(len(wb)):
  grass.run_command('g.rename', rast=wb[i] + ',wb_coarse_' + age[i], overwrite=True)

# This is the non-interpolating way
grass.run_command('g.region', flags='p', rast='cellsize_meters2')
for i in range(len(wb)):
  grass.run_command('r.resample', input='wb_coarse_' + age[i], output=wb[i], overwrite=True)


# Then mapcalc with the Mississippi drainage basins to get the meteoric contribution to the Mississippi
# And sum to get the time series to go along with the ice melt one
Mississippi = []
Meteoric = []
for i in range(len(age)):
  Mississippi.append("Mississippi_" + age[i])
  Meteoric.append("Meteoric_" + age[i])

# wb[i] = [mm/s]; cellsize_meters2 = [m**2]; Mississippi[i] = [-]; / 1000 turns [mm] --> [m] --> so [m**3/s]
for i in range(len(age)):
  mapcalc(Meteoric[i] + ' = ' + Mississippi[i] + ' * ' + wb[i] + ' * cellsize_meters2 / 1000' )

Qmeteoric = []
for i in range(len(age)):
  Qmeteoric.append(float(grass.parse_command('r.sum', rast=Meteoric[i])['SUM']))

# Output with ages

# Ages as ints
age_a = []
extra_zeros = '00'
for i in range(len(wb)):
  age_ka = age[i][:2]
  age_ca = age[i][-2]
  age_a_str = age_ka + age_ca + extra_zeros # string
  age_a.append(int(age_a_str))

mean_age_int = []
mean_age = []

for i in range(1,len(wb)):
  mean_age_int.append( (age_a[i-1] + age_a[i]) / 2 )

dt_int = np.diff(age_a)

dt = []
for i in range(len(dt_int)):
  dt.append(str(dt_int[i]))

for i in range(0,len(wb)-1):
  mean_age.append('%05d' %mean_age_int[i] )

# And save it

tQ = np.zeros((len(Qmeteoric),2))
tQ[:,0] = age_a
tQ[:,1] = Qmeteoric

np.savetxt('Mississippi_meteoric_discharge.txt',tQ,fmt='%.6f')


# Import Tom Marchitto's data, as corrected to calendar years by me
marchitto=np.genfromtxt('marchitto_tab_calyrs_estimate.txt',skip_header=0)
p_reworked=np.array([marchitto[:,5]*100]).transpose()
marchitto=np.concatenate((marchitto,p_reworked),axis=1)


# Plot against meltwater discharge
tQ = np.loadtxt('Mississippi_meteoric_discharge.txt')
Qmelt = np.loadtxt('Mississippi_meltwater_discharge.txt')

from scipy.interpolate import interp1d

f = interp1d(tQ[:,0], tQ[:,1])
Qmeteoric_mid = f(Qmelt[:,0])

from matplotlib import pyplot as plt

# INCORPORATES MARCHITTO, RITTENOUR, AND ME

# SUGGESTIONS:
# RED AXIS ON RHS - DONE
# SHADE LIGHTLY UNDER TOTAL DISCHARGE?
# POINTS ON LEGEND (FAKE IT BY ADDING SOMETHING WAY OFF THE SCREEN TO AX1)
# (MY THOUGHT) Make points bigger?

fig = plt.figure(1,figsize=(12,6))
ax1 = fig.add_subplot(111)

ax1.set_xlabel("Age [kyr BP]", fontsize=16)
ax1.set_ylabel("Mississippi River discharge [m$^3$/s]", fontsize=16)

ax2 = ax1.twinx()
ax2.tick_params(axis='y', colors='red')
ax2.spines['right'].set_color('red')
ax1.spines['right'].set_color('red')
ax2.plot(marchitto[:,1]/1000.,marchitto[:,6],'r.', linewidth=1, markersize=15, markeredgecolor = 'k', markerfacecolor='r', label="Reworked Nannofossil %")
ax2.set_ylabel('Reworked Calcareous Nannofossil %',color='r',fontsize=16)

ax1.plot(Qmelt[:,0]/1000,Qmelt[:,1],'b-', color='DarkTurquoise', label="Meltwater discharge", linewidth=2)
ax1.plot(tQ[:,0]/1000,tQ[:,1],'g-', label="Meteoric discharge", linewidth=2)
ax1.plot(Qmelt[:,0]/1000,Qmelt[:,1]+Qmeteoric_mid,'k-', label="Total discharge", linewidth=4)
ax1.plot(-1000000,50,'r.', linewidth=1, markersize=15, markeredgecolor = 'k', markerfacecolor='r', label="Reworked Nannofossil %") # Hidden: to put label on ax1 legend

# LMR incision around MWP 1a
rittenour = ax1.axvspan(12.4, 15.0, facecolor='y', alpha=0.7, label="Lower Mississippi River Incision")

plt.xlim((0,21))

ax1.legend(loc=2,numpoints=1)

plt.show()



# Get drainage basin area
Mississippi_cell_areas_km2 = []
for i in range(len(age)):
  Mississippi_cell_areas_km2.append('Mississippi_cell_areas_km2_' + age[i])

for i in range(len(age)):
  print age[i]
  mapcalc(Mississippi_cell_areas_km2[i] + ' = ' + Mississippi[i] + ' * cellsize_km2')

Area_km2 = []
for i in range(len(age)):
  Area_km2.append(float(grass.parse_command('r.sum', rast=Mississippi_cell_areas_km2[i])['SUM']))

np.savetxt('/home/awickert/Documents/AGU2011/Drainage/Area_km2.txt',Area_km2)


# Get ice coverage
Mississippi_ice_coverage = []
for i in range(len(age)):
  Mississippi_ice_coverage.append('Mississippi_ice_coverage_' + age[i])

ice = []
for i in range(len(wb)):
  ice.append('ice_' + age[i])

for i in range(len(age)):
  print age[i]
  mapcalc(Mississippi_ice_coverage[i] + ' = ' + Mississippi_cell_areas_km2[i] + ' * (' + ice[i] + ' > 0 )')

Area_km2_ice = []
for i in range(len(age)):
  Area_km2_ice.append(float(grass.parse_command('r.sum', rast=Mississippi_ice_coverage[i])['SUM']))

# Save it as well
np.savetxt('/home/awickert/Documents/AGU2011/Drainage/Area_km2_ice.txt',Area_km2_ice)




# Get temperature - just pick Iowa as averaging cold Rockies and warm lower Mississippi
# I did this with gcm_nc.py
# Picked 42.5N, 93W
TS = np.loadtxt('/home/awickert/Documents/AGU2011/Drainage/GCM_42.5N-93.0E/TS.txt')
tTS = np.loadtxt('/home/awickert/Documents/AGU2011/Drainage/GCM_42.5N-93.0E/time.txt')

tTSa = []
for tTSi in tTS:
  tTSa.append(round(-1000 * tTSi))

# Will throw an error, but works
TK = []
i = 0
j = -1
for tTSi in tTSa:
  if tTSi == age_a[j]:
    TK.append(TS[i])
    j -= 1
  i += 1

np.savetxt('/home/awickert/Documents/AGU2011/Drainage/T_Iowa.txt',T)


# DO QART (B=1) sediment load prediction and one that includes glacial area
  # Need to interpolate to midpoints for ice melt
  # "f" defined above

# Time on time steps
t_ts = tQ[:,0]
# (t is time at midpoints)

# Ignoring ice (warm vs. cold base?), lithology (don't have map... and ice interactions?)
# No human influence
B = 1
# With ice cover
fice_ts = np.array(Area_km2_ice) / np.array(Area_km2)
f = interp1d(t_ts, fice_ts)
fice = f(t) # proportion and interpolation
Bice = 1 + 9 * fice # Now an array; * 9 instead of .09 b/c I have 0-1 fract instead of 0-100 %

# Total discharge
Q = Qmeteoric_mid + Qmelt[:,1]

# Total drainage area
f = interp1d(t_ts, np.array(Area_km2))
A = f(t)

# Relief
R = 2500./1000 # km post-division, sort of representative for the super-high parts of the Missouri

# Temperature
T = np.array(TK) - 273.15 # to deg C
f = interp1d(t_ts, T)
T = f(t)

# Coefficient
omega = 0.02 # kg/s

# For equivalent dpeosit volume
omega_deprate = omega
lambda_p = 0.65
omega_deprate *= (1./2650.) * (1./(1.-lambda_p)) # m**3/s deposit
omega_deprate /= 1E9 # for km**3
omega_deprate *= 31556926. # to years
#omega_deprate /= 300000. # Mississippi fan minimum area, per Bouma 
#omega_deprate /= 40000. # Approx for large sedimentation area of recent Mississippi fan lobe, km**2
omega_deprate /= 80000. # Approx for whole area of recent Mississippi fan lobe, km**2
omega_deprate *= 1000. # Converting heights back to meters

# THE EQUATION
# T stuff is to make sure that it is "1" for T<2 and T for T>2:
Qs = omega * B * Q**0.31 * A**0.5 * R * (T * (T>=2) + (T<2))
etadot = omega_deprate * B * Q**0.31 * A**0.5 * R * (T * (T>=2) + (T<2))
# With ice
Qs_ice = omega * Bice * Q**0.31 * A**0.5 * R * (T * (T>=2) + (T<2))
etadot_ice = omega_deprate * Bice * Q**0.31 * A**0.5 * R * (T * (T>=2) + (T<2))


# Get cumulative deposition from 21 ka to present
cumdep_increments = etadot * dt_int
cumdep_ice_increments = etadot_ice * dt_int

cumdep = []
trev = []
for i in range(len(age)-1):
  cumdep.append(sum(cumdep_increments[-i-1:]))
  trev.append(t[-i-1])
  
cumulative_deposit = np.zeros((38,2))
cumulative_deposit[:,0] = trev
cumulative_deposit[:,0] = cumdep

cumdep_ice = []
for i in range(len(age)-1):
  cumdep_ice.append(sum(cumdep_ice_increments[-i-1:]))

# Plotting

trev = np.array(trev)

plt.figure(1,figsize=(12,6))

plt.plot(trev/1000., cumdep, 'k.-', label="No ice considered in sediment input", linewidth=2)
plt.plot(trev/1000., cumdep_ice, 'b.-', label="Alpine glacier style erosion", linewidth=2)

plt.xlim((0,21))

plt.xlabel("Age [kyr BP]", fontsize=16)
plt.ylabel("Cumulative Mississippi Offshore Aggradation [m]", fontsize=16)

plt.legend(loc=0)

plt.show()




t = Qmelt[:,0]
plt.plot(t,Qs)
plt.show()





# This way includes basic interpolation

# Step I: rename current files to wb_coarse
for i in range(len(wb)):
  grass.run_command('g.rename', rast=wb[i] + ',wb_coarse_' + age[i], overwrite=True)

# Step II: iteratively interpolate
for i in range(len(wb)):
  # First, set up the region resolution to approximately that of the climate model output
  grass.run_command('g.region', flags='p', rast='cellsize_meters2')
  grass.run_command('g.region', flags='p', res=3.75) # Climate model native resolution
  # Then resample and run r.neighbors as a "primer"
  grass.run_command('g.copy', rast='wb_coarse_' + age[i] + ',' + wb[i], overwrite=True)
  # Subsequent iterations with nearest-neighbor interpolation until we get past the size of the topography
  res = 3.75
  next_res = res / 3.
  while next_res > reg['nsres']: # Up to the semi-final needed resolution
    res /= 3.
    next_res /= 3.
    grass.run_command('g.region', flags='p', res=res) # Climate model native resolution
    grass.run_command('r.resample', input=wb[i], output=wb[i], overwrite=True)
    grass.run_command('r.neighbors', input=wb[i], output=wb[i], size=3, overwrite=True)
  # Final resolution
  grass.run_command('g.region', flags='p', rast='cellsize_meters2')
  grass.run_command('r.resample', input=wb[i], output=wb[i], overwrite=True)


# AFTER THIS, CALCULATE THE RUNOFF GENERATED PER CELL

runoff = []
for a in age:
  runoff.append('runoff_' + a)

# Forgot colon earlier, so will just skip to the bigger step below
for i in range(len(wb)):
  mapcalc(runoff[i] + " = " + wb[i] + " * cellsize_meters2")


# NOW, RE-DO DRAINAGE CALCULATIONS; ACCUMULATION MAP WILL BE TRUE WATER
# In units of mm * m**2, so divide by 1000 and you get m**3 in Q from that cell 
# From meteoric water

# Although, perhaps I should g.move rast=runoff,runoff_meteoric, and add 
# some estimate of ice discharge to these as well.
# Yes - calc ice mass balance for whole range (go in middle of times)
# And then average between middles to get back to nodes
# And then do water discharge
# With each of those cells generating their negative ice (converted to water) 
# volume balance in discharge

# OK - that's the plan! And it's generalizably applicable!
# ... so long as I don't have overflow OR round to 0 problems...

ice = []
for i in range(len(wb)):
  ice.append('ice_' + age[i])

# Ages as ints
age_a = []
extra_zeros = '00'
for i in range(len(wb)):
  age_ka = age[i][:2]
  age_ca = age[i][-2]
  age_a_str = age_ka + age_ca + extra_zeros # string
  age_a.append(int(age_a_str))

mean_age_int = []
mean_age = []

for i in range(1,len(wb)):
  mean_age_int.append( (age_a[i-1] + age_a[i]) / 2 )

dt_int = np.diff(age_a)

dt = []
for i in range(len(dt_int)):
  dt.append(str(dt_int[i]))

for i in range(0,len(wb)-1):
  mean_age.append('%05d' %mean_age_int[i] )

# VOLUME BALANCES INCORPORATE A LOSS-POSITIVE SIGN CONVENTION

dHi_dt = []
dVi_dt = []
for i in range(len(dt)):
  dHi_dt.append('dHi_dt_' + mean_age[i])
  dVi_dt.append('dVi_dt_' + mean_age[i])

for i in range(len(dt)):
  # Meters per year
  print dHi_dt[i]
  mapcalc(dHi_dt[i] + ' = (' + ice[i+1] + ' - ' + ice[i] + ') / ' + dt[i])

for i in range(len(dt)):
  # m * km**2 / yr
  print dVi_dt[i]
  mapcalc(dVi_dt[i] + ' = ' + dHi_dt[i] + ' * cellsize_km2')

dVi_dt_ts = []
for i in range(len(age)-1):
  dVi_dt_ts.append("dVi_dt_" + age[i])

# Use averaging to get back to our good old time step values to allow for intercomparison
for i in range(len(dVi_dt_ts)-1):
  # Weight the midpoint
  midpoint = age_a[i+1]
  recent = float(mean_age[i])
  older = float(mean_age[i+1])
  iweight = (midpoint - recent) / (older - recent)
  iplusweight = (older - midpoint) / (older - recent)
  if iweight != iplusweight:
    print dVi_dt_ts[i+1]
    mapcalc( dVi_dt_ts[i+1] + ' = (' + dVi_dt[i] + ' * ' + str(iweight) + ') + (' + dVi_dt[i+1] + ' * ' + str(iplusweight) + ')' )

# 250 and 0 a should both have no ICE-5G loss (or ice at all!) in AK
# (model ice, that is - real ice doesn't count)
# So we can just copy to add on the modern day
grass.run_command('g.copy', rast=dVi_dt[0] + ',' + dVi_dt_ts[0])

# Now I will have to convert ice into mm * m**2 / s and sum with surface water to get full runoff grid
# Also convert from ice volume to water volume
for i in range(len(dt)):
  # m * km**2 / yr
  print runoff[i]
  mapcalc( runoff[i] + " = (" + wb[i] + " * cellsize_meters2) + (" + dVi_dt_ts[i] + " * 1000 * 1000 * 1000 * 0.917 / 31556926)" )


# Note, tpi, tpin already defined from previous calcs, so I can just refrence them here

tpi = []
tpin = []
tpin_coastline = []
flooded = []
direction = []
accumulation = []
stream = []
continent_plus_coast = []
notocean_only = []
for i in range(len(age)): # "len dt" b/c I have no map to dictate overland flow input at LGM, and len(dt) is 1 less than full length
  # Actually, can do full "age": needed for consistent move (below) and doesn't hurt
  tpi.append("topo_plus_ice_" + age[i])
  tpin.append("topo_plus_ice_noocean_" + age[i])
  tpin_coastline.append("topo_plus_ice_coastline_" + age[i]) # tpin with water cells along shore included
  flooded.append("tpin_flooded_" + age[i]) # closed basins filled - USED ONLY FOR TERRAFLOW, BUT NEEDED FOR LAKES, UNLESS I USE ANOTHER ALGORITHM (r.flooded or something?)
  direction.append("flow_direction_" + age[i])
  accumulation.append("flow_accum_" + age[i])
  stream.append("stream_" + age[i])
  continent_plus_coast.append("continent_plus_coast_" + age[i])
  notocean_only.append("notocean_only_" + age[i])

threshold = '10000.0' # 10 m**3 / s threshold; I wonder if putting it with the decimal will make it float and remove floor division
# Seems that it is desired to be an int, though...
# 10000 b/c my units are mm m**2, not m**3

# tpin_coastline already calced, so don't have to do it
for i in range(len(dt)):
  mapcalc(tpin_coastline[i] = $tpi * $continent_plus_coast" # Now have shore-touching water cells included; r.watershed doesn't dive into nulls like r.terraflow did

# Move old files to "AREA_" for based only on drainage area
ages = age
for i in range(len(ages)):
  grass.run_command('g.rename', rast=accumulation[i] + ',AREA_' + accumulation[i])

for i in range(len(ages)):
  grass.run_command('g.rename', rast=direction[i] + ',AREA_' + direction[i])

for i in range(len(ages)):
  grass.run_command('g.rename', rast=stream[i] + ',AREA_' + stream[i])
  
# Calculate drainage here
import time
for i in range(len(dt)):
  print "Drainage calculations starting " + age[i]
  tstart = time.time()
  grass.run_command('r.watershed', elevation=tpin_coastline[i], flow=runoff[i], accumulation=accumulation[i], drainage=direction[i], stream=stream[i], threshold=threshold, overwrite=True) # Takes ~80 seconds on desktop, without $stream / $threshold portion, which (it turns out) doesn't take a whole lot more time... though outputting the maps might. SFD requires 'flag="s"' to be added if run in GRASS 7
  telapsed = time.time() - tstart
  print "%.1f" %telapsed + " seconds elapsed"


# Stream vector characteristics - junctions, coastlines
#########################################################

# Move old map - just 21 k
ages = age
for i in range(len(ages)):
  grass.run_command('g.rename', vect=tributary_junctions[i] + ',AREA_' + tributary_junctions[i])

ages = age
for i in range(len(ages)):
  grass.run_command('g.rename', vect=stream[i] + ',AREA_' + stream[i])
grass.run_command('g.remove', vect='AREA_stream_20_0k') # Already copied over this

stream_thinned = []
for i in range(len(ages)):
  stream_thinned.append("stream_thinned_" + age[i])

stream_nodes_raw = []
for i in range(len(ages)):
  stream_nodes_raw.append("stream_nodes_raw_" + age[i])

stream_nodes = []
for i in range(len(ages)):
  stream_nodes.append("stream_nodes_" + age[i])

tributary_junctions = []
for i in range(len(ages)):
  tributary_junctions.append('tributary_junctions_' + ages[i])

topo = []
for i in range(len(ages)):
  topo.append('topo_' + ages[i])



# Stream to vector
for i in range(len(ages)-1):
  print ages[i]
  grass.run_command('r.thin', input=stream[i], output=stream_thinned[i])
  grass.run_command('r.to.vect', input=stream_thinned[i], output=stream[i], overwrite=True)

# FIND TRIBUTARY JUNCTIONS
for i in range(len(ages)-1):
  # Get nodes at start and end of all segments
  print ages[i]
  grass.run_command('v.to.points', input=stream[i], output=stream_nodes_raw[i], flags='n', overwrite=True)
  
# Clean map: this removes duplicates (and could, w/ the right settings, I think, be used to limit the number of points in a particular area)
for i in range(len(ages)-1):
  print ages[i]
  grass.run_command('v.clean', input=stream_nodes_raw[i], output=stream_nodes[i], tool='rmdupl')

# Drop the tables in prepr to re-add: category problems
for i in range(len(ages)-1):
  grass.run_command('v.db.droptable', map=stream_nodes[i], layer=2, flags='f')
  grass.run_command('v.db.droptable', map=stream_nodes[i], layer=1, flags='f')

for i in range(len(ages)-1):
  grass.run_command('v.category', input=stream_nodes[i], opt='del', output='tmp')
  grass.run_command('v.category', input='tmp', opt='add', output=stream_nodes[i], overwrite=True)
  grass.run_command('g.remove', vect='tmp')

# Remove the nodes that are at the starts and ends of the river segments
for i in range(len(ages)-1):
  grass.run_command('v.db.addtable', map=stream_nodes[i], col="cat int, Accumulation double precision, FD_Elev double precision")
  # Upstream
  grass.run_command('v.what.rast', vect=stream_nodes[i], rast=accumulation[i], col='Accumulation', layer=1)
  # Downstream
  grass.run_command('v.what.rast', vect=stream_nodes[i], rast=tpi[i], col='FD_Elev')

# Then extract only those parts that are desired
for i in range(len(ages)-1):
  # basin_threshold in input file
  basin_threshold = threshold # nomenclature inconsistency above
  basin_threshold_plus_float = float(basin_threshold) + float(basin_threshold)/2. # guess number to add based on scaling in basin_threshold (10x km^2 for Beringia) and cell size and margin of error - though now just have multiple
  basin_threshold_plus = str(basin_threshold_plus_float)

wherestr = "Accumulation > " + basin_threshold_plus + " and FD_Elev > 0"
for i in range(len(ages)-1):
  grass.run_command('v.extract', input=stream_nodes[i], output=tributary_junctions[i], type='point', where=wherestr)

for i in range(len(ages)-1):
  grass.run_command('r.colors', map=topo[i], color='etopo2')

# Plot the tributary junctions
grass.run_command('d.mon', start='x3')
grass.run_command('d.rast', map=topo[i])
grass.run_command('d.vect', map=stream[i])
#d.vect $channel_nodes
grass.run_command('d.vect', map=tributary_junctions[i], color='blue')


# Coastline

# "Shore" already created
shore = []
shore_thinned = []
for i in range(len(ages)):
  shore.append('shore_' + ages[i])
  shore_thinned.append('shore_thinned_' + ages[i])

# Thin
for i in range(len(ages)):
  grass.run_command('r.thin', input=shore[i], output=shore_thinned[i])

# Vectorize
for i in range(len(ages)):
  grass.run_command('r.to.vect', input=shore_thinned[i], output=shore[i])

# Find river outlets with a v.clean trick
# See http://www.surfaces.co.il/?p=241
# Copied and pastied, then modified, part of their example
stream_shore = []
stream_toolong = []
pplocs = []
for i in range(len(ages)):
  stream_shore.append('stream_shore_' + ages[i])
  pplocs.append('pplocs_' + ages[i])
  stream_toolong.append('stream_toolong_' + ages[i])

# Guess they weren't created!
for i in range(len(ages)):
  grass.run_command('g.rename', vect=pplocs[i] + ',AREA_' + pplocs[i])

# Streams have dangles on their ends that need to be removed... 
# DON'T DO THIS NOW
#for i in range(len(ages)):
#  grass.run_command('g.rename', vect=stream[i] + ',' + stream_toolong[i])

for i in range(len(ages)-1):
  grass.run_command('v.patch', input=stream[i] + ',' + shore[i], output=stream_shore[i])
  
for i in range(len(ages)-1):
  grass.run_command('v.clean', input=stream_shore[i], output='tmp', tool='break', error=pplocs[i], overwrite=True)
  grass.run_command('g.remove', vect='tmp')


# Now, get the discharges at river mouths
for i in range(len(ages)-1):
  grass.run_command('v.db.addtable', map=pplocs[i], col="cat int, Accumulation double precision, Discharge double precision")

for i in range(len(ages)-1):
  grass.run_command('g.copy', vect=pplocs[i] + ',tmp', overwrite=True)
  grass.run_command('v.category', input='tmp', output=pplocs[i], overwrite=True)
  
for i in range(len(ages)-1):
  grass.run_command('v.to.db', map=pplocs[i], columns="cat", option="cat", type='point')
  
for i in range(len(ages)-1):
  print ages[i]
  grass.run_command('v.what.rast', vect=pplocs[i], rast=accumulation[i], col='Accumulation', layer=1)

for i in range(len(ages)-1):
  print ages[i]
  grass.write_command('db.execute', stdin = "UPDATE " + pplocs[i] + " SET Discharge=Accumulation/1000" )


# As well as discharges at confluences
for i in range(len(ages)-1):
  print ages[i]
  grass.run_command('v.db.addcol', map=tributary_junctions[i], col="Discharge double precision")

for i in range(len(ages)-1):
  print ages[i]
  grass.write_command('db.execute', stdin = "UPDATE " + tributary_junctions[i] + " SET Discharge=Accumulation/1000" )



# After this, need to create a set of points on all of the channel segments.
# These points will contain flow direction and discharge, allowing me to map 
# out water discharge down streams by exporting these in order
# How to get in order? Not sure yet...



# r.fillnulls input=wb_orig_00_0k output=wb_00_0k --o


# INTERPOLATE TO HIGH RES
# Some function goes here to resample and interpolate
# based on my ICE-5G interpolation
def resamp_interp:
  pass


# Maybe create threshold for the start of streams based on water 
# balance being greater than some value - not just drainage area 
# anymore

# Make vars
for age in self.ages:
  flow_contribution.append('flow_contribution_' + age)

# Written for Beringia, I think
# I had cellsize_km2 * 10
# My wb values range from 10^-5-ish to 10^-4 in really wet places
# So I think that I won't overflow
# Plus, this is MUCH smaller than all of North America!
# wb in [mm/s], so this is in units of [mm m^2 / s]
for i in range(len(self.ages)):
  mapcalc(flow_contribution[i] + ' = cellsize_meters2 * ' + wb[i])
  
# Use 10 mm*3/s as threshold!

# Accum is negative in areas that go offmap, should turn these into
# NULL for the rest of my analysis

# Nope, see above, need to do it before "accum" step!
# Accum has area_size already built in, roughly
accum_scale = '10' # Accumulation scale 
mapcalc(accumulation * 

# Within each basin, I will want to do some sort of water balance
for age in ages:
  wb = wb_[i] + age
  
  
  
