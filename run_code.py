# Just make a set of plots
import os
os.chdir('/home/awickert/Dropbox/Papers/InProgress/DrainageMethods/program')
import drainage
d = drainage.Drainage(ICE=None) # Instantiate class
d.modelSetup()
d.modelDrainageSetup()
d.modelOutput()

# Standard (non-interactive) mode
import os
os.chdir('/home/awickert/Dropbox/Papers/InProgress/DrainageMethods/program')
import drainage
d = drainage.Drainage(ICE=None) # Instantiate class
#d.modelSetup(n=3500000, s=-3500000, w=-3500000, e=3500000, res=1000)
d.modelSetup(n=85, s=-60, w=-180, e=180, res='00:00:30')
d.modelClimate()
d.truncate_ages()
d.modelDrainageSetup()
d.modelDrainageSurface() # Often gets commented out on re-runs
d.modelDrainageFlowAccumulation()
d.modelDrainageLargeRivers()
d.modelDrainageBasins(rebuild_discharge_at_mouths=False)
d.modelDrainageBasinDischarge()
d.modelOutput()
"""

# Patch application a posteriori
for age in self.ages:
  grass.run_command('r.null', map='ocean_binary_'+age, null='0') # Turns ocean to null
self.compute_drainage.grow_ocean(self)
self.compute_drainage.vectorize_ocean_plus_shore(self)
self.compute_drainage.mouths(self)
self.compute_drainage.discharge_at_mouths(self)
d.modelDrainageBasins(rebuild_discharge_at_mouths=False)
d.modelDrainageBasinDischarge()
d.modelOutput(save=False)
"""

 # For interactive mode
reload(drainage)
#reload(d.initialize)
reload(d.setup)
reload(d.climate)
reload(d.compute_drainage)
reload(d.output)

# Starting interactive mode
import os
os.chdir('/home/awickert/Dropbox/Papers/InProgress/DrainageMethods/program')
import drainage
d = drainage.Drainage(ICE=None)
self = d
import sys
import re
from glob import glob
import numpy as np
import grass.script as grass
from grass.script import mapcalc
from grass.script import array as garray
from grass.script import db
from grass.script import vector as vect
from matplotlib import pyplot as plt
#d.modelSetup(n=85, s=-60, w=-180, e=180, res='00:00:30')
d.modelSetup()
d.climate.define_some_lists(d)
d.truncate_ages() # IMPORTANT! Should have this inside the functions, but hasn't been too much of a problem out here yet
#d.output.basin_discharge_save(d)

"""
d.ages = d.ages[1:] # IMPORTANT! Should have this inside the functions, but hasn't been too much of a problem out here yet
d.ages_numeric = self.ages_numeric[1:]
d.dt_numeric = np.diff(d.ages_numeric[::-1])[::-1] # update
d.dt = d.dt_numeric.astype(str)
d.midpoint_age = d.midpoint_age[1:] # first ts removed here as well
d.compute_drainage.printstart(d)
"""

# Python notes (because I was learning it while writing this)
# (and therefore this only uses classes half-heartedly)
# 1st "d" is for module fcn definition
# 2nd "d" is to pass the class to the function

#####################################################
# START OUT IN A REGION WITH INPUT DEMs PRE-DEFINED #
#####################################################

"""
for age in self.ages:
  grass.run_command('r.null', map='ice_'+age, null=0)
# Setup: Import the list of time steps
d.setup.start(d)
d.setup.setConstants(d)
d.setup.generateAges(d) # Used for the kk_h ages -- or removes padding (newer versions)
d.climate.define_some_lists(d)
"""

# Ice resampling, if needed:
d.climate.resample_ice(d)

# Water balance to produce input grids for runoff calculations
d.climate.runoff_input_meteoric(d)
d.climate.runoff_input_ice(d)
d.climate.runoff_input_total(d)

# Compute drainage
# Removing first time-step because we cannot calculate drainage routing at first ice time-step
# (No earlier dH_i/dt; and assuming for modern that dH_i/dt is constant... which is becoming
# not true, but not yet on a glacial cycle magnitude
#d.ages = d.ages[1:] # IMPORTANT! Should have this inside the functions, but hasn't been too much of a problem out here yet
#d.ages_numeric = self.ages_numeric[1:]
d.truncate_ages(d)
d.dt_numeric = np.diff(d.ages_numeric[::-1])[::-1] # update
d.dt = d.dt_numeric.astype(str)
d.midpoint_age = d.midpoint_age[1:] # first ts removed here as well
d.compute_drainage.printstart(d)
# Generate flow-routing grid (Z_r)
d.compute_drainage.flow_routing_grid_withocean(d, subglacial=False)
d.compute_drainage.separate_oceans(d)
d.compute_drainage.flow_routing_grid(d, subglacial=False)
d.compute_drainage.apply_etopo2_colormap(d)
# Use r.watershed to route flow and accumulate it for meteoric water, ice melt, and both
# My use of it here is not efficient, but is robust
d.compute_drainage.flow_routing_r_watershed(d)
d.compute_drainage.accum_nulls(d)
d.compute_drainage.flow_accum_ice(d)
d.compute_drainage.flow_accum_meteoric(d)
# Start to define rivers and drainage basins
d.compute_drainage.big_rivers(d)
d.compute_drainage.vectorize_streams(d)
d.compute_drainage.grow_ocean(d)
d.compute_drainage.vectorize_ocean_plus_shore(d)
d.compute_drainage.mouths(d)
#     Start here if river_mouth_regions have been modified -- remember to update "ages" (above)!
d.compute_drainage.discharge_at_mouths(d)
d.compute_drainage.build_basin_outlets(d)
d.compute_drainage.check_for_duplicate_outlets(d)
d.compute_drainage.build_basins_rast(d)#, river_name='Mackenzie')
d.compute_drainage.basins_to_null_int(d)
d.compute_drainage.build_basins_vect(d)
d.compute_drainage.add_basins_rast(d)
d.compute_drainage.build_basins_vect(d) # Have to do this twice -- interaction between rast & vect stored on disk.
#     Discharge and ice volumes in basins
d.compute_drainage.basin_discharge(d)
d.compute_drainage.build_basins_vect(d) # For some reason, must re-do this
d.compute_drainage.basin_discharge(d)
#d.compute_drainage.basin_ice_volume(d)
#save=False; show=True; legend=False; bigtitle=False; alltimes=False; onefig=True; ICE=None
d.output.basin_discharge_save(d)
d.output.basin_discharge_plots(d, save=True, show=True, legend=False, bigtitle=False, alltimes=False, onefig=True)

d.compute_drainage.ice_volume_sum(d)

"""
# Typically not used
#d.compute_drainage.drainage_basins_AND(d)
#d.compute_drainage.Q_i_unsmoothed(d)
#d.compute_drainage.rivers_by_area(d)
#d.compute_drainage.flood_basins(d)
#d.compute_drainage.basins(d)
#d.compute_drainage.lakes_ice_eq(d)
#d.compute_drainage.lakes_volume(d)
#d.compute_drainage.lakes_volume_sum(d)
# ADDED MARCH 2011 FOR ICE!
#d.compute_drainage.ice_volume(d)
#d.compute_drainage.ice_volume_sum(d)

# Update geophysics



# Update lakes again - make sure that it works



# Finalize
"""

