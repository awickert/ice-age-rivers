#! /usr/bin/python

# Code to connect lakes to ice model code, somewhat crudely but 
# quickly (in human-time)
# Started 14 Nov 2011 by ADW

# IMPORT LIBRARIES
import os
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

# Running with GRASS GIS: SET ENVIRONMENT VARIALBLES
######################################################

# Don't need to do this if running from inside a GRASS
# session already
# Which I plan to do, at first at least...

"""
os.environ['GISBASE'] = "/usr/local/grass-6.4.2svn/"
path = 
export PATH="$PATH:$GISBASE/bin:$GISBASE/scripts"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$GISBASE/lib"
# for parallel session management, we use process ID (PID) as lock file number:
export GIS_LOCK=$$
# path to GRASS settings file
export GISRC="$HOME/.grassrc6"
# PROBABLY NEED SOMETHING LINKING TO THE MAPSET AS WELL
"""

class Drainage(object):
  """
  The "Drainage" object acts to access the modules containing the analysis 
  functions and to serve as a container for class-wide scoped variables.
  """
  
  # IMPORT MODULES
  import setup
  import climate
  import compute_drainage
  import output
  
  def __init__(self, ICE=None):
    self.ICE = ICE

  def modelSetup(self, n=None, s=None, w=None, e=None, res=None, before_GCM=False):
    # n,s,w,e,res defaults are for North America and are set in "setup.start"
    self.setup.start(self, n=n, s=s, w=w, e=e, res=res)
    self.setup.setConstants(self)
    self.setup.generateAges(self, before_GCM=before_GCM) # Used for the kk_h ages -- or removes padding (newer versions)

  def modelClimate(self):
    self.climate.define_some_lists(self)
    self.climate.resample_ice(self)
    self.climate.runoff_input_meteoric(self)
    self.climate.runoff_input_ice(self)
    self.climate.runoff_input_total(self)

  def truncate_ages(self):
    """
    Removing first time-step because we cannot calculate drainage routing at first ice time-step
    (No earlier dH_i/dt; and assuming for modern that dH_i/dt is constant... which is becoming
    not true, but not yet on a glacial cycle magnitude
    """
    # Trim first entry
    self.ages = self.ages[1:] # IMPORTANT! Should have this inside the functions, but hasn't been too much of a problem out here yet
    self.ages_numeric = self.ages_numeric[1:]
    
  def modelDrainageSetup(self):
    """
    Compute drainage
    """
    self.dt_numeric = np.diff(self.ages_numeric[::-1])[::-1] # update
    self.dt = self.dt_numeric.astype(str)
    self.midpoint_age = self.midpoint_age[1:] # first ts removed here as well
    # And then print that we are starting
    self.compute_drainage.printstart(self)

  def modelDrainageSurface(self):
    """
    Generate flow-routing grid (Z_r)
    """
    self.compute_drainage.flow_routing_grid_withocean(self, subglacial=False)
    self.compute_drainage.separate_oceans(self)
    self.compute_drainage.flow_routing_grid(self, subglacial=False)
    self.compute_drainage.apply_etopo2_colormap(self)
    
  def modelDrainageFlowAccumulation(self):
    """
    Use r.watershed to route flow and accumulate it for meteoric water, ice melt, and both
    My use of it here is not efficient, but is robust
    """
    self.compute_drainage.flow_routing_r_watershed(self)
    self.compute_drainage.accum_nulls(self)
    self.compute_drainage.flow_accum_ice(self)
    self.compute_drainage.flow_accum_meteoric(self)

  def modelDrainageLargeRivers(self):
    """
    Build large rivers
    """
    self.compute_drainage.big_rivers(self)
    self.compute_drainage.vectorize_streams(self)
    self.compute_drainage.grow_ocean(self)
    self.compute_drainage.vectorize_ocean_plus_shore(self)
    self.compute_drainage.mouths(self)
    self.compute_drainage.discharge_at_mouths(self)

  def modelDrainageBasins(self, rebuild_discharge_at_mouths=True):
    """
    Build drainage basins
    """
    if rebuild_discharge_at_mouths:
      # Do this again here by default in case river mouth regions have changed
      self.compute_drainage.discharge_at_mouths(self)
    self.compute_drainage.build_basin_outlets(self)
    self.compute_drainage.check_for_duplicate_outlets(self)
    self.compute_drainage.build_basins_rast(self)#, river_name='Mackenzie')
    self.compute_drainage.basins_to_null_int(self)
    self.compute_drainage.build_basins_vect(self)
    self.compute_drainage.add_basins_rast(self)
    self.compute_drainage.build_basins_vect(self) # Have to do this twice -- interaction between rast & vect stored on disk.
#     Discharge and ice volumes in basins

  def modelDrainageBasinDischarge(self):
    self.compute_drainage.basin_discharge(self)
    #self.compute_drainage.basin_ice_volume(self)
    
  def modelOutput(self, save=True):
    self.ICE = grass.parse_command('g.gisenv')['LOCATION_NAME']
    self.ICE = re.findall("[a-zA-Z0-9]+", self.ICE)[0]
    self.output.basin_discharge_plots(self, save=save)

