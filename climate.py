#! /usr/bin/python

"""
Module to do the water balance% import and interpolation
"""

from drainage import *

def define_some_lists(self):
  self.wb = []
  self.ice = []
  self.ice_raw_import = []
  self.runoff_input_ice = []
  self.runoff_input_meteoric = []
  self.runoff_input = []
  for age in self.ages:
    self.wb.append('wb_' + age)
    self.ice.append('ice_' + age)
    self.ice_raw_import.append('ice_raw_' + age)
    self.runoff_input_meteoric.append('runoff_input_meteoric_' + age)
  # 1: because this is before I start ignoring the first time-step
  for age in self.ages[1:]:
    self.runoff_input_ice.append('runoff_input_ice_' + age)
    self.runoff_input.append('runoff_input_' + age)

def resample_ice(self):
  print ""
  print "Resampling ice maps from global to regional higher-resolution, if necessary"
  print ""
  grass.run_command('g.region', region='default')
  # Iteratively interpolate - assume "ages" is a match
  for i in range(len(self.ice)):
    print self.ages[i]
    exists = len(grass.parse_command('g.list', type='raster', pattern=self.ice[i]))
    if exists:
      print "  Resampled ice map already exists; skipping."
    else:
      # First, set up the region resolution to approximately that of the climate model output
      grass.run_command('g.region', region='default') # Ice model native resolution, with some padding
      reg = grass.region() # Get original region parameters - original for the region, final for the climate model
      grass.run_command('g.region', res=1) # Ice model native resolution, with some padding # <------------------------------ NEED TO GENERALIZE IF NOT 1 DEGREE!
      grass.run_command('g.region', n=reg['n']+1, s=reg['s']-1)
      if reg.w > -179:
        grass.run_command('g.region', w=reg['w']-1)
      else:
        grass.run_command('g.region', w=-180)
      if reg.e < 179:
        grass.run_command('g.region', e=reg['e']+1)
      else:
        grass.run_command('g.region', e=180)
      # Then resample and run r.neighbors as a "primer"
      grass.run_command('g.copy', rast=self.ice_raw_import[i] + ',' + self.ice[i], overwrite=True)
      # Subsequent iterations with nearest-neighbor interpolation until we get past the size of the topography
      res = 1.
      next_res = res / 2.
      while next_res > reg['nsres']: # Up to the semi-final needed resolution
        res /= 2.
        next_res /= 2.
        grass.run_command('g.region', res=res) # Climate model native resolution
        print "resolution:", grass.region()['nsres']
        grass.run_command('r.resample', input=self.ice[i], output=self.ice[i], overwrite=True)
        grass.run_command('r.neighbors', input=self.ice[i], output=self.ice[i], size=7, overwrite=True)
      # Final resolution
      grass.run_command('g.region', region='default')
      grass.run_command('r.resample', input=self.ice[i], output=self.ice[i], overwrite=True)
      
def runoff_input_meteoric(self):
  """
  Input water balance, multiplied by area, for runoff model.
  This one is simply the meteoric balance.
  """
  # m/yr to m**3/s
  for i in range(len(self.ages)):
    print self.runoff_input_meteoric[i]
    mcstr = self.runoff_input_meteoric[i]+" = "+self.wb[i]+" * cellsize_meters2 / 31557600."
    mapcalc(mcstr, overwrite=True)
    print ''

def dVi_dt(self, i):
  """
  (ice_mass_balance_midpoint, but now as a function instead of a grid-builder)
  Utility to compute an equation for the on-midpoint without saving a grid
  (and therefore saving time and HDD space)
  This gives dH_ice [m] * cellsize [m**2] / dt [yr]
  """
  mcstr_midpoints = '( (' + self.ice[i+1] + ' - ' + self.ice[i] + ') / ' + self.dt[i] + ') * cellsize_meters2'
  return mcstr_midpoints

def runoff_input_ice(self):
  """
  Weight ice mass balance at a time-step based on its position next to
  the midpoints of the mass-balance differences (i.e this produces
  numerical diffusion)
  """
  # might make more sense for 1, len(), but grandfathered in
  for i in range(len(self.ages)-2):
    agenow = self.ages_numeric[i+1]
    print self.runoff_input_ice[i]
    print "Ice mass balance centered around", agenow, "ka"
    younger = float(self.midpoint_age[i+1])
    older = float(self.midpoint_age[i])
    iweight = (agenow - younger) / (older - younger) # weight on older
    iplusweight = (older - agenow) / (older - younger) # weight on younger
    print iweight, iplusweight
    # Change in ice volume per time
    mcstr_dVi_dt_ts = '( (' + self.climate.dVi_dt(self, i) + ' * ' + str(iweight) + ') + (' + self.climate.dVi_dt(self, i+1) + ' * ' + str(iplusweight) + ') )'
    # Converted to runoff: take the negative (loss of ice is gain of runoff),
    # multiply by ice density / water density, and convert from m**3/yr to m**3/s
    mcstr_runoff_input_ice_conversion = " * -1 * "+str(self.rho_ice/self.rho_water)+" / "+str(self.seconds_in_year)
    # Put it all together
    mcstr = self.runoff_input_ice[i]+" = "+ mcstr_dVi_dt_ts + mcstr_runoff_input_ice_conversion
    print mcstr
    mapcalc( mcstr , overwrite=True)
    print ''
      
  #i = len(self.ages)-2
  # 0 a and the next time before it should both have ~0 ice loss, or at least be about the same
  print self.runoff_input_ice[i+1]+' = '+self.runoff_input_ice[i]
  grass.run_command('g.copy', rast=self.runoff_input_ice[i] + ',' + self.runoff_input_ice[i+1], overwrite=True)

  # No value for 21000 - at least, with that being my first time step (ts)

def runoff_input_total(self):
  # Sum meltwater and meteoric runoff to get full runoff grid
  # Skip first "ages" b/c not doing this for 21 ka
  for i in range(1,len(self.ages)):
    # m * km**2 / yr
    print self.runoff_input[i-1]
    mcstring = self.runoff_input[i-1] + " = " + self.runoff_input_meteoric[i] + " + " + self.runoff_input_ice[i-1]
    print mcstring
    mapcalc( mcstring , overwrite=True)
    print ''


