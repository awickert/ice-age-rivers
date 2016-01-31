# Module to keep function namespaces separate but variable namespaces together
# Originally for lakes.py
# Then for drainage.py
# Started 15 NOV 2011 by ADW

# Import all of the built-in modules that "drainage" imports before its class 
# definition
from drainage import *

def start(self, n=85, s=15, w=-170, e=-40, res='0:0:30'):
  """
  Prints a starting message and sets the region.
  The region defaults to North America at 30-arcsecond resolution
  """
  print '####################################################################'
  print ' SETUP STEP, TO IMPORT AND CONFIGURE DATA AND MODELS '
  print '####################################################################'    
  print ''
  grass.run_command('g.region', n=n, s=s, w=w, e=e, res=res, flags='p')

def setConstants(self):
  # Constants
  self.rho_ice = 917.
  self.rho_water = 1000.
  self.seconds_in_year = 31556926.

def cellArea(self):
  """
  Create a map of cell areas in square meters
  Already done as part of standard initial import
  """  
  grass.mapcalc("area_meters2 = ( 111133. * nsres() ) * ( ewres() * (3.14159/180.) * 6371000. * cos(y()) )")

def adjacentAverage(self, array):
  """
  Utility function to find midpoints in a grid
  """
  return (array[1:] + array[:-1]) / 2.

def generateAges(self, LGM=False, ICE='ANU', custom=None):
  """
  1. Generate a list of ages from time-step file names
     - strings
     - numeric
     - ages on midpoints
     - dt between timesteps
  """
  print "Generating ages from time-step file-names."
  print ''
  # icemaps from old to young
  icemaps = sorted( grass.parse_command('g.list', type='raster', pattern='ice_??????').keys() )[::-1]
  # on time-steps
  self.ages = []
  self.ages_numeric = []
  for icemap in icemaps:
    # on time-step
    self.ages.append(re.findall('\d+', icemap)[0])
    # Use float instead of int to avoid undesired
    # floor division issues
    self.ages_numeric.append(float(self.ages[-1]))
  self.ages = np.array(self.ages)
  self.ages_numeric = np.array(self.ages_numeric)
  # on_midpoints
  self.midpoint_age_numeric = self.setup.adjacentAverage(self, self.ages_numeric)
  self.midpoint_age = []
  for age in self.midpoint_age_numeric:
    self.midpoint_age.append('%06d' %age)
  # dt
  self.dt_numeric = np.diff(self.ages_numeric[::-1])[::-1]
  self.dt = self.dt_numeric.astype(str) # no need for padding here

