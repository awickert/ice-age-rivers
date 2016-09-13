#! /usr/bin/python

import os
os.chdir('/home/awickert/Dropbox/Papers/InProgress/DrainageMethods/GRASSplot_active_editing')

# Limit plotted rivers to those that are being shown in the plots
wherestr = "river = 'Mississippi' OR river = 'Colorado' OR river = 'Rio Grande' OR river = 'Susquehanna' OR river = 'Hudson' OR river = 'Saint Lawrence' OR river = 'Hudson Strait' OR river = 'Mackenzie' OR river = 'Columbia'"

#import grassplot2 as gp
import cartopy
from matplotlib import pyplot as plt
import numpy as np
from grass import script as grass
from matplotlib.colors import Normalize
from matplotlib import colors
import matplotlib
from grass.script import array as garray
from matplotlib import cm
from matplotlib.colors import Normalize, LinearSegmentedColormap
import re
from scipy.interpolate import interp1d
from grass.pygrass import gis

grass.run_command('g.region', rast='ice_000000')

def read_vector_lines(vect, area=True, cats=None, wherestr=None):
  print vect
  # First the data
  # Parse the vertices from v.out.ascii
  all_lines_output = []
  if wherestr:
    grass.run_command('v.extract', input=vect, where=wherestr, output='tmp', overwrite=True)
    vertices_raw = grass.read_command('v.out.ascii', input='tmp', output='-', type='line,boundary', format='wkt', cats=cats)
  else:
    vertices_raw = grass.read_command('v.out.ascii', input=vect, output='-', type='line,boundary', format='wkt', cats=cats)
  vector_lines = vertices_raw.split('\n')
  for vector_line in vector_lines:
    if vector_line != '': # Last line should be empty, this will remove it safely
      vertices_output = []
      # strips parentheses and text, and then separates out coordiante pairs
      vertex_list = re.sub("[A-Z]|\(|\)", "", vector_line).split(', ')
      # Turn coordiante pairs into a numpy array and add to the output list
      all_lines_output.append( np.array([vertex.split() for vertex in vertex_list]).astype(float) )
  # And then the other attributes to go along with them
  return all_lines_output

def get_time_steps():
  indexmaps = grass.parse_command('g.list', type='vect', patt='drainage_basins_??????').keys()
  indexmaps = sorted(indexmaps)[::-1]
  ages = []
  for indexmap in indexmaps:
    ages.append( str(indexmap.split('_')[-1]) )
  #ages = ages[-25:] # TEMPORARY HACK-EY TO FINISH G12
  return ages

# From http://matplotlib.org/users/colormapnorms.html
class MidpointNormalize(colors.Normalize):
  def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
    self.midpoint = midpoint
    colors.Normalize.__init__(self, vmin, vmax, clip)

  def __call__(self, value, clip=None):
    # I'm ignoring masked values and all kinds of edge cases to make a
    # simple example...
    x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
    return np.ma.masked_array(np.interp(value, x, y))

etopo2 = np.genfromtxt('GRASScolors/etopo2', skip_footer=1)
z = etopo2[:,0].astype(int)
etopo2 = etopo2[(z > -6000) * (z < 6000)]
z = etopo2[:,0].astype(int)
r = etopo2[:,1].astype(float)
g = etopo2[:,2].astype(float)
b = etopo2[:,3].astype(float)
from scipy.interpolate import interp1d
ri = interp1d(z, r)
gi = interp1d(z, g)
bi = interp1d(z, b)
low_elev = np.min(z)
high_elev = np.max(z)
znew = np.linspace(low_elev, high_elev, 512)
znew = np.concatenate(( znew[znew<-1], [-1, 0], znew[znew>0])) # make sure key SL transition is intact!
rnew = ri(znew)
gnew = gi(znew)
bnew = bi(znew)
clscaled = np.linspace(0, 1, len(znew))
cdr = []
cdg = []
cdb = []
for i in range(len(znew)):
  cdr.append([clscaled[i], rnew[i]/255., rnew[i]/255.])
  cdg.append([clscaled[i], gnew[i]/255., gnew[i]/255.])
  cdb.append([clscaled[i], bnew[i]/255., bnew[i]/255.])
cdict = {'red': cdr, 'green': cdg, 'blue': cdb}
cm_etopo2 = LinearSegmentedColormap('etopo2',cdict,4096)

def midpoints(invar):
  return (invar[1:] + invar[:-1]) / 2

grass.run_command('g.region', res='.1')
toponow = garray.array()
toponow.read('toponow')

ages = get_time_steps()

for age in ages:

  grass.run_command('g.region', res='.1')
  reg = grass.region()
  s = reg['s']
  n = reg['n']
  w = reg['w']
  e = reg['e']
  nlats = reg['rows']
  nlons = reg['cols']

  lats = midpoints( np.linspace(s, n, nlats+1) )
  lons = midpoints( np.linspace(w, e, nlons+1) )

  z = garray.array()
  z.read('topo_'+age)

  grass.run_command('g.region', res='1')
  reg = grass.region()
  nlats = reg['rows']
  nlons = reg['cols']
  icelats = midpoints( np.linspace(s, n, nlats+1) )
  icelons = midpoints( np.linspace(w, e, nlons+1) )
  ice = garray.array()
  ice.read('ice_'+age)
  ice = ice[::-1]
  #ice[ice < 50] = np.nan

  drainage_basins = read_vector_lines('drainage_basins_'+age, wherestr=wherestr)
  rivers_1000cumec = read_vector_lines('rivers_1000cumec_'+age, area=False)
  shore = read_vector_lines('notocean_only_'+age)

  # EVENTUALLY UPDATE THIS TO BE FOR BASEMAP IN TOPO, GIA
  for basemap in ['topo']: #, 'GIA']:
    print ""
    print "***"
    print age
    print "***"
    print ""
    figsize = (8.4, 7.8)
    fig = plt.figure(figsize=figsize)
    ax = plt.axes( projection=cartopy.crs.LambertConformal(central_longitude=-100., central_latitude=52.0, standard_parallels = (44, 60)) )
    if basemap == 'topo':
      outpath = '/home/awickert/Dropbox/Papers/InProgress/DrainageMethods/figures/maps_alltimes_topo_20160808/map_'+gis.Location().name+'_'+age+'.png'
    elif basemap == 'GIA':
      outpath = '/home/awickert/Dropbox/Papers/InProgress/DrainageMethods/figures/maps_alltimes_gia_20160808/map_'+gis.Location().name+'_'+age+'.png'
    if os.path.exists(outpath):
      pass
    else:
      if basemap == 'topo':
        elev = ax.contourf(lons, lats, z[::-1], 150, cmap=cm_etopo2, transform=cartopy.crs.PlateCarree())
        elev.set_clim(vmin=low_elev, vmax=high_elev)
        ax.contourf(icelons, icelats, ice, levels=[40, 8000], transform=cartopy.crs.PlateCarree(), colors='w', alpha=.8)
      elif basemap == 'GIA':
        giamap = ax.contourf(lons, lats, z[::-1]-toponow[::-1], 20, cmap=cm.seismic, transform=cartopy.crs.PlateCarree(), norm=MidpointNormalize(midpoint=0.), zorder=0)
        #ax.add_patch(
        #    matplotlib.patches.Rectangle(
        #        (0.05, 0.05),   # (x,y)
        #        0.03,           # width
        #        0.5,            # height
        #        facecolor="white", zorder=1, transform=ax.transAxes))
        cbaxes = fig.add_axes([0.05, 0.05, 0.03, 0.5], zorder=2)
        cb = plt.colorbar(giamap, cax=cbaxes)  
        ax.contour(icelons, icelats, ice, levels=[40, 8000], linewidths=3, transform=cartopy.crs.PlateCarree(), colors='w', alpha=.8)
      ax.set_xlim((-3500000, 3500000))
      ax.set_ylim((-3500000, 3000000))
      ax.contour(icelons, icelats, ice, levels=[40], transform=cartopy.crs.PlateCarree(), colors='w', linewidth=4)
      #for line in border:
      #  ax.plot(line[:,0], line[:,1], linewidth=4, color='k', transform=cartopy.crs.PlateCarre())
      for line in shore:
        ax.plot(line[:,0], line[:,1], linewidth=.5, color='.2', transform=cartopy.crs.PlateCarree())
      for line in drainage_basins:
        ax.plot(line[:,0], line[:,1], color='0', linewidth=2, transform=cartopy.crs.PlateCarree())
      for line in rivers_1000cumec:
        ax.plot(line[:,0], line[:,1], color='blue', linewidth=1, transform=cartopy.crs.PlateCarree())
      plt.tight_layout()
      plt.subplots_adjust(left=0., right=1., top=1., bottom=0.)
      plt.savefig(outpath)
      plt.close()

