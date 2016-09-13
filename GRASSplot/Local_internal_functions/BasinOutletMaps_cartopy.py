from glob import glob # to get a list of files of interest
import numpy as np
from matplotlib import pyplot as plt
import os
from mpl_toolkits.basemap import Basemap, cm
from grass import script as grass
import matplotlib.path as mpath
import matplotlib.patches as mpatches 
#import grassplot as gp
import cartopy
import cartopy.feature as cfeature
import re

wherestr = "river = 'Mississippi' OR river = 'Colorado' OR river = 'Rio Grande' OR river = 'Susquehanna' OR river = 'Hudson' OR river = 'Saint Lawrence' OR river = 'Hudson Strait' OR river = 'Mackenzie' OR river = 'Columbia'"

# Import files
lat_lon_name = np.loadtxt('site_list_lat_lon_name', delimiter=', ', dtype=str)
lat = lat_lon_name[:,0].astype(float)
lon = lat_lon_name[:,1].astype(float)
rivername = lat_lon_name[:,2]
rivername[1] = 'Hudson Strait\nand Hudson Bay'
rivername[2] = 'St. Lawrence'

def read_vector_lines(vect, area=True, cats=None, wherestr=None):
  print vect
  print wherestr
  all_lines_output = []
  # "Where" command not working!
  if wherestr:
    grass.run_command('v.extract', input=vect, where=wherestr, output='tmp', overwrite=True)
    vertices_raw = grass.read_command('v.out.ascii', input='tmp', output='-', type='line,boundary', format='wkt', cats=cats)
  else:
    vertices_raw = grass.read_command('v.out.ascii', input=vect, output='-', type='line,boundary', format='wkt', cats=cats)
  vector_lines = vertices_raw.split('\n')
  for vector_line in vector_lines:
    if vector_line != '': # Last line should be empty, this will remove it safely
      vertices_output = []
      vertex_list = re.sub("[A-Z]|\(|\)", "", vector_line).split(', ')
      all_lines_output.append( np.array([vertex.split() for vertex in vertex_list]).astype(float) )
  return all_lines_output

# River mouth regions - run this inside NA30asBristol
rmr = read_vector_lines('river_mouth_regions', wherestr=wherestr)

# Offsets for rivers
xoffset = [0, 2E5, 0, 0, 0, -2E5, 0, 1E5, -2E5]
yoffset = [3.2E5, -6.5E5, -2.5E5, 8E4, -3.5E4, 1.5E5, -2.5E5, 5E4, -2E5]
horizontalalignment = ['center', 'right', 'right', 'left', 'center', 'center', 'center', 'left', 'right']
rotation = [0, 0, 30, 0, -55, 0, 0, 0, 0]

fig = plt.figure(1, figsize=(8,7.43))
ax = plt.axes( projection=cartopy.crs.LambertConformal(central_longitude=-100., central_latitude=52.0, standard_parallels = (44, 60)) )
ax.coastlines(color='.7')
ax.add_feature(cartopy.feature.BORDERS, linestyle='-', color='.7')
ax.add_feature(cartopy.feature.LAKES, linestyle='-', color='.7')
#m.drawstates(color='.7')
#x, y = m(lon, lat)
ax.plot(lon, lat, 'kd', transform=cartopy.crs.PlateCarree())
#m.plot(x, y, 'ko', markersize=12, marker='d')
# river mouth regions
for line in rmr:
 ax.plot(line[:,0], line[:,1], color='.3', linewidth=4, transform=cartopy.crs.PlateCarree())
#for label, xpt, ypt, xoff, yoff, halign, rot in zip(rivername, x, y, xoffset, yoffset, horizontalalignment, rotation):
#  plt.text(xpt+xoff, ypt+yoff, label, fontsize=16, fontweight='bold', rotation=rot, horizontalalignment=halign)
ax.set_xlim((-3500000, 3500000))
ax.set_ylim((-3500000, 3000000))
plt.tight_layout()
plt.show()

#plt.savefig('../rivers_and_river_mouth_regions.svg')
#plt.savefig('../rivers_and_river_mouth_regions.pdf')
#plt.savefig('../rivers_and_river_mouth_regions.png')

#plt.show()


"""
  fig = plt.figure(1)
  ax = fig.add_subplot(1,1,1)
  verts = []
  for l in line:
    verts.append(tuple(l))
  codes = []
  for i in range(len(verts)):
    if i == 0:
      codes.append(mpath.Path.MOVETO)
    elif i == (len(verts) - 1):
      codes.append(mpath.Path.CLOSEPOLY)
    else:
      codes.append(mpath.Path.LINETO)
  path = mpath.Path(verts, codes)
  patch = mpatches.PathPatch(path, facecolor='blue', edgecolor='black')
  ax.add_patch(patch)

fig = plt.figure(1)
ax = fig.add_subplot(1,1,1)
patch = mpatches.PathPatch(path, facecolor='blue', edgecolor='black')
ax.add_patch(patch)
plt.xlim((-120, -50))
plt.ylim((50, 75))
plt.show()
"""


