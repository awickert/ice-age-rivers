from glob import glob # to get a list of files of interest
import numpy as np
from matplotlib import pyplot as plt
import os
from mpl_toolkits.basemap import Basemap, cm
from grass import script as grass
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import grassplot as gp

# v.extract in=river_mouth_regions out=river_mouth_regions_ESurf where="river = 'Mississippi' OR river = 'Rio Grande' OR river = 'Colorado' OR river = 'Columbia' OR river = 'Mackenzie' OR river = 'Hudson Strait' OR river = 'Saint Lawrence' OR river = 'Hudson' OR river = 'Susquehanna'"

# Import files
lat_lon_name = np.loadtxt('site_list_lat_lon_name', delimiter=', ', dtype=str)
lat = lat_lon_name[:,0].astype(float)
lon = lat_lon_name[:,1].astype(float)
rivername = lat_lon_name[:,2]
rivername[1] = 'Hudson Strait\nand Hudson Bay'
rivername[2] = 'St. Lawrence'

# River mouth regions - run this inside NA30asBristol
#rmr = gp.read_vector_lines('river_mouth_regions_qsr')
rmr = gp.read_vector_lines('river_mouth_regions_ESurf')

# Offsets for rivers
xoffset = [0, -2E5, -1E5, 1E5, -2E5, -2E5, 0, 1E5, -3E5]
yoffset = [3E5, -1E5, -5E4, 5E4, -2E5, 1E5, -3E5, 5E4, -3E5]
horizontalalignment = ['center', 'right', 'right', 'left', 'right', 'center', 'center', 'left', 'right']
rotation = [0, 30, 30, 0, 0, 0, 0, 0, 0]


# Location map
m = Basemap(width=7000000,height=6500000,
            resolution='l',projection='laea',\
            lat_ts=52,lat_0=52,lon_0=-100.)


# Offsets for rivers
xoffset = [0, 2E5, 0, 0, 0, -2E5, 0, 1E5, -2E5]
yoffset = [3.2E5, -6.5E5, -2.5E5, 8E4, -3.5E4, 1.5E5, -2.5E5, 5E4, -2E5]
horizontalalignment = ['center', 'right', 'right', 'left', 'center', 'center', 'center', 'left', 'right']
rotation = [0, 0, 30, 0, -55, 0, 0, 0, 0]

fig = plt.figure(1, figsize=(8,7.43))
ax = fig.add_subplot(1,1,1)
m.drawcoastlines(color='.7')
m.drawcountries(color='.7')
m.drawstates(color='.7')
x, y = m(lon, lat)
m.plot(x, y, 'ko', markersize=12, marker='d')
# river mouth regions
fig = plt.figure(1)
ax = fig.add_subplot(1,1,1)
for line in rmr:
  m.plot(line[:,0], line[:,1], color='.3', linewidth=4, latlon=True)
for label, xpt, ypt, xoff, yoff, halign, rot in zip(rivername, x, y, xoffset, yoffset, horizontalalignment, rotation):
  plt.text(xpt+xoff, ypt+yoff, label, fontsize=16, fontweight='bold', rotation=rot, horizontalalignment=halign)
plt.tight_layout()

plt.savefig('../rivers_and_river_mouth_regions.svg')
plt.savefig('../rivers_and_river_mouth_regions.pdf')
plt.savefig('../rivers_and_river_mouth_regions.png')

plt.show()


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


