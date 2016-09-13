from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import matplotlib.pyplot as plt
from grass import script as grass
from grass.script import array as garray
import re
#import matplotlib.mpl as mpl
#import matplotlib.scale as scale
from matplotlib.colors import Normalize

# run in global database 30as

imshow = True

m = Basemap(width=7000000,height=6500000,
            resolution='l',projection='laea',\
            lat_ts=52,lat_0=52,lon_0=-100.)

grass.run_command('g.region', rast='wb_000000')
nx = grass.region()['cols']
ny = grass.region()['rows']

# Now the actual data
WB_000000 = garray.array()
WB_000000.read("wb_000000", null=np.nan)

WB = garray.array()
WB.read("wb_021000", null=np.nan)

WBdiff = np.flipud(WB - WB_000000) * 1000 # mm/yr


# Colorbar is bipolar:
# http://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib

class myNorm(Normalize):    
  def __init__(self,linthresh,vmin=None,vmax=None,clip=False):
    Normalize.__init__(self,vmin,vmax,clip)
    self.linthresh=linthresh
    self.vmin, self.vmax = vmin, vmax
    
  def __call__(self, value, clip=None):
    if clip is None:
      clip = self.clip

    result, is_scalar = self.process_value(value)

    self.autoscale_None(result)
    vmin, vmax = self.vmin, self.vmax
    if vmin > 0:
      raise ValueError("minvalue must be less than 0")
    if vmax < 0:
      raise ValueError("maxvalue must be more than 0")            
    elif vmin == vmax:
      result.fill(0) # Or should it be all masked? Or 0.5?
    else:
      vmin = float(vmin)
      vmax = float(vmax)
      if clip:
        mask = ma.getmask(result)
        result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                          mask=mask)
      # ma division is very slow; we can take a shortcut
      resdat = result.data

      resdat[resdat>0] /= vmax
      resdat[resdat<0] /= -vmin
      resdat=resdat/2.+0.5
      result = np.ma.array(resdat, mask=result.mask, copy=False)

    if is_scalar:
        result = result[0]

    return result

  def inverse(self, value):
    if not self.scaled():
      raise ValueError("Not invertible until scaled")
    vmin, vmax = self.vmin, self.vmax

    if cbook.iterable(value):
      val = ma.asarray(value)
      val=2*(val-0.5) 
      val[val>0]*=vmax
      val[val<0]*=-vmin
      return val
    else:
      if val<0.5: 
        return  2*val*(-vmin)
      else:
        return val*vmax

  def makeTickLabels(self, nlabels):
    #proportion_min = np.abs(self.vmin - self.linthresh) / ( np.abs(self.vmin - self.linthresh) + np.abs(self.vmax - self.linthresh) )
    #nlabels_min = np.round((nlabels - 1) * proportion_min) # save the last label for the midpoint
    #nlabels_max = nlabels - 1 - nlabels_min
    # Will always add a point at the middle
    ticks = np.concatenate(( np.linspace(0, 0.5, nlabels/2+1), np.linspace(.5, 1, nlabels/2+1)[1:]))
    tick_labels = np.concatenate(( np.linspace(self.vmin, self.linthresh, nlabels/2 + 1), np.linspace(self.linthresh, self.vmax, nlabels/2 + 1)[1:] ))
    tick_labels = list(tick_labels)
    for i in range(len(tick_labels)):
      tick_labels[i] = '%.2f' %tick_labels[i]
    return ticks, tick_labels
  
  
# imshow
maps = [WBdiff]
colormaps = [cm.GMT_polar_r, cm.GMT_polar, cm.GMT_polar]
map_titles = ['P$-$ET (21 ka) $-$ P$-$ET (0 ka)']
units = ['Precip. minus Evapotransp. difference [mm yr$^{-1}$]']
#fig, axes = plt.subplots(1, 3, figsize=(24,6))
fig = plt.figure(figsize=(8,6))
for i in range(len(maps)):
  ax = fig.add_subplot(len(maps), 1, i+1)
  lons = np.linspace(-180+.125, 180-.125, 1440)
  lats = np.linspace(-90+.125, 90-.125, 720)
  transformed = m.transform_scalar(maps[i],lons,lats,nx,ny)
  maxvalue = np.max(transformed)
  minvalue = np.min(transformed)
  mn = myNorm(0, minvalue, maxvalue)
  transformed_scaled = mn(transformed)
  ticks, labels = mn.makeTickLabels(11)
  labels = np.round(np.array(labels).astype(float)).astype(int)
  #ax = axes[i]
  im = m.imshow(transformed_scaled, colormaps[i], interpolation='nearest')#, vmin=-1*absmaxvalue, vmax=absmaxvalue)
  m.drawcoastlines(color='k')
  m.drawcountries(color='k')
  m.drawstates(color='k')
  cbar = m.colorbar(im, ticks=ticks)
  cbar.ax.set_yticklabels(labels)
  cbar.set_label(units[i], fontsize=16)
  ax.set_title(map_titles[i], fontsize=24)

fig.tight_layout()
plt.savefig('../figures/P_ET_T.pdf')
#plt.savefig('../figures/P_ET_T.png')
#plt.show()
plt.close()

"""
Tdiff_transformed = m.transform_scalar(Tdiff,lons,lats,nx,ny)
absmaxvalue = np.max(np.abs(Tdiff_transformed))
maxvalue = np.max(Tdiff_transformed)
minvalue = np.min(Tdiff_transformed)

mn = myNorm(0, minvalue, maxvalue)

Tdiff_transformed_scaled = mn(Tdiff_transformed)
ticks, labels = mn.makeTickLabels(11)

fig, ax = plt.subplots(figsize=(8,6))
im = m.imshow(Tdiff_transformed_scaled, cm.GMT_polar, interpolation='nearest')#, vmin=-1*absmaxvalue, vmax=absmaxvalue)

m.drawcoastlines(color='k')
m.drawcountries(color='k')
m.drawstates(color='k')

#cbar = fig.colorbar(im, ticks=ticks, orientation='horizontal')
#cbar.ax.set_xticklabels(labels)
cbar = m.colorbar(im, ticks=ticks)
cbar.ax.set_yticklabels(labels)
cbar.set_label(r'Temperature difference [$^\circ$C]', fontsize=16)

fig.tight_layout()

plt.show()
"""

