#! /usr/bin/python

import grassplot as gp
from mpl_toolkits.basemap import Basemap, cm
from matplotlib import pyplot as plt
import numpy as np
from grass import script as grass
from matplotlib.colors import Normalize
from matplotlib import colors
import matplotlib

m = Basemap(width=7000000,height=6500000,
            resolution='l',projection='laea',\
            lat_ts=52,lat_0=52,lon_0=-100.)

#axsize = (6.5,7)
axsize = (16, 10)

p = gp.grassplot(m)
etopo2low, etopo2high, etopo2cm = p.make_GRASS_etopo2_colormap()

files = sorted(grass.parse_command('g.list', type='rast', pattern='notocean*').keys())
ages = []
for f in files:
  ages.append(f[-6:]) # ICE5G, G12
  #ages.append(f[-8:]) # ICE3G
  #ages.append(f[-7:]) # ANU

# If needed:
#for age in ages:
#  print age
  # grass.run_command('r.shaded.relief', zmult=5, altitude=34, input=file, output='shaded_'+file[-6:], units='meters', overwrite=True) # ICE5G, G12
  #grass.run_command('r.shaded.relief', zmult=5, altitude=34, input='topo_'+age, output='shaded_'+age, overwrite=True) # ICE3G
#  grass.run_command('r.shaded.relief', zmult=5, altitude=34, input='topo_'+age, output='shaded_'+age, overwrite=False) # ANU

"""
for age in ages:
  grass.run_command('r.contour', input='ice_'+age, output='ice_50m_'+age, levels=50)
"""


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
  



toponow = p.rastprep('topo_'+ages[0])

plt.figure(figsize=axsize)
for age in ages[1:]:
  try:
    grass.run_command('g.region', rast='topo_000000')
    print "shaded relief"
    grass.run_command('r.relief', zscale=3, altitude=34, input='topo_'+age, output='shaded_'+age, overwrite=False) # ANU
  except:
    pass #probably already made
  #if age[-2] == '0':
  plt.clf()
  print "***", age, "***"
  shaded = p.rastprep('shaded_'+age)
  topo = p.rastprep('topo_'+age)
  ice_raw = p.rastprep('ice_'+age)
  tpi = topo + ice_raw
  #ice = ice_raw.copy()
  ##ice[ice_raw <= 40] = np.nan
  #ice[ice_raw <= 100] = np.nan # axis size
  drainage_basins = gp.read_vector_lines('drainage_basins_'+age)
  rivers_1000cumec = gp.read_vector_lines('rivers_1000cumec_'+age, area=False)
  shore = gp.read_vector_lines('notocean_only_'+age)
  ice = gp.read_vector_lines('ice_50m_'+age)
  i=0
  for line in shore:
    dist = ( (line[1:,0] - line[:-1,0])**2 + (line[1:,1] - line[:-1,1])**2 )**.5
    overshoot = (dist > 1.).nonzero()[0]
    for j in range(len(overshoot)):
      shore.append(line[:overshoot[j],:])
      line = line[overshoot[j]:,:]
      overshoot -= overshoot[j]
    shore[i] = line
    i+=1

  shademap = m.imshow(shaded, cmap=plt.cm.Greys)

  dtopo = topo - toponow
  """
  maxvalue = np.nanmax(dtopo)
  minvalue = np.nanmin(dtopo)
  mn = myNorm(0, minvalue, maxvalue)
  scaled = mn(dtopo)
  ticks, labels = mn.makeTickLabels(10)
  dtopomap = m.imshow(scaled, cmap=cm.GMT_polar_r, alpha=.6)
  cbar = m.colorbar(dtopomap, ticks=ticks)
  cbar.ax.set_yticklabels(labels)
  cbar.set_label('Elevation Difference [m]', fontsize=16)
  """

  #X, Y = np.meshgrid(p.lons[1:-1], p.lats[1:-1])
  #x, y = m(X, Y)
  #icemap = m.contour(x, y, ice)#, levels=[40.], colors='0')
  """
  icemargin = (ice < 50) * (ice > 30)
  icemargin = icemargin.astype(float)
  icemargin[icemargin == 0] = np.nan
  icecmap = colors.ListedColormap(['.9'])
  icemap = m.imshow(icemargin, cmap=icecmap)
  """
  
  

  colswitch = -np.nanmin(dtopo) / (np.nanmax(dtopo) - np.nanmin(dtopo))
  cdict = {'red': ((0.0, 0.0, 0.0),
                   (colswitch, 1.0, 1.0),
                   (1.0, 1.0, 1.0)),
           'green': ((0.0, 0.0, 0.0),
                     (colswitch, 1.0, 1.0),
                     (1.0, 0.0, 0.0)),
           'blue': ((0.0, 1.0, 1.0),
                    (colswitch, 1.0, 1.0),
                    (1.0, 0.0, 0.0))}
  my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,2048)
  dtopomap = m.imshow(dtopo, cmap=my_cmap, alpha=.6)
  cbar = plt.colorbar(ticks = np.linspace(np.nanmin(dtopo), np.nanmax(dtopo), 9))
  cbar.set_label('Elevation Difference [m]', fontsize=32)
  cbar.ax.tick_params(labelsize=24)

  for line in ice:
    m.plot(line[:,0], line[:,1], color='.75', linewidth=2, latlon=True)
  for line in shore:
    m.plot(line[1:-1,0], line[1:-1,1], linewidth=.5, color='.2', latlon=True)
  for line in drainage_basins:
    m.plot(line[:,0], line[:,1], color='0', linewidth=1, latlon=True)
  for line in rivers_1000cumec:
    m.plot(line[:,0], line[:,1], color='blue', linewidth=1, latlon=True)

  plt.title(str(int(age)/1000.)+' ka', fontsize=48)
  # G12
  #plt.savefig(str(int(age)/1000.)+' ka', fontsize=48)
  # 5G
  #plt.savefig('/home/awickert/Dropbox/Papers/InProgress/DrainageMethods/figures/topodiff_mapsICE5G/drainage_'+age+'.png')
  #plt.savefig('/home/awickert/Dropbox/Papers/InProgress/DrainageMethods/figures/topodiff_mapsICE5G/drainage_'+age+'.pdf')
  # 3G
  plt.savefig('/home/awickert/Dropbox/Papers/InProgress/DrainageMethods/figures/topodiff_mapsICE3G/drainage_'+age+'.png')
  plt.savefig('/home/awickert/Dropbox/Papers/InProgress/DrainageMethods/figures/topodiff_mapsICE3G/drainage_'+age+'.pdf')
  # ANU
  #plt.savefig('/home/awickert/Dropbox/Papers/PaleohydrologyNorthAmerica/figures/topodiff_mapsANU/drainage_'+age+'.png')
  # 6G
  #plt.savefig('/home/awickert/Dropbox/Papers/InProgress/DrainageMethods/figures/topodiff_mapsICE6G/drainage_'+age+'.png')
  #plt.savefig('/home/awickert/Dropbox/Papers/InProgress/DrainageMethods/figures/topodiff_mapsICE6G/drainage_'+age+'.pdf')

  #plt.savefig('figures5G/drainage_'+age+'.svg')
  #plt.savefig('figures5G/drainage_'+age+'.pdf')

