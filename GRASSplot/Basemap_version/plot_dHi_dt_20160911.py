from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import matplotlib.pyplot as plt
from grass import script as grass
from grass.script import array as garray
import re
import matplotlib
import grassplot as gp
from matplotlib.colors import Normalize

imshow = True

m = Basemap(width=7000000,height=6500000,
            resolution='l',projection='laea',\
            lat_ts=52,lat_0=52,lon_0=-100.)
p = gp.grassplot(m)

# Now the actual data
axsize=(8,6)
dHi_dt = p.rastprep('dHi_dt_014250', figsize=axsize, resolution=150) 

# Use red--blue colormap for loss or gain -- red loss, blue gain
# WHOA - IS THIS JUST GMT_POLAR??!! -- ah, with myNorm
colswitch = -np.nanmin(dHi_dt) / (np.nanmax(dHi_dt) - np.nanmin(dHi_dt))
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
  

# Region of calculation box
def read_vector_lines(vect):
  # Parse the vertices from v.out.ascii
  all_lines_output = []
  vertices_raw = grass.read_command('v.out.ascii', input=vect, output='-', type='line,boundary', format='wkt')
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

# Region of calculation box
def read_vector_points(vect):
  # Parse the vertices from v.out.ascii
  all_points_output = []
  vertices_raw = grass.read_command('v.out.ascii', input=vect, output='-', type='point', format='wkt')
  vector_points = vertices_raw.split('\nPOINT')
  for vector_point in vector_points:
    if vector_point != '': # Last line should be empty, this will remove it safely
      vertices_output = []
      # strips parentheses and text, and then separates out coordiante pairs
      vertex_list = re.sub("[A-Z]|\(|\)", "", vector_point).split(', ')
      # Turn coordiante pairs into a numpy array and add to the output list
      all_points_output.append( np.array([vertex.split() for vertex in vertex_list][0]).astype(float) )
  # And then the other attributes to go along with them
  return np.array(all_points_output)

# used v.in.region and v.db.addtable <-- not sure if latter is n
grass.run_command('v.in.region', output='bounding_box', type='line', overwrite=True)
grass.run_command('v.to.points', input='bounding_box', output='bounding_box_points', dmax=.25, overwrite=True)
#bounding_box = read_vector_lines('bounding_box')[0]
bounding_box = read_vector_points('bounding_box_points')
# imshow
maps = [dHi_dt]
colormaps = [cm.GMT_polar_r]
map_titles = ['ICE-6G: 13.5 ka to 13.0 ka']
units = ['Rate of change in ice thickness [m yr$^{-1}$]']
fig = plt.figure(figsize=axsize)
for i in range(len(maps)):
  ax = fig.add_subplot(len(maps), 1, i+1)
  maxvalue = np.nanmax(dHi_dt)
  minvalue = np.nanmin(dHi_dt)
  mn = myNorm(0, minvalue, maxvalue)
  scaled = mn(dHi_dt)
  ticks, labels = mn.makeTickLabels(16)
  #ax = axes[i]
  im = m.imshow(scaled, colormaps[i], interpolation='nearest')#, vmin=-1*absmaxvalue, vmax=absmaxvalue)
  m.drawcoastlines(color='k')
  m.drawcountries(color='k')
  m.drawstates(color='k')
  x, y = m(bounding_box[:,0], bounding_box[:,1])
  m.plot(x, y, 'k', linewidth=2)
  cbar = m.colorbar(im, ticks=ticks)
  cbar.ax.set_yticklabels(labels)
  cbar.set_label(units[i], fontsize=16)
  ax.set_title(map_titles[i], fontsize=24)

fig.tight_layout()
plt.savefig('/home/awickert/Dropbox/Papers/InProgress/DrainageMethods/figures/dHi_dt.pdf')
#plt.show()


