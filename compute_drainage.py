#! /usr/bin/python

# Module to keep function namespaces separate but variable namespaces together
# For lakes.py
# Started 15 NOV 2011 by ADW

# Import all of the built-in modules that "lakes" imports before its class 
# definition
from drainage import *

def printstart(self):
  print '####################################################################'
  print ' UPDATE DRAINAGE STEP '
  print '####################################################################'    

# 3. Create flow routing grid: topo plus ice with ocean
def flow_routing_grid_withocean(self, subglacial=False):
  print "*************"
  print "Calculating grids of topo and ice to route surface water flow:"
  print "First pass includes routing across the ocean floor"
  if subglacial:
    print "*** Subglacial flow routing selected. ***"
  for age in self.ages:
    print age
    ice = 'ice_' + age
    topo = 'topo_' + age
    tpi = 'tpi_' + age # "topo plus ice"
    try:
      # Routing subglacially, hence the rho_i / rho_w
      # mapcalc(tpi + ' = ' + topo + ' + ' + str(ice_over_water) + ' * ' + ice)
      # Routing supraglacially - much of LIS is cold-based
      if subglacial:
        mapcalc(tpi + ' = ' + topo + ' + 0.917 * ' + ice, overwrite=True)
      else:
        mapcalc(tpi + ' = ' + topo + ' + ' + ice, overwrite=True)
    except:
      pass

# 2.5. Separate continent and ocean
# I am going for the shoreline approach, which works so long
# as every ocean touches the edge of my lat/lon map.
# Also, I should note that this works because there are no
# major lakes on the map margins or on the # 180 meridian
def separate_oceans(self):
  """
  Uses vector logic to separate 
  """
  print "*************"
  print "Setting ocean (defined by tpi 'shoreline') to NULL"

  #print "Changing region minutely to create edges for vector map" # N and S set to deal with mapcalc error
  #grass.run_command('g.region',flags='p',w='-179.9999',e='180',n='75',s='-68',rows='256',cols='512')

  for age in self.ages:
    print age
    tpi = 'tpi_' + age # "topo plus ice"
    sl_binary="sl_binary_" + age
    ocean_binary="ocean_binary_" + age
    notocean_only="notocean_only_" + age

    print "Creating binary above/below sea level map and converting to vector..."
    # Make raster
    try:
      mapcalc (sl_binary + ' = ' + tpi + ' > 0', overwrite=True)
    except:
      pass
    # Convert raster to vector; "value" column has 1 (land) and 0 (water)
    grass.run_command('r.to.vect' , input=sl_binary , output=sl_binary , type='area' , overwrite=True, quiet=True)
    #grass.run_command('v.to.rast' , input=sl_binary , output='sl_binary_echo' , use='attr', attrcolumn='value' , overwrite=True)

    print "Finding values at edges of vector map..."
    # Get rid of old layer, in case it is there and in the way
    grass.run_command("g.remove" , type='vector', name=ocean_binary, flags='f')
    # add categories for boundaries of the input vector map, in layer 2:
    grass.run_command('v.category' , input=sl_binary , output=ocean_binary , layer='2' , type='boundary' ,  option='add', quiet=True)
    # add a table with columns named "left" and "right" to layer 2 of the input
    # vector map:
    grass.run_command('v.db.addtable', map=ocean_binary, layer='2', columns="left integer,right integer", quiet=True)
    # find and upload categories of left and right areas:
    grass.run_command('v.to.db', map=ocean_binary, option='sides', columns='left,right', layer='2', quiet=True)
    # Add a column to Layer 1 ("notocean")
    grass.run_command('v.db.addcolumn', map=ocean_binary, layer='1', columns="notocean integer", quiet=True)
    # Set default value in this column to "1" (i.e. land or intracontinental below SL areas)
    grass.write_command('db.execute', stdin='UPDATE ' + ocean_binary + ' SET notocean=1', input='-', quiet=True)

    print "Finding polygons that touch edge; those that do are ocean..."
    # db_select must have changed: causing problems. So need to convert to list via numpy squeeze
    edgel = np.squeeze(np.array(db.db_select(table=ocean_binary, sql='SELECT right FROM ' + ocean_binary + '_2 WHERE left=-1')))
    if edgel.ndim == 0:
      edgel = np.array([str(edgel)])
    edger = np.squeeze(np.array(db.db_select(table=ocean_binary, sql='SELECT left FROM ' + ocean_binary + '_2 WHERE right=-1')))
    if edger.ndim == 0:
      edger = np.array([str(edger)])
    # Remove duplicates from edge.txt list (shaves off a few milliseconds :) )
    edge = list(set(list(edgel) + list(edger)))

    # Loop through lines in this file to set "notocean=0" for ocean areas
    for edgeval in edge: # edgeval = integer on line in "edge_sorted.txt"
      # ocean where value = water (0) and at edge
      grass.write_command('db.execute', stdin='UPDATE ' + ocean_binary + ' SET notocean=0 where cat=' + edgeval + ' AND value=0', input='-', quiet=True)

    print "Creating a binary raster of continent (1) and ocean (0)"
    grass.run_command('v.to.rast', input=ocean_binary, output=ocean_binary, use='attr', attrcolumn='notocean', layer='1', overwrite=True, quiet=True)
    
    print "Producing topography with ocean removed..."
    # NULL values for ocean
    # leave the binary raster intact
    grass.run_command('g.copy', rast=ocean_binary + ',' + notocean_only, overwrite=True)
    grass.run_command('r.null', map=notocean_only, setnull='0') # Turns ocean to null
    # Just in case the ocean binary is somehow lacking ocean -- change this to 0!
    # This is really just a patch for something that came up while doing this on a global scale
    for age in self.ages:
      grass.run_command('r.null', map='ocean_binary_'+age, null='0') # Turns ocean to null

    print "Converting notocean_only maps to vector..."

  for age in self.ages:
    print age
    notocean_only="notocean_only_" + age
    grass.run_command('r.to.vect', input=notocean_only, output=notocean_only, type_='area', overwrite=True)

    print "Done."

# 3. Create flow routing grid: topo plus ice without ocean
def flow_routing_grid(self, subglacial=False):
  print "*************"
  print "Calculating grids of topo and ice to route surface water flow"
  for age in self.ages:
    print age
    notocean_only="notocean_only_" + age
    ice = 'ice_' + age
    topo = 'topo_' + age
    tpin = 'tpin_' + age # "topo plus ice, no-ocean"
    try:
      # Routing subglacially, hence the rho_i / rho_w
      # mapcalc(tpin + ' = (' + topo + ' + 0.917 * ' + ice + ') * ' + notocean_only)
      # Routing supraglacially - much of LIS is cold-based
      if subglacial:
        mapcalc(tpin + ' = (' + topo + ' + 0.917 * ' + ice + ') * ' + notocean_only, overwrite=True)
      else:
        mapcalc(tpin + ' = (' + topo + ' + ' + ice + ') * ' + notocean_only, overwrite=True)
    except:
      pass
#    mapcalc(tpi + ' = ' + topo + ' + ' + str(ice_over_water) + ' * ' + ice)

def apply_etopo2_colormap(self):
  """
  Give a real topographic/bathymetric colormap to tpi and tpin for easier 
  viewing
  """
  print "*************"
  print "Applying etopo2 colormap to flow-routing grids"
  print ""
  for age in self.ages:
    tpi = 'tpi_' + age # "topo plus ice"
    tpin = 'tpin_' + age # "topo plus ice, no-ocean"
    grass.run_command('r.colors', map=tpi, color='etopo2')
    grass.run_command('r.colors', map=tpin, color='etopo2')

# 3.5: flow routing
def flow_routing_r_watershed(self):
  print "*************"
  print "Using r.watershed to build flow accumulation grid and drainage networks"
  for age in self.ages:
    print age
    tpin = 'tpin_' + age # "topo plus ice, no-ocean"
    runoff_input = 'runoff_input_' + age # runoff per cell
    accumulation = 'accumulation_' + age
    flowdir = 'flowdir_' + age
    # SFD flow accumulation
    grass.run_command('r.watershed', elevation=tpin, flow=runoff_input, accumulation=accumulation, drainage=flowdir, overwrite=True, flags="s")
  print ''
      
# 3.6 Need to replace NaN in flow accumulation grid with 0... or maybe values
# from runoff generation grid
def accum_nulls(self):
  for age in self.ages:
    print age
    accumulation = 'accumulation_' + age
    accumulation_with_nulls = 'accumulation_with_nulls_' + age
    grass.run_command('g.copy', rast=accumulation+','+accumulation_with_nulls, overwrite=True)
    grass.run_command('r.null', map=accumulation, null=0)
    # or maybe a mapcalc with accumulation + runoff * isnull(accumulation)
    # So probably rename accumulation to accumulation_with_nulls, and then 
    # mapcalc a new "accumulation"
    # 0 probably better: the nulls seem to border areas of negative and positive accumulation
    # Yes - this makes sense - NaN happens when there is no flow to accumulate

# 3.62: Flow accumulation from meltwater inputs
def flow_accum_ice(self):
  print "*************"
  print "Using r.watershed to build flow accumulation grid and drainage networks"
  for age in self.ages:
    print age
    notocean_only="notocean_only_" + age
    tpin = 'tpin_' + age # "topo plus ice, no-ocean"
    runoff_input_ice = 'runoff_input_ice_' + age # runoff per cell
    accumulation_ice = 'accumulation_ice_' + age
    grass.run_command('r.watershed', elevation=tpin, flow=runoff_input_ice, accumulation=accumulation_ice, flags="s", overwrite=True)

# 3.62: Flow accumulation from meteoric inputs
def flow_accum_meteoric(self):
  print "*************"
  print "Using r.watershed to build flow accumulation grid and drainage networks"
  for age in self.ages:
    print age
    notocean_only="notocean_only_" + age
    tpin = 'tpin_' + age # "topo plus meteoric, no-ocean"
    runoff_input_meteoric = 'runoff_input_meteoric_' + age # runoff per cell
    accumulation_meteoric = 'accumulation_meteoric_' + age
    grass.run_command('r.watershed', elevation=tpin, flow=runoff_input_meteoric, accumulation=accumulation_meteoric, flags="s", overwrite=True)

# 3.65 Now calculate rivers: everything with Q > 1000 cumecs
def big_rivers(self):
  for age in self.ages:
    print age
    accumulation = 'accumulation_' + age
    rivers_1000cumec = 'rivers_1000cumec_' + age
    rivers_1000cumec_r_watershed = 'rivers_1000cumec_r_watershed_' + age
    #grass.run_command('g.rename', rast=rivers_1000cumec+','+rivers_1000cumec_r_watershed, overwrite=False)
    # Abs because sometimes the big rivers go negative because of how r.watershed
    # works with the nulls. (It probably thinks there is offmap flow.
    # My negative runoff areas get nowhere near -1000 cumecs - 
    # seem to have problems joining up together, as might be expected
    mapcalc(rivers_1000cumec + ' = abs(' + accumulation + ') > 1000', overwrite=True)

"""
def tmp_multi(self):
  # Super experimental
  import multiprocessing as multi
   
  # Find number of workers that can be used on system. This variable could 
  # also be set manually.
  workers = multi.cpu_count()
   
  # Check if workers are already being used
  if workers is 1 and "WORKERS" in os.environ:
      workers = int(os.environ["WORKERS"])
  if workers < 1:
      workers = 1
   
  # Initialize process dictionary
  proc = {}

  nproc = 0
   
  # Loop over jobs
  for age in self.ages:
    rivers_1000cumec = 'rivers_1000cumec_' + age
    rivers_1000cumec_thinned = 'rivers_1000cumec_thinned_' + age
    # Stream to vector
    grass.run_command('r.null', map=rivers_1000cumec, setnull=0)
    # Insert job into dictinoary to keep track of it
    proc[i] = grass.start_command('r.slope.aspect',
                                  elevation='elev_' + str(i),
                                  slope='slope_' + str(i))
    # If the workers are used up, wait for all of them from the last group to
    # finish.
    if i % workers is 0:
      for j in range(workers):
        proc.[i - j].wait()
   
  # Make sure all workers are finished.
  for i in range(jobs):
    if proc[i].wait() is not 0:
      grass.fatal(_('Problem running analysis on evel_' + str(i) + '.')
"""

# 3.7. Vectorize streams
def vectorize_streams(self):
  print "*************"
  print "Thinning and vectorizing stream networks"
  for age in self.ages:
    print age
    rivers_1000cumec = 'rivers_1000cumec_' + age
    rivers_1000cumec_thinned = 'rivers_1000cumec_thinned_' + age
    # Stream to vector
    grass.run_command('r.null', map=rivers_1000cumec, setnull=0)
    grass.run_command('r.thin', input=rivers_1000cumec, output=rivers_1000cumec_thinned, overwrite=True)
    grass.run_command('r.to.vect', input=rivers_1000cumec_thinned, output='tmp', type='line', overwrite=True)
    grass.run_command('v.clean', input='tmp', layer=1, output=rivers_1000cumec, tool='rmdangle', thresh='0.01,0.01', overwrite=True) # Removed "snap" - seems a bad idea and was causing problems
    print ''


# 3.75 Grow (buffer incl. old stuff) ocean for shore
def grow_ocean(self):
  for age in self.ages:
    print age
    ocean_binary="ocean_binary_" + age
    ocean_only="ocean_only_" + age
    ocean_plus_shore="ocean_plus_shore_" + age
    grass.run_command('g.copy', rast=ocean_binary + ',' + ocean_only, overwrite=True)
    grass.run_command('r.null', map=ocean_only, setnull='1') # Turns continent to null
    grass.run_command('r.grow', input=ocean_only, output=ocean_plus_shore, radius=1.5, overwrite=True) # Radius big enough to catch angled approaches (sqrt(2) = 1.41)

# And vectorize:
def vectorize_ocean_plus_shore(self):
  for age in self.ages:
    print age
    ocean_plus_shore="ocean_plus_shore_" + age
    grass.run_command('r.to.vect', input=ocean_plus_shore, output=ocean_plus_shore, type='area', overwrite=True)


# 3.8 Obtain locations of mouths of streams
def mouths(self):
  print "*************"
  print "Extracting river mouths"
  for age in self.ages:
    print age
    rivers_1000cumec = 'rivers_1000cumec_' + age
    stream_nodes = 'stream_nodes_' + age
    river_mouths = 'river_mouths_' + age
    ocean_plus_shore="ocean_plus_shore_" + age
    grass.run_command('v.to.points', input=rivers_1000cumec, output='tmp', use='node', flags='t', type='line', overwrite=True)
    grass.run_command('v.category', input='tmp', output='tmp2', opt='del', cat=-1, layer=2, overwrite=True)
    grass.run_command('v.category', input='tmp2', output='tmp', opt='del', cat=-1, layer=1, overwrite=True)
    grass.run_command('v.clean', input='tmp', output='tmp2', tool='rmdupl', overwrite=True)
    grass.run_command('v.category', input='tmp2', output=stream_nodes, opt='add', layer=1, overwrite=True)
    grass.run_command('v.db.addtable', map=stream_nodes, columns='shore DOUBLE PRECISION')
    grass.run_command('v.what.vect', map=stream_nodes, qmap=ocean_plus_shore, column='shore', qcolumn='value')
    grass.run_command('v.extract', input=stream_nodes, output=river_mouths, where='shore = 0', overwrite=True)    

# Now find everything in particular drainage basins
# (I made river_mouth_regions by hand)

def discharge_at_mouths(self):
# First, add river names to column
  for age in self.ages:
    print age
    river_mouths = 'river_mouths_' + age
    grass.run_command('v.db.addcolumn', map=river_mouths, columns='river VARCHAR')
    grass.run_command('v.what.vect', map=river_mouths, qmap='river_mouth_regions', column='river', qcolumn='river')

# Then add discharge - this may be negative
  for age in self.ages:
    print age
    river_mouths = 'river_mouths_' + age
    Q = 'accumulation_' + age
    grass.run_command('v.db.addcolumn', map=river_mouths, columns='Q DOUBLE PRECISION')
    grass.run_command('v.what.rast', map=river_mouths, rast=Q, column='Q') # MAY BE NEGATIVE
    # SO FIX IT!
    grass.write_command('db.execute', stdin='UPDATE ' + river_mouths + ' SET Q=abs(Q)', input='-')

# Then add the partitioned discharge
  for age in self.ages:
    print age
    river_mouths = 'river_mouths_' + age
    Q_ice = 'accumulation_ice_' + age
    Q_meteoric = 'accumulation_meteoric_' + age
    grass.run_command('v.db.addcolumn', map=river_mouths, columns='Q_ice DOUBLE PRECISION, Q_meteoric DOUBLE PRECISION')
    grass.run_command('v.what.rast', map=river_mouths, rast=Q_meteoric, column='Q_meteoric') # MAY BE NEGATIVE
    # SO FIX IT!
    grass.write_command('db.execute', stdin='UPDATE ' + river_mouths + ' SET Q_meteoric=abs(Q_meteoric)', input='-')
    # grass.run_command('v.what.rast', map=river_mouths, rast=Q_ice, column='Q_ice') # MAY BE NEGATIVE
    # Q_ice hugely underestimates - I am almost positive that this is becasue of the round to 0 issues in the GRASS code!
    # Calculate this instead based on the other columns.
    grass.write_command('db.execute', stdin='UPDATE ' + river_mouths + ' SET Q_ice = Q - Q_meteoric', input='-')

"""
# Then add lat and lon with v.to.db
  for age in self.ages:
    print age
    river_mouths = 'river_mouths_' + age
    grass.run_command('v.db.addcolumn', map=river_mouths, columns='lon DOUBLE PRECISION, lat DOUBLE PRECISION')
    grass.run_command('v.what.vect', map=river_mouths, qmap='river_mouth_regions', column='river', qcolumn='river')
"""

def discharge_hist(self, rivername):

  self.rivername = rivername # for manual plotting

  from matplotlib import pyplot as plt

  self.Q_t = [] # total
  self.Q_i = [] # ice melt
  self.Q_m = [] # meteoric

  for age in self.ages:
    Q_t = grass.parse_command('v.db.select', map='river_mouths_' + age, columns='Q', where="river='" + rivername + "'", flags="c")
    Q_t = list(Q_t)
    Q_i = grass.parse_command('v.db.select', map='river_mouths_' + age, columns='Q_ice', where="river='" + rivername + "'", flags="c")
    Q_i = list(Q_i)
    Q_m = grass.parse_command('v.db.select', map='river_mouths_' + age, columns='Q_meteoric', where="river='" + rivername + "'", flags="c")
    Q_m = list(Q_m)

    i=0
    for item in Q_t:
      Q_t[i] = float(item)
      i+=1

    i=0
    for item in Q_i:
      try:
        Q_i[i] = float(item)
      except:
        Q_i[i] = 0
      i+=1

    i=0
    for item in Q_m:
      try:
        Q_m[i] = float(item)
      except:
        Q_m[i] = 0
      i+=1

    Q_t = sum(Q_t) # 0 if no items - perfect!
    Q_i = sum(Q_i) # 0 if no items - perfect!
    Q_m = sum(Q_m) # 0 if no items - perfect!

    self.Q_t.append(Q_t)
    self.Q_i.append(Q_i)
    self.Q_m.append(Q_m)

  fig = plt.figure(1)
  ax = fig.add_subplot(111)
  ax.plot(self.ages_numeric, self.Q_t, 'k-', label="Total", linewidth=3)
  ax.plot(self.ages_numeric, self.Q_i, 'b-', label="Ice sheet", linewidth=2)
  ax.plot(self.ages_numeric, self.Q_m, 'g-', label="Meteoric", linewidth=2)
  ax.set_title(self.rivername, fontsize=16)
  ax.set_xlabel("Time [ka BP]", fontsize=16)
  ax.set_ylabel(r"Discharge [m$^3$/s]", fontsize=16)
  ax.legend(loc='upper left')
  plt.show()

# Build basin outlets
def build_basin_outlets(self):
  print "Building basin outlets."
  for age in self.ages:
    print age
    ocean_plus_shore = 'ocean_plus_shore_' + age
    inshore_line = 'inshore_line_' + age
    rivers_1000cumec = 'rivers_1000cumec_' + age
    rivers_1000cumec_brkint = 'rivers_1000cumec_brkint_' + age
    drainage_outlets = 'drainage_outlets_' + age
    # Create a line for overlay
    grass.run_command("v.category", input=ocean_plus_shore, output="tmplinecats", opt="add", type="line,boundary", overwrite=True, quiet=True)
    grass.run_command("v.extract", input="tmplinecats", output="tmpbdry", type="boundary", overwrite=True, quiet=True)
    grass.run_command("v.type", input="tmpbdry", output=inshore_line, from_type="boundary", to_type="line", overwrite=True, quiet=True)
    # Break intersections internal to rivers
    grass.run_command("v.clean", input=rivers_1000cumec, output=rivers_1000cumec_brkint, tool="break", overwrite=True, quiet=True)
    # Then use the v.clean trick to get river mouths
    grass.run_command("v.patch", input=rivers_1000cumec_brkint+","+inshore_line, output="tmp_patch", overwrite=True, quiet=True)
    grass.run_command("v.clean", input="tmp_patch", output="tmp", tool="break", error="tmp_crossings_nocats", overwrite=True, quiet=True)
    # Add database table and populate it with lon, lat, and river outlet names
    grass.run_command("v.category", input="tmp_crossings_nocats", output=drainage_outlets, option="add", overwrite=True, quiet=True)
    grass.run_command("v.db.addtable", map=drainage_outlets, columns="lon DOUBLE PRECISION, lat DOUBLE PRECISION, river VARCHAR", quiet=True)
    grass.run_command("v.to.db", map=drainage_outlets, option="coor", columns="lon,lat", units="degrees", quiet=True) 
    grass.run_command("v.what.vect", map=drainage_outlets, qmap="river_mouth_regions", column="river", qcolumn="river", quiet=True)

def check_for_duplicate_outlets(self):
  for age in self.ages:
    print age
    rivers_1000cumec = 'rivers_1000cumec_' + age
    drainage_outlets = 'drainage_outlets_' + age
    Q = 'accumulation_' + age
    grass.run_command("v.db.addcolumn", map=drainage_outlets, columns="rivercat INT", quiet=True)
    grass.run_command("v.db.addcolumn", map=drainage_outlets, columns="Q DOUBLE PRECISION", quiet=True)
    grass.run_command("v.what.vect", map=drainage_outlets, qmap=rivers_1000cumec, column="rivercat", qcolumn="cat", dmax=0.1, quiet=True)
    grass.run_command('v.what.rast', map=drainage_outlets, rast=Q, column='Q', quiet=True)

    # Create a line for overlay

# Defined a grid called "zeros": mapcalc(drainage_basins + ' = 0 * flowdir_000_0k')
# huh - shouldn't it be: mapcalc('zeros = 0 * topo_000_000k') # topo b/c that already exists
def build_basins_rast(self, river_name=None):
  # Before anything, define an array object
  #dbrast = garray.array()
  # First get list of all possible basins and assign numbers to them
  self.rivers = np.array(list(set(list(grass.parse_command('v.db.select', map='river_mouth_regions', columns='river', flags="c")))))
  self.rnum = np.arange(1,len(self.rivers)+1)
  grass.run_command('g.region', rast='ice_000000') # Just in case it was changed while running the code before
  # Then build basins
  for age in self.ages:
    print ""
    print "***", age, "***"
    print ""
    # Lat/lon
    drainage_outlets = 'drainage_outlets_' + age
    drainage_basins = 'drainage_basins_' + age
    flowdir = 'flowdir_' + age
    accumulation = 'accumulation_' + age
    # Ordering gets messed up, so need to "select" in 1 step
    dodata = list(grass.parse_command('v.db.select', map=drainage_outlets, where="river NOT NULL", columns='lon,lat,river,rivercat,Q', flags="c"))
    lon = []
    lat = []
    river = []
    rivercat = []
    Q = []
    for line in dodata:
      linesplit = line.split('|')
      lon.append(linesplit[0])
      lat.append(linesplit[1])
      river.append(linesplit[2])
      rivercat.append(linesplit[3])
      Q.append(linesplit[4])
    # Check for duplicates and zap them
    rc = np.array(rivercat)
    Qtmp = []
    indextmp = []
    for i in range(len(rivercat)):
      # if duplicates, append to new list
      if sum(rc == rivercat[i]) > 1:
        indextmp.append(rivercat[i])
    # And then condense that list to have one entry per category
    indextmp = list(set(indextmp))
    # Then for each of these, find the one with the largest Q and pop the rest
    # out of all of the lists
    Qnum = []
    for Qi in Q:
      Qnum.append(abs(float(Qi)))
    Qa = np.array(Qnum)
    for iti in indextmp:
      Qs = Qa[rc == iti]
      Qmax = np.max(Qs)
      Qout = Qs[Qs != Qmax]
      for qoi in Qout:
        try:
          ind = Q.index(str(qoi))
        except:
          try:
            ind = Q.index(str(-1*qoi))
          except:
            Q2 = np.array(Q, dtype='float64')
            Q3 = []
            for i in range(len(Q2)):
              Q3.append( round(Q2[i],4) )
            try:
              ind = Q3.index(round(qoi, 4))
            except:
              ind = Q3.index(round(-1*qoi, 4))
        # Remove them!
        lon.pop(ind)
        lat.pop(ind)
        river.pop(ind)
        rivercat.pop(ind)
        Q.pop(ind)
    # Build grid
    grass.mapcalc(drainage_basins+' = 0', overwrite=True) # Start with zeros
    reg = grass.region()
    res = reg['nsres'] # NS and EW the same
    dmax = res/2. # not really used anymore
    for i in range(len(lon)):
      print river[i]
      #if river[i] == 'Rio Grande':
      basin_ncells = 0 # start with this
      first_loop = True # flag for first time through while loop
      basinhere = False # Start this as false for first time through loop
      broken = False # Don't create a basin if you break out of the loop
      iter_counter = 0
      # Need to make sure that basin actually covers its area - otherwise 
      # march upstream until it does
      clon, clat = float(lon[i]), float(lat[i]) # Center lon and lat; they update w/ time
      lonc, latc = lon[i], lat[i]
      while (basin_ncells < 10000) and (basinhere == False) and (iter_counter < 10):
        # Might in the future want to write error log for too mnay iterations
        # But if it can't get above 10000 cells (~5000-10,000 km **2) then way too small
        # to produce 1000 m3/s of water (10 cm/s of rain on basin??)
        # Code without this has also let us march to hte headwaters.... and get stuck there
        if first_loop == False:
          # If still in loop at this point, clearly some problem! Need to step upstream.
          # (we start at a point between cells, which may be part of the problem)
          iter_counter += 1
          print "Marching upstream, iteration ", iter_counter, ":", river[i]
          grass.run_command('g.region', n=str(np.floor((clat+3.*res)/res)*res), s=str(np.ceil((clat-3.*res)/res)*res), e=str(np.ceil((clon+3.*res)/res)*res), w=str(np.floor((clon-3.*res)/res)*res), quiet=True)
          # 1.5 is too slow. Need to march several cells, at least for test
          #grass.run_command('g.region', n=str(float(lat[i])+1.5*res), s=str(float(lat[i])-1.5*res), e=str(float(lon[i])+1.5*res), w=str(float(lon[i])-1.5*res))
          # Placing point on center of next available upstream cell - find by accumulation
          grass.run_command('r.to.vect', input=accumulation, type='point', column='accum', out='tmp', overwrite=True, quiet=True)
          grass.run_command('v.db.addcolumn', map='tmp', columns="lon DOUBLE PRECISION, lat DOUBLE PRECISION, isriver INT", quiet=True)
          grass.run_command('v.to.db', map='tmp', option='coor', columns='lon,lat', units='degrees', quiet=True)
          #grass.run_command('v.what.vect', map='tmp', qmap='rivers_1000cumec_020_0k', column='isriver', qcolumn='value', dmax=dmax)
          # Need a dmax because of imperfect overlap
          # Now moving center point
          print "Moving pour point to channel cell center:"
          # Find the correct cell - SQL MAX() is hard to figure out in GRASS, so just do it in Python
          accum_str_list = list(grass.parse_command('v.db.select', map='tmp', flags='c', columns='accum'))#, where='isriver NOT NULL'))
          accum_list = []
          for item in accum_str_list:
            accum_list.append(abs(float(item))) # abs b/c of offmap negatives
          accum_array = np.array(accum_list)
          # Min should be nothing close to max - choose lowest (i.e. most upstream) that is > 1/2 of max
          # 1/2 should be good for picking the right side of trib junctions too, but if we pass those
          # there are additional problems...
          accum_selected = str(np.min(accum_array[accum_array > np.max(accum_array)/2]))
          print "Accum selected:", accum_selected
          # Now select the high cell out of those that have the river going through them
          llcent = list(grass.parse_command('v.db.select', flags='c', map='tmp', columns='lon,lat', where="accum="+accum_selected + " OR " + "accum=-"+accum_selected ))# where="isriver NOT NULL AND accum="+accum_selected))
          lonc, latc = llcent[0].split('|') # Got 'em!
          # And use c on the other side as an ill-advised var name for the float versions
          clon, clat = float(lonc), float(latc)
          # back to standard region
          grass.run_command('g.region', rast='ice_000000', quiet=True)
          # Now finally onto the basins themselves
        # Very first, to save more time, check if specific 
        if (river[i] == river_name) or (river_name == None):
          # First, check if there is already a basin here from a previous run
          # Hey, I can do better than below - direct r.what!
          basinqueery = grass.parse_command('r.what', map=drainage_basins, coordinates=lonc+','+latc)
          basinhere = bool(int(basinqueery.keys()[0].split('|')[-1]))
          if basinhere:
            print "Basin already exists at this location:", river[i]+':', latc+'N;', lonc+'E'
          else:
            # NO ADDITIONAL "ELSE" NEEDED: defaults to the chosen mid-grid pour point
            print "Computing basin:   ", river[i]+':', latc+'N;', lonc+'E'
            grass.run_command('r.water.outlet', input=flowdir, output='rwotmp', coordinates=lonc+','+latc, overwrite=True, quiet=True)
            bac = list(grass.parse_command('r.stats', flags='c', input='rwotmp', output='-'))
            # "basin" is "1"
            for line in bac:
              linesplit = line.split(' ')
              if linesplit[0] == '1':
                basin_ncells = int(linesplit[-1])
                print "Computed basin size in pixels:", basin_ncells
                break
            # At this point, should either go back through or continue on
            first_loop = False
          # Out of "while" loop here
        else:
          break # Break out of loop and go on to next pour point
                # if we don't have the proper river system
          broken=True
      # Sloppy double-if, but I think it will do the trick - nested statements
      # Yes! Now marching upstream properly - needs to be outside so the marching
      # happens, but still contained in the river_name if (which should probably
      # be further outside instead of being inside twice)
      if (basinhere == False) and (broken == False):
        if (river[i] == river_name) or (river_name == None):
          print "Found acceptable basin. Amalgamating it with the others."
          # Following command unimportant when using numpy arrays
          #start = time.time()
          grass.run_command('r.null', map='rwotmp', null=0) # Set all nulls to 0 so we don't create a map of nulls
          river_number = int(self.rnum[self.rivers==river[i]])
          
          mapcalc(drainage_basins +' = '+ drainage_basins +' + rwotmp * '+ str(river_number), overwrite=True) # Still allows you to overwrite itself - good.
          #print time.time() - start
          # NOT FASTER! Replace mapcalc with garray command to build array in memory - faster!
          #start = time.time()
          #tmp = dbrast[...].copy()
          #dbrast.read('rwotmp')
          #dbrast[...] += tmp
          #print time.time() - start
    # After loop, set all 0's back to null before vector map creation
    # Should move this into next step now
#  for age in self.ages:
#    print age
#    drainage_basins = 'drainage_basins_' + age
#    grass.run_command('r.null', map=drainage_basins, setnull=0)


def basins_to_null_int(self):
  # Here because on the 3G run, these rasters were not recognized as integers, and therefore there were category problems
  for age in self.ages:
    print age
    drainage_basins = 'drainage_basins_' + age
    grass.run_command('g.rename', rast=drainage_basins+',tmp', overwrite=True)
    mcstr = drainage_basins+' = int(tmp)'
    print mcstr
    mapcalc(mcstr)
    grass.run_command('r.null', map=drainage_basins, setnull=0)


def add_basins_rast(self):
  from grass.script import vector as vect
  # Just get highest flow accum cell in area if there is nothing otherwise
  # (keep Q=0 away).
  # Next 3 lines copied from below
  self.rivers = np.array(list(set(list(grass.parse_command('v.db.select', map='river_mouth_regions', columns='river', flags="c")))))
  self.rnum = np.arange(1,len(self.rivers)+1)
  # Prepare a new river mouth regions raster
  grass.run_command('v.db.addcolumn', map='river_mouth_regions', col='rnum double precision')
  for i in range(len(self.rivers)):
    grass.run_command('v.db.update', map='river_mouth_regions', column='rnum', value=self.rnum[i], where='river == "'+self.rivers[i]+'"')
  grass.run_command('v.to.rast', input='river_mouth_regions', output='river_mouth_regions_by_number', use='attr', attrcolumn='rnum', overwrite=True)
  for age in self.ages: #['008500', '006500', '005500', '005000', '003500', '003000', '002500']:
    #self.ages:#[21:23]: #self.ages[16:20]:#
    print ""
    print "***", age, "***"
    print ""
    drainage_outlets = 'drainage_outlets_' + age
    drainage_basins = 'drainage_basins_' + age
    flowdir = 'flowdir_' + age
    accumulation = 'accumulation_' + age
    # Get an array of rivers that have drainage basins provided already w/ the 1000 cumec cutoff
    drainage_basins = 'drainage_basins_' + age
    drainage_basins_cols = vect.vector_db_select(drainage_basins)['columns']
    drainage_basins_lol = vect.vector_db_select(drainage_basins)['values'].values()
    river_name_column_db = (np.array(drainage_basins_cols) == 'river').nonzero()[0][0]
    self.drainage_basins_rivernames = np.array(drainage_basins_lol)[:,river_name_column_db]
    first_round=True
    for rivername in self.rivers:
      if any(self.drainage_basins_rivernames == rivername):
        print rivername, "is all set."
      else:
        print rivername, "basin doesn't have a big enough river (1000 cumec in GCM + ICE-5G)."
        print "Picking the cell with highest discharge inside its mouth region's bounding box."
        if first_round:
          print "Changing nulls to 0 to amalgamate basins"
          grass.run_command('r.null', map=drainage_basins, null=0)
          first_round = False
        # Extract the portion of river_mouth_regions for this river
        grass.run_command('v.extract', input='river_mouth_regions', output='tmp', where="river='" + rivername + "'", overwrite=True, quiet=True)
        # Get full region for full resolution
        fullres = grass.region()['nsres'] # Assume NS and EW res are the same
        # Then find the largest discharge inside this bounding box
        grass.run_command('g.region', vect='tmp', quiet=True)
        # Change bounding box edges so cells are sure to fall in right place
        smreg = grass.region()
        fixreg = grass.region()
        mult = 1./fullres
        for bound in 'n','s','w','e':
          tmpbound = smreg[bound]
          #tmpbound -= np.floor(tmpbound)
          tmpbound *= mult # turns it into an integer problem w/ internal units
          tmpbound = np.round(tmpbound)
          tmpbound /= mult
          fixreg[bound] = tmpbound
          #print tmpbound, 'degrees'
        grass.run_command('g.region', n=fixreg['n'], s=fixreg['s'], w=fixreg['w'], e=fixreg['e'])
        # mapcalc( "tmp = (river_mouth_regions" == "+this_river_number+") * "+accumulation, overwrite=True)
        rmr_specific = garray.array()
        rmr_specific.read(accumulation)
        # Compute inside GRASS the region occupied by the river_mouth_regions vector
        rmr_more_specific = garray.array()
        rmr_more_specific.read('river_mouth_regions_by_number') # 1 in defined regions, 0 (or NULL) elsewhere (generated ouside of this script)
        not_other_basin = garray.array()
        not_other_basin.read(drainage_basins)
        not_other_basin[self.rnum[self.rivers == rivername]] = 0
        rmr_specific[...] *= (rmr_more_specific == self.rnum[self.rivers == rivername]) * (not_other_basin == 0)
        maxval = float(np.max(np.abs(rmr_specific)))
        if maxval > 0:
          rmr_specific[...] = 1. * (np.abs(rmr_specific) == maxval)
          rmr_specific.write('tmp_mouth', overwrite=True)
          # Obtain a point and its coordinates, from which drainage basin will be built
          grass.run_command('r.null', map='tmp_mouth', setnull=0, quiet=True)
          print "Creating point at mouth"
          grass.run_command('r.to.vect', input='tmp_mouth', output='tmp_mouth', type='point', overwrite=True, quiet=True)
          print "Setting region"
          grass.run_command('g.region', rast='ice_000000', quiet=True)
          print "Database..."
          grass.run_command("v.db.addcolumn", map='tmp_mouth', columns="lon DOUBLE PRECISION, lat DOUBLE PRECISION, cellval INTEGER")
          grass.run_command("v.to.db", map='tmp_mouth', option="coor", columns="lon,lat", units="degrees", quiet=True) 
          # Check that I haven't used this side method to build a basin for this before
          grass.run_command('v.what.rast', map='tmp_mouth', column='cellval', raster=drainage_basins, quiet=True)
          tmp_mouth_cols = vect.vector_db_select('tmp_mouth')['columns']
          tmp_mouth_lol = vect.vector_db_select('tmp_mouth')['values'].values()
          cvcol = (np.array(tmp_mouth_cols) == 'cellval').nonzero()[0][0]
          basin_status = tmp_mouth_lol[0][cvcol]
          #basin_already_there = bool(int(tmp_mouth_lol[0][cvcol]))
          # In case already changed to nulls (so won't pick up a 0):
          basin_already_there = True # Start with this, so always defined
          if basin_status == '' or basin_status == '0':
            basin_already_there = False
          if basin_already_there:
            print tmp_mouth_cols
            print tmp_mouth_lol
            print "You already have a basin in this location: either from a prior"
            print "run of this step of from something really weird... hope it's the"
            print "first!"
          else:
            # Get these coordinates into Python
            loncol = (np.array(tmp_mouth_cols) == 'lon').nonzero()[0][0]
            latcol = (np.array(tmp_mouth_cols) == 'lat').nonzero()[0][0]
            lon = tmp_mouth_lol[0][loncol]
            lat = tmp_mouth_lol[0][latcol]
            # And then build the drainage basin
            latc = lat
            lonc = lon
            print "Computing basin:   ", rivername+':', latc+'N;', lonc+'E'
            grass.run_command('r.water.outlet', input=flowdir, output='rwotmp', coordinates=lonc+','+latc, overwrite=True)
            bac = list(grass.parse_command('r.stats', flags='c', input='rwotmp', output='-'))
            # "basin" is "1"
            for line in bac:
              linesplit = line.split(' ')
              if linesplit[0] == '1':
                basin_ncells = int(linesplit[-1])
                print "Computed basin size in pixels:", basin_ncells
                break
            print "Found acceptable basin. Amalgamating it with the others."
            grass.run_command('r.null', map='rwotmp', null=0) # Set all nulls to 0 so we don't create a map of nulls
            river_number = int(self.rnum[self.rivers==rivername])
            mapcalc(drainage_basins +' = '+ drainage_basins +' + rwotmp * '+ str(river_number), overwrite=True) # Still allows you to overwrite itself - 
        else:
          print "No drainage in " + rivername + " at " + age
    if first_round:
      pass
    else:
      print "Returning NULL values to maps."
      grass.run_command('r.null', map=drainage_basins, setnull=0)

def build_basins_vect(self):
  self.rivers = np.array(list(set(list(grass.parse_command('v.db.select', map='river_mouth_regions', columns='river', flags="c")))))
  self.rnum = np.arange(1,len(self.rivers)+1)
  f = open('/tmp/rivers', 'w')
  # Build rules file
  for i in range(len(self.rnum)):
    f.write( str(self.rnum[i]) + '\t' + self.rivers[i] + '\n')
  f.close()
  for age in self.ages:
    print age
    drainage_basins = 'drainage_basins_' + age
    # Attach categories and convert to vector
    grass.run_command('r.category', map=drainage_basins, rules='/tmp/rivers')
    # Will create extra areas - so just use the "v" flag
    grass.run_command('r.to.vect', flags='v', input=drainage_basins, output=drainage_basins, type='area', overwrite=True)#, column='river_number')
    grass.run_command('v.db.renamecolumn', map=drainage_basins, column='label,river')


def rivers_by_area(self):
  for age in self.ages:
    print age
    tpin = 'tpin_' + age
    rivers_40000km2 = 'rivers_40000km2_' + age
    rivers_40000km2_thinned = 'rivers_40000km2_thinned_' + age
    accum_by_area = 'accum_by_area_' + age
    grass.run_command('r.watershed', flags='s', elev=tpin, accum=accum_by_area, stream=rivers_40000km2, thresh=40000, flow='cellsize_km2')
    grass.run_command('r.thin', input=rivers_40000km2, output=rivers_40000km2_thinned)
    grass.run_command('r.to.vect', input=rivers_40000km2_thinned, output=rivers_40000km2, type='line')

def basin_discharge(self):

  # Get list of river names
  # Actually, this just recreates "rivers"
  rmr_cols = vect.vector_db_select('river_mouth_regions')['columns']
  rmr_lol = vect.vector_db_select('river_mouth_regions')['values'].values()
  river_name_column = (np.array(rmr_cols) == 'river').nonzero()[0][0]
  self.river_names = list(np.array(rmr_lol)[:,river_name_column])
  # Include self.rivers - redundant but used (and rnum is needed)
  self.rivers = np.array(list(set(list(grass.parse_command('v.db.select', map='river_mouth_regions', columns='river', flags="c")))))
  self.rnum = np.arange(1, len(self.rivers)+1)

  # Rasterize the river mouth regions for later extraction of discharge data
  # grass.run_command('v.to.rast', 
  # Nah, just find largest in bounding box

  # Add columns to all of the basin areas
  for age in self.ages:#[1:4]:#[21:23]:
    print age
    drainage_basins = 'drainage_basins_' + age
    grass.run_command('v.db.addcolumn', map=drainage_basins, columns='Q_ice DOUBLE PRECISION, Q_meteoric_uncorrected DOUBLE PRECISION, Q_modern DOUBLE PRECISION, Q_meteoric DOUBLE PRECISION, Q_total DOUBLE PRECISION, Area_km2 DOUBLE PRECISION')

  # Mississippi is including the Atchafalaya--Red
  # Generate dictionary of modern discharges [cumecs]
  # Colorado and Rio Grande are pre-modern discharges
  # Rio Grande is a bit of a guess, based on:
  # "Observed historical discharge data from major rivers for climate model validation"
  # From the Max Planck Institute, 2000
  # Colorado is this, USGS, and following some reasoning for downstream rivers
  # Modern Hudson Strait is the Koksoak River
  # And note that I am assuming a similar bias in Q for the regions even as drainage area changes
  # Rivers with no value will just get calculated outputs
  # This column is the first-round normalization, but is not the same thing
  # that I am using now. In fact, it's not even the same thing that I am using
  # later -- and look for Q_meteoric_simulated_now.txt.
  Qmodern = dict(Mississippi=18434., Mackenzie=9910., Colorado=610., Rio_Grande=200., Hudson=606., Columbia=7500., Susquehanna=318., Saint_Lawrence=16800., Hudson_Strait=2800., Yukon=6428.)
  Qnormalization = dict(Mississippi=1., Mackenzie=1., Colorado=1., Rio_Grande=1., Hudson=1., Columbia=1., Susquehanna=1., Saint_Lawrence=1., Hudson_Strait=1., Yukon=1., Other_North=1., Barents=1., White_Sea=1., Kara_Northwest=1., Kara_Southeast=1., Daugava_and_Baltic=1., Vefsnfjord=1., Puget_Sound=1., Volga=1., Canadian_Cordillera=1., Newfoundland_and_Labrador=1., Skeena=1.) # Start with this but be sure to update every item!
  # Go from present backwards, so I can get the modern normalization in place 
  # first
  # No Q_t (Q_total) here yet because of the need to re-normalize discharge
  # Something I wrote here prevents rivers outside the domain from being included in the table...
  # ... but not sure what it is! Maybe one of the try/excepts.
  for age in self.ages[::-1]:#(self.ages[1:4] + [self.ages[-1]])[::-1]:
    print age
    drainage_basins = 'drainage_basins_' + age
    river_mouths = 'river_mouths_' + age
    accumulation = 'accumulation_' + age
    # Easiest to do first is to get the areas of all of the basins
    grass.run_command('v.to.db', map=drainage_basins, option='area', units='kilometers', columns="Area_km2")  
    # Insert modern discharge
    for rivername in self.river_names:
      rnu = rivername.replace(' ', '_') # Underscore for dict
      print rivername
      #if grass.parse_command('v.db.select', map=drainage_basins, column=
      try:
        modern_discharge = Qmodern[rnu]
        grass.run_command('v.db.update', map=drainage_basins, column='Q_modern', where="river="+'"'+rivername+'"', value=Qmodern[rnu])
      except:
        pass
      self.Q_i = [] # ice melt
      self.Q_m = [] # meteoric
      # Acquire and sum discharges
      Q_i = grass.parse_command('v.db.select', map=river_mouths, columns='Q_ice', where="river='" + rivername + "'", flags="c")
      Q_i = list(Q_i)
      Q_m = grass.parse_command('v.db.select', map=river_mouths, columns='Q_meteoric', where="river='" + rivername + "'", flags="c")
      Q_m = list(Q_m)

      i=0
      for item in Q_i:
        try:
          Q_i[i] = float(item)
        except:
          Q_i[i] = 0
        i+=1

      i=0
      for item in Q_m:
        try:
          Q_m[i] = float(item)
        except:
          Q_m[i] = 0
        i+=1
        
      if len(Q_m) == 0:
        print rivername, "too small; using its highest Q cell"
        # rnum = river number
        #this_river_number = str(self.rnum[(self.rivers)==rivername][0]) # not really needed anymore
        # Extract the portion of river_mouth_regions for this river
        grass.run_command('v.extract', input='river_mouth_regions', output='tmp', where="river='" + rivername + "'", overwrite=True, quiet=True)
        # Then find the largest discharge inside this bounding box
        grass.run_command('g.region', vect='tmp', quiet=True)
        # mapcalc( "tmp = (river_mouth_regions" == "+this_river_number+") * "+accumulation, overwrite=True)
        rmr_specific = garray.array()
        rmr_specific.read(accumulation)
        # Assuming that the tiny basins must be meteoric-only; not generally true, but almost always is
        Q_m.append(float(np.max(rmr_specific)))
        grass.run_command('g.region', rast='ice_000000', quiet=True)

      Q_i = sum(Q_i) # 0 if no items - perfect!
      Q_m = sum(Q_m) # 0 if no items - perfect!

      # Get ratio between Q_m and modern
      # This seems to be the cause for the sudden drop to the modern time-step!
      if age == '000_0k' or age == '000_00k' or age == '000_000k' or age == '000000':
        try:
          Qnormalization[rnu] = Qmodern[rnu] / Q_m
        except:
          # Keep the normalization value at 1: don't correct
          pass

      # Update proper columns with discharge values
      # Maybe I should have non-normalized values like these here
      grass.run_command('v.db.update', map=drainage_basins, column='Q_ice', where="river="+'"'+rivername+'"', value=Q_i)
      grass.run_command('v.db.update', map=drainage_basins, column='Q_meteoric_uncorrected', where="river="+'"'+rivername+'"', value=Q_m)
      #grass.run_command('v.db.update', map=drainage_basins, column='Q_meteoric', where="river="+'"'+rivername+'"', value=Q_m)
      #grass.run_command('v.db.update', map=drainage_basins, column='Q_total', where="river="+'"'+rivername+'"', value=Q_m + Q_i)
      grass.run_command('v.db.update', map=drainage_basins, column='Q_meteoric', where="river="+'"'+rivername+'"', value=Q_m * Qnormalization[rnu])
      grass.run_command('v.db.update', map=drainage_basins, column='Q_total', where="river="+'"'+rivername+'"', value=Q_m * Qnormalization[rnu] + Q_i)
      
      print ''

def basin_ice_volume(self):
  """
  Ice volumes in each basin - for cumulative sea level curves without needing 
  to average Q over adjacent time steps (i.e. avoiding numerical diffusion)
  """

  print ""
  print "***"  
  print "Initializing arrays for basin ice volume."
  print "***"  
  print ""

  # Get list of river names
  self.rivers = np.array(list(set(list(grass.parse_command('v.db.select', map='river_mouth_regions', columns='river', flags="c")))))
  self.rnum = np.arange(1,len(self.rivers)+1)

  grass.run_command("g.region", rast="cellsize_meters2") # make sure that region is full
  km_squared = garray.array()
  ice_thickness_array = garray.array()
  drainage_basins_array = garray.array()
  km_squared.read("cellsize_km2")  
  for age in self.ages:
    print ""
    print "***", age, "***"
    print ""
    ice = 'ice_' + age
    drainage_basins = 'drainage_basins_' + age
    grass.run_command('v.db.addcolumn', map=drainage_basins, columns='V_ice_km3 DOUBLE PRECISION')
    ice_thickness_array.read(ice)
    drainage_basins_array.read(drainage_basins)
    ice_volume = ice_thickness_array/1E3 * km_squared # km^3
    for river in self.rivers:
      rnu = river.replace(' ', '_') # Underscore for dict
      print river
      river_number = self.rnum[self.rivers == river]
      V_ice_km3_in_basin = np.sum(ice_volume * (drainage_basins_array == river_number)) # km^3
      grass.run_command('v.db.update', map=drainage_basins, column='V_ice_km3', where="river="+'"'+river+'"', value=V_ice_km3_in_basin)

def drainage_basins_AND(self):
  """
  Creates drainage basins at time step midpoints using an AND operator, mostly 
  to generate output meltwater discharges without numerical diffusion
  """
  # Get list of river names
  self.rivers = np.array(list(set(list(grass.parse_command('v.db.select', map='river_mouth_regions', columns='river', flags="c")))))
  self.rnum = np.arange(1,len(self.rivers)+1)
  for i in range(1, len(self.mean_age)):
    print ""
    print "***", self.mean_age[i], "***"
    print self.ages[i], self.ages[i+1]
    print ""
    drainage_basins_prev = 'drainage_basins_' + self.ages[i]
    drainage_basins_next = 'drainage_basins_' + self.ages[i+1]
    drainage_basins_AND_midpoint = 'drainage_basins_AND_' + self.mean_age[i]
    # Limit to self-overlaps and remove scraps
    grass.run_command('v.overlay', ainput=drainage_basins_prev, binput=drainage_basins_next, operator='and', output='tmp', overwrite=True, quiet=True)
    grass.run_command('v.db.addcolumn', map='tmp', columns='Area_km2', quiet=True)
    grass.run_command('v.to.db', map='tmp', option='area', units='kilometers', columns="Area_km2", quiet=True)  
    grass.run_command('v.extract', input='tmp', output=drainage_basins_AND_midpoint, where='a_cat = b_cat AND Area_km2 > 100', overwrite=True, quiet=True)
    # Clean up
    grass.run_command('v.db.dropcolumn', map=drainage_basins_AND_midpoint, columns='a_Q_ice,a_Q_meteoric_uncorrected,a_Q_modern,a_Q_meteoric,a_Q_total,a_Area_km2,a_V_ice_km3,b_cat,b_river,b_Q_ice,b_Q_meteoric_uncorrected,b_Q_modern,b_Q_meteoric,b_Q_total,b_Area_km2,b_V_ice_km3', quiet=True)
    grass.run_command('v.db.renamecolumn', map=drainage_basins_AND_midpoint, column='a_river,river', quiet=True)
    # Old categories
    grass.run_command('v.reclass', input=drainage_basins_AND_midpoint, output='tmp', column='a_cat', overwrite=True, quiet=True)
    grass.run_command('v.db.addtable', map='tmp', columns="river varchar", quiet=True)
    for i in range(len(self.rivers)):
      grass.run_command('v.db.update', map='tmp', column='river', value=self.rivers[i], where='cat='+str(i+1), quiet=True)
    grass.run_command('g.rename', vect='tmp,'+drainage_basins_AND_midpoint, overwrite=True, quiet=True)

def Q_i_unsmoothed(self):
  """
  Uses dVi_dt on time step midpoints to calculate discharge, avaeraging this 
  over the drainage basin areas before and after the time step midpoint, and 
  therefore smoothing in drainage basin area rather than in discharge, and 
  retaining peaks in overall meltwater discharge from the ice models
  """
  # Get list of river names
  self.rivers = np.array(list(set(list(grass.parse_command('v.db.select', map='river_mouth_regions', columns='river', flags="c")))))
  self.rnum = np.arange(1,len(self.rivers)+1)
  # Do calculations in memory in Python
  grass.run_command("g.region", rast="cellsize_meters2") # make sure that region is full
  km_squared = garray.array()
  dVi_dt_array = garray.array()
  drainage_basins_prev_array = garray.array()
  drainage_basins_next_array = garray.array()
  km_squared.read("cellsize_km2")  
  for i in range(3, len(self.mean_age)):
    print ""
    print "***", self.mean_age[i], "***"
    print self.ages[i], self.ages[i+1]
    print ""
    # Drainage basin and dVi_dt load
    drainage_basins_AND_midpoint = 'drainage_basins_AND_' + self.mean_age[i]
    drainage_basins_prev = 'drainage_basins_' + self.ages[i]
    drainage_basins_next = 'drainage_basins_' + self.ages[i+1]
    dVi_dt = self.dVi_dt[i]
    if i == 1:
      drainage_basins_prev_array.read(drainage_basins_prev)
    else:
      drainage_basins_prev_array = drainage_basins_next_array.copy()
    drainage_basins_next_array.read(drainage_basins_next)
    dVi_dt_array.read(dVi_dt)
    # Create a water discharge in cumecs from the ice inputs in km3 ice / yr
    Q_i_array = dVi_dt_array * -1. * 0.917 / 31556926.
    # Add a column to the table to accept these values
    grass.run_command('v.db.addcolumn', map=drainage_basins_AND_midpoint, columns='Q_ice DOUBLE PRECISION', quiet=True)
    for rnum in self.rnum:
      print "  "+self.rivers[rnum-1]
      Q_i_in_basin = float(np.sum( Q_i_array * ( 0.5 * (drainage_basins_prev_array == rnum) + 0.5 * (drainage_basins_next_array == rnum)) ) )
      grass.run_command('v.db.update', map=drainage_basins_AND_midpoint, column='Q_ice', value=Q_i_in_basin, where='cat='+str(rnum), quiet=True)

def write_basin_discharge_to_list(self, rivername='Mississippi'):
  from grass.script import db_select
  # Get list of river names
  self.rivers = np.array(list(set(list(grass.parse_command('v.db.select', map='river_mouth_regions', columns='river', flags="c")))))
  self.rnum = np.arange(1,len(self.rivers)+1)
  # Compile ages and discharges
  age = []
  Q_m = []
  Q_i = []
  for i in range(1, len(self.mean_age)):
    print ""
    print "***", self.mean_age[i], "***"
    print self.ages[i], self.ages[i+1]
    print ""
    river_mouths_prev = 'river_mouths_' + self.ages[i]
    river_mouths_next = 'river_mouths_' + self.ages[i+1]
    drainage_basins_AND_midpoint = 'drainage_basins_AND_' + self.mean_age[i]
    # Acquire and sum meteoric discharges
    # qm = db_select(sql = 'SELECT river,Q_meteoric FROM ',+river_mouths_prev+' where river <> ""') # For all
    qm_prev = db_select(sql = 'SELECT Q_meteoric FROM '+river_mouths_prev+' where river = "'+rivername+'"')
    qm_next = db_select(sql = 'SELECT Q_meteoric FROM '+river_mouths_next+' where river = "'+rivername+'"')
    qm_prev = np.sum(np.array(qm_prev, dtype=float))
    qm_next = np.sum(np.array(qm_next, dtype=float))
    qm_averaged = np.average( (qm_prev,qm_next) )
    Q_m.append(qm_averaged)
    # Ice
    qi = db_select(sql = 'SELECT Q_ice FROM '+drainage_basins_AND_midpoint+' where river = "'+rivername+'"')
    try:
      Q_i.append(float(np.array(qi, dtype=float)))
    except:
      Q_i.append(0.) # In this case, there is no entry, so no significant ice
    # Time
    age.append(float(self.mean_age[i]))
  # dt (because here we have the first time step and therefore all dt)
  dt_mid = list(-1*np.diff(np.array(self.mean_age, dtype=float)))
  # Create array and save output
  discharge_outputs = np.array([age, Q_m, Q_i, dt_mid]).transpose()
  np.savetxt('/home/awickert/Desktop/'+model+'_'+rivername+'_discharge_unsmoothed.txt', discharge_outputs, fmt='%d')
  np.savetxt('/home/awickert/Documents/geology_docs/papers/Working copies/MississippiBasinGlacial/figures/discharge_at_mouth/'+model+'_'+rivername+'_discharge_unsmoothed.txt', discharge_outputs, fmt='%d')

def get_basin_ice_volume_sea_level(self, basin='Mississippi'):
  """
  DO NOT USE! PRETENDS THAT A SHRINKING BASIN IS THE SAME AS LOSS OF ICE. OOPS!
  """
  print "Getting time series basin ice volume!"
  print basin, 'basin'
  print ""
  from grass.script import vector as vect
  self.V_ice_km3 = []
  for age in self.ages:
    print age
    drainage_basins = 'drainage_basins_' + age
    drainage_basins_cols = vect.vector_db_select(drainage_basins)['columns']
    river_name_column = (np.array(drainage_basins_cols) == 'river').nonzero()[0][0]
    V_ice_km3_column = (np.array(drainage_basins_cols) == 'V_ice_km3').nonzero()[0][0]
    drainage_basins_lol = vect.vector_db_select(drainage_basins)['values'].values()
    drainage_basins_array = np.array(drainage_basins_lol)
    self.V_ice_km3.append(float(drainage_basins_array[drainage_basins_array[:,river_name_column] == basin, V_ice_km3_column][0]))
  print "Converting to approximate sea level equivalent [m]"
  print "using modern sea surface area"
  self.SLeq_m = np.array(self.V_ice_km3) * 1000 / 361000000

# 4. Flood DEM depressions
# r.fill.dir in old/update_lakes_old3.py
def flood_basins(self):
  print "*************"

  for age in self.ages:
    print age
    
    tpin = 'tpin_' + age
    flooded = 'floded_' + age
    flowdir = 'flowdir_' + age # Save flow direction - maybe useful
    accum = 'accum_' + age # Actually save flow accumulation too - why not
    
    print "Initializing flooded raster"
    
    # "s" flag for single flow direction
    grass.run_command('r.terraflow', elevation=tpin, filled=flooded, direction=flowdir, swatershed='swatershed_tmp', accumulation=accum, tci='tci_tmp', flags='s', overwrite=True)
    grass.run_command('r.null', map=flooded, null='-9999') # Replace nulls with elev = -9999: < SL, which is what matters
    # mapcalc('tmp = ' + flooded + ' > 0')
    
  # Once we break out of that loop, cleanup:
  grass.run_command('g.remove', rast='swatershed_tmp')
  grass.run_command('g.remove', rast='tci_tmp')

# 5. Get depths of enclosed basins
# Must remove oceans and other areas that are flooded by JXM group model
def basins(self):
  print "*************"
  print "Obtaining unaccounted-for depths of enclosed basins"

  for age in self.ages:
    print age
    flooded = 'flooded_' + age
    tpin = 'tpin_' + age
    modelWaterDepth = 'modelWaterDepth_' + age
    basins = 'basins_' + age # Will be basins not counting the JXM group model lake areas
    # Don't allow any negative areas
    # This gets rid of oceans and already-accounted-for lakes
    mapcalc(basins + ' = max(((' + flooded + ' - ' + tpin + ') - ' + modelWaterDepth + '), 0)')

# 6. Calculate lake depths by subtraction
def lakes_ice_eq(self):
  print "*************"
  print "Calculating lake depths compared to present"

  # Start with modern lakes - need to subtract these out
  # In part to get rid of non-lake basins
  # In part to get rid of present lakes which are not considered loads because 
  # they are part of our frame of reference
  for age in self.ages:
    print age
    basins = 'basins_' + age
    dlakes_ice = 'dlakes_ice_' + age # Lake depth difference
    mapcalc(dlakes_ice + ' = 1.09051254 * (' + basins + ' - basins_000_0k)')
    # Floating point math issue:
    #mapcalc(dlakes_ice + ' = (' + str(self.rho_water) + ' / ' + str(self.rho_ice) + ') * (' + basins + ' - ' + lakes_now + ')')

# 7. Multiply these by cell sizes to get lake volumes in km**3
# Not necessary for coupling, but of scientific interest
def lakes_volume(self):
  print "*************"
  print "Calculating lake volumes compared to present"

  for age in self.ages:
    print age
    dlakes_ice = 'dlakes_ice_' + age # Lake depth difference
    dVlakes = 'dVlakes_' + age # Lake volume difference, km**3
    mapcalc(dVlakes + ' = 0.917 * (' + dlakes_ice + ' / 1000) * cellsize_km2')

# 8. Export new "ice" ( = lake water), summed with ICE-5G, to file
# This is for next run of JXM group's model
def export_new_ice(self):
  print "*************"
  print 'Making and exporting new "ice" files for next iteration with geophysical model'
  for age in self.ages:
    print age
    ice_plus_lake = 'ice_plus_lake' + age
    ice = 'ice_' + age
    dlakes_ice = 'dlakes_ice' + age # Lakes not present now
    mapcalc(ice_plus_lake + ' = ' + ice + ' + ' + dlakes_ice)
    #grass.run_command('r.out.ascii' , FILL OUT THE REST OF THIS FOR INPUTS AND PATH

# 9. Sum map to get full lake volume
# Not necessary until last step
def lakes_volume_sum(self):
  print "*************"
  print "Summing lake volumes through time"
  self.dVlakes_tot = []
  for age in self.ages:
    print age
    dVlakes = 'dVlakes_' + age # Lake volume difference, km**3
    self.dVlakes_tot.append(float(grass.parse_command('r.sum', rast=dVlakes)['SUM']))
  #self.dz_ocean_eq = (np.array(dVlakes_tot) / (4. * np.pi * (2./3.) * 6378.**2)) * 1000 # [m] Ocean ~ 2/3 of land surface
  # Marshall and Clarke value (converted to km**2, then * 1000 for m)
  self.dz_lake_ocean_eq = (np.array(self.dVlakes_tot) / 3.62E8) * 1000 # [m] km**3 / km**2; km--> m
  print self.dz_lake_ocean_eq
  import matplotlib.pyplot as plt
  plt.plot(self.ages_numeric, self.dz_lake_ocean_eq,'k-')
  plt.xlabel('Age [ka]', fontsize=16)
  plt.ylabel('Sea level equivalent lake volume [m]', fontsize=16)
  plt.show()

# 10. Get global ice volume with time (km**3)
def ice_volume(self):
  print "*************"
  print "Calculating ice volumes compared to present"
  for age in self.ages:
    print age
    ice = 'ice_' + age # Ice thickness difference
    dVice = 'dVice_' + age # Ice volume difference, km**3
    mapcalc(dVice + ' = (' + ice + ' - ice_000_0k) / 1000 * cellsize_km2')


def ice_volume_sum(self):
  grass.run_command('g.region', rast='cellsize_meters2')
  area_m2 = garray.array()
  area_m2.read('cellsize_meters2')
  h_i = garray.array()
  Vice_list = []
  for age in self.ages:
    print age
    ice = 'ice_' + age # Lake volume difference, km**3
    h_i.read(ice)
    Vice = np.sum(h_i * area_m2)
    Vice_list.append(Vice)
    print '  ', Vice/1E9, 'km**3'
  Vice_array = np.array(Vice_list)
  # Marshall and Clarke value for ocean area
  # (converted to km**2, then * 1000 for m)
  SLE = Vice_array * (0.917 / (3.62E8 * 1E6))
  print SLE
  out = np.vstack((self.ages_numeric, Vice_array, SLE)).transpose()
  loc = grass.gisenv()['LOCATION_NAME']
  outname = 'SLE_out/' + loc + '_SLE.txt'
  np.savetxt(outname, out)

# 11. Sum ice volume and plot with lakes
# Not necessary until last step
def ice_volume_sum_old(self):
  print "*************"
  print "Summing ice volume through time"
  self.Vice_list = []
  for age in self.ages:
    print age
    Vice = 'Vice_' + age # Lake volume difference, km**3
    univar_out = grass.parse_command('r.univar', map=ice)['SUM']
    kw = np.array( univar_out.keys()[0].split('|') )
    val = np.array( univar_out.keys()[1].split('|') )
    Vice = val[kw == 'sum']
    self.Vice_list.append(float())
  # Marshall and Clarke value for ocean area (converted to km**2, then * 1000 for m)
  # Also converted from ice volume to water volume
  self.Vice_tot = np.array(self.Vice_tot)
  self.z_ice_ocean_eq = (self.Vice_tot * 0.917 / 3.62E8) * 1000 # [m] km**3 / km**2; km--> m
  print self.z_ice_ocean_eq
  import matplotlib.pyplot as plt
  plt.plot(self.ages_numeric, self.dz_ice_ocean_eq,'k-')
  plt.xlabel('Age [ka]', fontsize=16)
  plt.ylabel('Sea level equivalent ice volume [m]', fontsize=16)
  plt.show()

# 12. Plot ice and lake volumes together
#  "and plotting with lake volume history"
