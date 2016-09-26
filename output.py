# All plotting, etc, here
# ADW, 25 NOV 2015

#save=False; show=True; legend=False; bigtitle=False; alltimes=False; onefig=True; ICE=None

# Import all of the built-in modules that "lakes" imports before its class 
# definition
from drainage import *


def basin_discharge_save(self):
  """
  Clone of basin_discharge_plots, for the pre-plotting part
  """
  # OLD: Qmodern = dict(Mississippi=18430., Mackenzie=9910., Colorado=665., Rio_Grande=150., Hudson=620., Columbia=7500., Susquehanna=1082., Saint_Lawrence=12000., Hudson_Strait=30900., Yukon=6428.)
  Qmodern = dict(Mississippi=18430., Mackenzie=9910., Colorado=665., Rio_Grande=570., Hudson=590., Columbia=7500., Susquehanna=1082., Saint_Lawrence=14400., Hudson_Strait=30900., Yukon=6428.)

  try:
    ICE = self.ICE
  except:
    # hack-ey solution here to plotting
    ICE = grass.parse_command('g.gisenv')['LOCATION_NAME']
    ICE = re.findall("[a-zA-Z0-9]+", ICE)[0]
  # Check if None-type too -- also hack-ey, there is so much cleaning to do.
  if ICE:
    pass
  else:
    # hack-ey solution here to plotting
    ICE = grass.parse_command('g.gisenv')['LOCATION_NAME']
    ICE = re.findall("[a-zA-Z0-9]+", ICE)[0]
  # LATER -- GET RID OF ICE DUMMY VAR
  self.ICE = ICE

  # River list
  self.riversAll = np.array(list(set(list(grass.parse_command('v.db.select', map='river_mouth_regions', columns='river', flags="c")))))
  self.rnumAll = np.arange(1,len(self.riversAll)+1) # Is this OK? what about database table rnum?

  # Create dicts to hold lists of Q_i, Q_m, Q_t
  # WHICH RIVERS WE USE ARE PRE-DEFINED HERE; CLUNKY TO CHANGE!
  self.Q_i = dict(Mississippi=[], Mackenzie=[], Colorado=[], Rio_Grande=[], Hudson=[], Columbia=[], Susquehanna=[], Saint_Lawrence=[], Hudson_Strait=[], Yukon=[])
  self.Q_m = dict(Mississippi=[], Mackenzie=[], Colorado=[], Rio_Grande=[], Hudson=[], Columbia=[], Susquehanna=[], Saint_Lawrence=[], Hudson_Strait=[], Yukon=[])
  self.Q_t = dict(Mississippi=[], Mackenzie=[], Colorado=[], Rio_Grande=[], Hudson=[], Columbia=[], Susquehanna=[], Saint_Lawrence=[], Hudson_Strait=[], Yukon=[])
  #shortRivers = np.array(['Mississippi River', 'Susquehanna River', 'Hudson River', 'Saint Lawrence River', 'Hudson Strait', 'Mackenzie River', 'Columbia River', 'Colorado River', 'Rio Grande'])
  #srnu = np.array(['Mississippi', 'Susquehanna', 'Hudson', 'Saint_Lawrence', 'Hudson_Strait', 'Mackenzie', 'Columbia', 'Colorado', 'Rio_Grande'])
  #self.Q_i = dict(Mississippi=[], Mackenzie=[], Columbia=[], Saint_Lawrence=[], Hudson_Strait=[], Other_North=[], Volga=[], White_Sea = [], Daugava_and_Baltic=[])
  #self.Q_m = dict(Mississippi=[], Mackenzie=[], Columbia=[], Saint_Lawrence=[], Hudson_Strait=[], Other_North=[], Volga=[], White_Sea = [], Daugava_and_Baltic=[])
  #self.Q_t = dict(Mississippi=[], Mackenzie=[], Columbia=[], Saint_Lawrence=[], Hudson_Strait=[], Other_North=[], Volga=[], White_Sea = [], Daugava_and_Baltic=[])
  srnu = self.Q_t.keys()
  shortRivers = []
  for i in range(len(srnu)): # space for list of names
    shortRivers.append(srnu[i].replace('_', ' '))

  # Update rivers and rnum to include only these desired outputs
  self.rivers = np.array(self.Q_t.keys())
  for i in range(len(self.rivers)): # space for list of names
    self.rivers[i] = self.rivers[i].replace('_', ' ')
  self.rnum = []
  for river in self.rivers:
    self.rnum.append(self.rnumAll[self.riversAll == river][0])
  
  for age in self.ages:
    print age
    drainage_basins = 'drainage_basins_' + age
    try:
      drainage_basins_cols = vect.vector_db_select(drainage_basins)['columns']
      drainage_basins_lol = vect.vector_db_select(drainage_basins)['values'].values()
      river_name_column = (np.array(drainage_basins_cols) == 'river').nonzero()[0][0]
      meltwater_discharge_column = (np.array(drainage_basins_cols) == 'Q_ice').nonzero()[0][0]
      meteoric_discharge_column = (np.array(drainage_basins_cols) == 'Q_meteoric_uncorrected').nonzero()[0][0]
      #total_discharge_column = (np.array(drainage_basins_cols) == 'Q_total').nonzero()[0][0]
    except:
      print 'Error at', age
      drainage_basins_cols = None
      for river in self.rivers:
        rnu = river.replace(' ', '_') # Underscore for dict
        # Some rivers not appearing in database table, so have to do this as a quick fix!
        # Maybe clashing rules files from running many of these at once?
        self.Q_i[rnu].append(np.nan)
        self.Q_m[rnu].append(np.nan)
        self.Q_t[rnu].append(np.nan)

    if drainage_basins_cols:
      for river in self.rivers:
        rnu = river.replace(' ', '_') # Underscore for dict
        # Some rivers not appearing in database table, so have to do this as a quick fix!
        # Maybe clashing rules files from running many of these at once?
        self.Q_i[rnu].append(np.nan)
        self.Q_m[rnu].append(np.nan)
        self.Q_t[rnu].append(np.nan)
        for row in drainage_basins_lol:
          if row[river_name_column] == river:
            #print row
            #if row[total_discharge_column] != '':
            self.Q_i[rnu][-1] = row[meltwater_discharge_column]
            self.Q_m[rnu][-1] = row[meteoric_discharge_column]
            
  self.Q_meteoric_simulated_now = []
  for river in self.rivers:
    rnu = river.replace(' ', '_') # Underscore for dict
    self.Q_m[rnu] = np.array(self.Q_m[rnu]).astype(float)
    self.Q_t[rnu] = np.array(self.Q_i[rnu]).astype(float)
    try:
      self.Q_meteoric_simulated_now.append(self.Q_m[rnu][-1])
      self.Q_m[rnu] *= float(Qmodern[rnu]) / np.mean(np.array(self.Q_m[rnu]).astype(float)[self.ages_numeric < 3000])
      self.Q_t[rnu] += np.array(self.Q_m[rnu]).astype(float) # now forced to match!
    except:
      self.Q_meteoric_simulated_now.append(np.nan)

  outvars = {'ages': self.ages, 'ages_numeric': self.ages_numeric, \
             'Q_i': self.Q_i, 'Q_m': self.Q_m, 'Q_t': self.Q_t}
  np.save(self.ICE+'.npy', outvars)


def basin_discharge_plots(self, save=False, show=True, legend=True, bigtitle=False, alltimes=False, onefig=True, ICE=None):
  from grass.script import vector as vect
  import matplotlib.pyplot as plt
  
  Qmodern = dict(Mississippi=18430., Mackenzie=9910., Colorado=665., Rio_Grande=150., Hudson=620., Columbia=7500., Susquehanna=1082., Saint_Lawrence=12000., Hudson_Strait=30900., Yukon=6428.)

  if ICE:
    pass
  else:
    try:
      ICE = self.ICE
    except:
      # hack-ey solution here to plotting
      ICE = grass.parse_command('g.gisenv')['LOCATION_NAME']
      ICE = re.findall("[a-zA-Z0-9]+", ICE)[0]
  self.ICE = ICE

  # River list
  self.riversAll = np.array(list(set(list(grass.parse_command('v.db.select', map='river_mouth_regions', columns='river', flags="c")))))
  self.rnumAll = np.arange(1,len(self.riversAll)+1) # Is this OK? what about database table rnum?

  # Create dicts to hold lists of Q_i, Q_m, Q_t
  # WHICH RIVERS WE USE ARE PRE-DEFINED HERE; CLUNKY TO CHANGE!
  self.Q_i = dict(Mississippi=[], Mackenzie=[], Colorado=[], Rio_Grande=[], Hudson=[], Columbia=[], Susquehanna=[], Saint_Lawrence=[], Hudson_Strait=[], Yukon=[])
  self.Q_m = dict(Mississippi=[], Mackenzie=[], Colorado=[], Rio_Grande=[], Hudson=[], Columbia=[], Susquehanna=[], Saint_Lawrence=[], Hudson_Strait=[], Yukon=[])
  self.Q_t = dict(Mississippi=[], Mackenzie=[], Colorado=[], Rio_Grande=[], Hudson=[], Columbia=[], Susquehanna=[], Saint_Lawrence=[], Hudson_Strait=[], Yukon=[])
  #shortRivers = np.array(['Mississippi River', 'Susquehanna River', 'Hudson River', 'Saint Lawrence River', 'Hudson Strait', 'Mackenzie River', 'Columbia River', 'Colorado River', 'Rio Grande'])
  #srnu = np.array(['Mississippi', 'Susquehanna', 'Hudson', 'Saint_Lawrence', 'Hudson_Strait', 'Mackenzie', 'Columbia', 'Colorado', 'Rio_Grande'])
  #self.Q_i = dict(Mississippi=[], Mackenzie=[], Columbia=[], Saint_Lawrence=[], Hudson_Strait=[], Other_North=[], Volga=[], White_Sea = [], Daugava_and_Baltic=[])
  #self.Q_m = dict(Mississippi=[], Mackenzie=[], Columbia=[], Saint_Lawrence=[], Hudson_Strait=[], Other_North=[], Volga=[], White_Sea = [], Daugava_and_Baltic=[])
  #self.Q_t = dict(Mississippi=[], Mackenzie=[], Columbia=[], Saint_Lawrence=[], Hudson_Strait=[], Other_North=[], Volga=[], White_Sea = [], Daugava_and_Baltic=[])
  srnu = self.Q_t.keys()
  shortRivers = []
  for i in range(len(srnu)): # space for list of names
    shortRivers.append(srnu[i].replace('_', ' '))

  # Update rivers and rnum to include only these desired outputs
  self.rivers = np.array(self.Q_t.keys())
  for i in range(len(self.rivers)): # space for list of names
    self.rivers[i] = self.rivers[i].replace('_', ' ')
  self.rnum = []
  for river in self.rivers:
    self.rnum.append(self.rnumAll[self.riversAll == river][0])
  
  for age in self.ages:
    print age
    drainage_basins = 'drainage_basins_' + age
    try:
      drainage_basins_cols = vect.vector_db_select(drainage_basins)['columns']
      drainage_basins_lol = vect.vector_db_select(drainage_basins)['values'].values()
      river_name_column = (np.array(drainage_basins_cols) == 'river').nonzero()[0][0]
      meltwater_discharge_column = (np.array(drainage_basins_cols) == 'Q_ice').nonzero()[0][0]
      meteoric_discharge_column = (np.array(drainage_basins_cols) == 'Q_meteoric_uncorrected').nonzero()[0][0]
      #total_discharge_column = (np.array(drainage_basins_cols) == 'Q_total').nonzero()[0][0]
    except:
      print 'Error at', age
      drainage_basins_cols = None
      for river in self.rivers:
        rnu = river.replace(' ', '_') # Underscore for dict
        # Some rivers not appearing in database table, so have to do this as a quick fix!
        # Maybe clashing rules files from running many of these at once?
        self.Q_i[rnu].append(np.nan)
        self.Q_m[rnu].append(np.nan)
        self.Q_t[rnu].append(np.nan)

    if drainage_basins_cols:
      for river in self.rivers:
        rnu = river.replace(' ', '_') # Underscore for dict
        # Some rivers not appearing in database table, so have to do this as a quick fix!
        # Maybe clashing rules files from running many of these at once?
        self.Q_i[rnu].append(np.nan)
        self.Q_m[rnu].append(np.nan)
        self.Q_t[rnu].append(np.nan)
        for row in drainage_basins_lol:
          if row[river_name_column] == river:
            #print row
            #if row[total_discharge_column] != '':
            self.Q_i[rnu][-1] = row[meltwater_discharge_column]
            self.Q_m[rnu][-1] = row[meteoric_discharge_column]
            
            """
            else:
              print age, '!!!', river
              self.Q_i[rnu].append(np.nan)
              self.Q_m[rnu].append(np.nan)
              self.Q_t[rnu].append(np.nan)
            if river == 'Hudson':
              print age, row[total_discharge_column]
            """
  
  self.Q_meteoric_simulated_now = []
  for river in self.rivers:
    rnu = river.replace(' ', '_') # Underscore for dict
    self.Q_m[rnu] = np.array(self.Q_m[rnu]).astype(float)
    self.Q_t[rnu] = np.array(self.Q_i[rnu]).astype(float)
    try:
      self.Q_meteoric_simulated_now.append(self.Q_m[rnu][-1])
      self.Q_m[rnu] *= float(Qmodern[rnu]) / np.mean(np.array(self.Q_m[rnu]).astype(float)[self.ages_numeric < 3000])
      self.Q_t[rnu] += np.array(self.Q_m[rnu]).astype(float) # now forced to match!
    except:
      self.Q_meteoric_simulated_now.append(np.nan)

  self.output.HighDischargeTimes(self)


  outvars = {'ages': self.ages, 'ages_numeric': self.ages_numeric, \
             'Q_i': self.Q_i, 'Q_m': self.Q_m, 'Q_t': self.Q_t}
  np.save(self.ICE+'.npy', outvars)

  sys.exit("Mostly, need to fix long vs. short names, but just use Qhist_all.py after running this")











  if alltimes:
    agei=0
    if onefig:
      print "OneFig doesn't function currently for all times"
    #for age in self.ages:
      for i in range(len(self.rivers)):
        river = self.rivers[i]
        print river
        rnu = river.replace(' ', '_') # Underscore for dict  
        plt.figure(i, figsize=(8,5))
        if bigtitle:
          plt.plot(self.ages_numeric/1000., self.Q_t[rnu],'k-', linewidth=10, label='Total')
          plt.plot(self.ages_numeric/1000., self.Q_m[rnu],'g-', linewidth=4, label='Meteoric')
          plt.plot(self.ages_numeric/1000., self.Q_i[rnu], 'b--', linewidth=4, label='Melt')
          plt.plot(self.ages_numeric/1000.[1+agei], self.Q_t[rnu][agei], 'ro', markersize=20)
        elif onefig:
          pass
        else:
          plt.plot(self.ages_numeric/1000., self.Q_t[rnu],'k-', linewidth=6, label='Total')
          plt.plot(self.ages_numeric/1000., self.Q_m[rnu],'g-', linewidth=2, label='Meteoric')
          plt.plot(self.ages_numeric/1000., self.Q_i[rnu], 'b--', linewidth=2, label='Melt')
          plt.plot(self.ages_numeric/1000.[1+agei], self.Q_t[rnu][agei], 'ro', markersize=12)
        if legend:
          plt.legend(loc='upper left')
        elif onefig:
          pass
        if bigtitle:
          plt.title(river, fontsize=48, fontweight='bold')
          plt.xlabel('Age [ka]', fontsize=24)
          plt.ylabel('Discharge [m$^3$/s]', fontsize=24)
        elif onefig:
          pass
        else:
          plt.title(river, fontsize=24, fontweight='bold')
          plt.xlabel('Age [ka]', fontsize=16)
          plt.ylabel('Discharge [m$^3$/s]', fontsize=16)
        fig = plt.gcf()
        if bigtitle:
          fig.subplots_adjust(top=0.84)
          for label in plt.gca().get_xticklabels() + plt.gca().get_yticklabels():
            label.set_fontsize(20)
          fig.subplots_adjust(left=0.22)
          fig.subplots_adjust(bottom=0.15)
        elif onefig:
          pass
        else:
          fig.subplots_adjust(left=0.16)      
          fig.subplots_adjust(bottom=0.12)
        if save:
          if bigtitle:
            plt.savefig('alltimes/dischargefigs_bigtitles/'+rnu+'_'+age)
          elif legend:
            plt.savefig('alltimes/dischargefigs_legends/'+rnu+'_'+age)
          else:
            plt.savefig('alltimes/dischargefigs/'+rnu+'_'+age)
          plt.clf()
      agei+=1
      if show:
        print "Not allowing hundreds of figures to be shown - poor computer!"
  

  else:
    if onefig:
      print "onefig ignores bigtitle"
      fig = plt.figure(figsize = (10, 14))
      fig.text((1+.17-.1)/2., 0.04, 'Age [ka]', ha='center', va='center', fontsize=24, fontweight='bold')
      fig.text(0.06, 0.5, 'Discharge [m$^3$s$^{-1}$]', ha='center', va='center', rotation='vertical', fontsize=24, fontweight='bold')
      #ax0 = fig.add_subplot(111)
      #ax0.spines['top'].set_color('none')
      #ax0.spines['bottom'].set_color('none')
      #ax0.spines['left'].set_color('none')
      #ax0.spines['right'].set_color('none')
      #ax0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
      #ax0.set_xlabel('Age [ka]', fontsize=16)
      #ax0.set_ylabel('Discharge [m$^3$s$^{-1}$]', fontsize=16)
      iof = 1
      for i in range(len(shortRivers)):
        print shortRivers[i]
        rnu = srnu[i].replace(' ', '_') # Underscore for dict  
        ax = fig.add_subplot(5, 2, iof)
        try:
          for row in self.EnhancedQ[shortRivers[i]]:
            ax.axvspan(row[0], row[1], facecolor='0.8', linewidth='0')
        except:
          print "No geological constraints entered for", shortRivers[i]
        ax.plot(self.ages_numeric/1000., self.Q_t[rnu],'k-', linewidth=6, label='Total discharge')
        ax.plot(self.ages_numeric/1000., self.Q_i[rnu], 'b--', linewidth=2, label='Meltwater discharge')
        ax.plot(self.ages_numeric/1000., self.Q_m[rnu],'g-', linewidth=2, label='Meteoric water\n(P-ET) discharge')
        ax.set_title(shortRivers[i], fontsize=16, fontweight='bold')
        ax.set_xlim((0, 20))
        ax.set_ylim((0, ax.get_ylim()[-1]))
        try:
          plt.plot(self.ages_numeric[-1]/1000., self.Q_meteoric_simulated_now[(self.rivers == rnu).nonzero()[0][0]], 'go', markersize=10, zorder=-100)
        except:
          print "No discharge constraints entered for", shortRivers[i]
        # set y-lim based on model outputs
        """
        if rnu == 'Mississippi':
          ax.set_ylim((0, 200000))
        elif rnu == 'Susquehanna':
          ax.set_ylim((0, 40000))
        elif rnu == 'Hudson':
          ax.set_ylim((0, 120000))
        elif rnu == 'Saint_Lawrence':
          ax.set_ylim((0, 160000))
        elif rnu == 'Hudson_Strait':
          ax.set_ylim((0, 250000))
        elif rnu == 'Mackenzie':
          ax.set_ylim((0, 160000))
        elif rnu == 'Columbia':
          ax.set_ylim((0, 18000))
        elif rnu == 'Colorado':
          ax.set_ylim((0, 2000))
        elif rnu == 'Rio_Grande':
          ax.set_ylim((0, 600))
        """
        if iof == 9:
          ax.legend(loc = 'center', bbox_to_anchor = (1.71, 0.5), fontsize=14)
        iof += 1
      fig.tight_layout()
      fig.subplots_adjust(left=.17, bottom=.08, right=None, top=None, wspace=None, hspace=.4)
      if save:
        plt.savefig('all_rivers_'+ICE+'.png')
        plt.savefig('all_rivers_'+ICE+'.pdf')
    else:
      for i in range(len(self.rivers)):
        river = self.rivers[i]
        rnu = river.replace(' ', '_') # Underscore for dict  
        plt.figure(i, figsize=(8,5))
        if bigtitle:
          plt.plot(self.ages_numeric/1000., self.Q_t[rnu],'k-', linewidth=10, label='Total')
          plt.plot(self.ages_numeric/1000., self.Q_m[rnu],'g-', linewidth=4, label='Meteoric')
          plt.plot(self.ages_numeric/1000., self.Q_i[rnu], 'b--', linewidth=4, label='Melt')
        else:
          plt.plot(self.ages_numeric/1000., self.Q_t[rnu],'k-', linewidth=6, label='Total')
          plt.plot(self.ages_numeric/1000., self.Q_m[rnu],'g-', linewidth=2, label='Meteoric')
          plt.plot(self.ages_numeric/1000., self.Q_i[rnu], 'b--', linewidth=2, label='Melt')
        if legend:
          plt.legend(loc='upper left')
        if bigtitle:
          plt.title(river, fontsize=48, fontweight='bold')
          plt.xlabel('Age [ka]', fontsize=24)
          plt.ylabel('Discharge [m$^3$/s]', fontsize=24)
        else:
          plt.title(river, fontsize=24, fontweight='bold')
          plt.xlabel('Age [ka]', fontsize=16)
          plt.ylabel('Discharge [m$^3$/s]', fontsize=16)
        fig = plt.gcf()
        if bigtitle:
          fig.subplots_adjust(top=0.84)
          for label in plt.gca().get_xticklabels() + plt.gca().get_yticklabels():
            label.set_fontsize(20)
            fig.subplots_adjust(left=0.22)
          fig.subplots_adjust(bottom=0.15)
        else:
          fig.subplots_adjust(left=0.16)      
          fig.subplots_adjust(bottom=0.12)
        if save:
          if bigtitle:
            plt.savefig('dischargefigs_bigtitles/'+rnu)
          elif legend:
            plt.savefig('dischargefigs_legends/'+rnu)
          else:
            plt.savefig('dischargefigs/'+rnu)
    if show:
      plt.show()
    else:
      if onefig:
        plt.close()
      else:
        for i in range(len(self.rivers)):
          plt.figure(i)
          #plt.clf()
          plt.close()
  
  # Has to be done only once
  # Though should in the future put into array above
  # np.savetxt('Q_meteoric_simulated_now.txt', np.array(self.Q_meteoric_simulated_now).astype(float))

def HighDischargeTimes(self):
  """
  Timing of enhanced (melt)water outflow
  """
  self.EnhancedQ = { 'Susquehanna River':          np.array([[26.5, 19.9],
                                                             [18.4, 15.9]]),

                     'Hudson River':               np.array([[19.9, 18.4],
                                                             [15.9, 13.0]]),

                     'Saint Lawrence River':       np.array([[12.8,  8.8]]),

                     'Mackenzie River':            np.array([[13.0,  9.3]]),

                     'Mississippi River':          np.array([[26.5, 18.9],
                                                             [18.4, 12.8],
                                                             [11.7, 11.6],
                                                             [10.8, 10.7]]),

                     'Hudson Strait':              np.array([[25.7, 24.0],
                                                             [18.2, 17.7],
                                                             [12.8, 12.3],
                                                             [11.4, 10.9],
                                                             [10.0,  9.4],
                                                             [8.65, 8.55],
                                                             [8.25, 8.15]]),

                     'Columbia River':             np.array([[36.1, 10.6]]),

                     'Colorado River':             np.array([[23.0, 15.0],
                                                             [12.9, 11.7]]),

                     'Rio Grande':                 np.array([[21.0, 13.5],
                                                             [12.9, 11.7]])
  }


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
  ax.plot(self.ages_numeric/1000., self.Q_t, 'k-', label="Total", linewidth=3)
  ax.plot(self.ages_numeric/1000., self.Q_i, 'b-', label="Ice sheet", linewidth=2)
  ax.plot(self.ages_numeric/1000., self.Q_m, 'g-', label="Meteoric", linewidth=2)
  ax.set_title(self.rivername, fontsize=16)
  ax.set_xlabel("Time [ka BP]", fontsize=16)
  ax.set_ylabel(r"Discharge [m$^3$/s]", fontsize=16)
  ax.legend(loc='upper left')
  plt.show()

def output_Q(self, rivername, outpath=''):
  # Run after the plotting step, for now
  # Also keeping on time-steps instead of trying to make more precise midpoints
  # like I did for MWP-1A paper
  # Also have to check on renormalization
  output = np.vstack((self.ages, self.Q_m[rivername], self.Q_i[rivername], self.Q_t[rivername])).transpose().astype(float)
  np.savetxt(outpath+'/'+ICE+'.txt', output)


