import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

EnhancedQ = {      'Susquehanna River':          np.array([[26.5, 19.9],
                                                           [18.4, 15.9]]),

                   'Hudson River':               np.array([[19.9, 18.4],
                                                           [15.9, 13.0]]),

                   'Saint Lawrence River':       np.array([[12.8,  8.8]]),

                   'Mackenzie River':            np.array([[13.6,  11.5]]),

                   'Mississippi River':          np.array([[26.5, 19.9],
                                                           [19.3, 18.8],
                                                           [18.4, 15.9],
                                                           [13.6, 12.8],
                                                           [11.8, 11.3],
                                                           [11.0, 10.6]]),

                   'Hudson Strait':              np.array([[16.8, 16.3],
                                                           [15.9, 15.7],
                                                           [11.5, 11.3],
                                                           [11.0, 10.6],
                                                           [9.7,  9.2],
                                                           [8.65, 8.55],
                                                           [8.25, 8.15]]),

                   'Columbia River':             np.array([[22.9, 16.4]]),

                   'Colorado River':             np.array([[23.0, 15.0],
                                                           [12.9, 11.7]]),

                   'Rio Grande':                 np.array([[21.0, 13.5],
                                                           [12.9, 11.7]])
}

# OLD H0 AGES (and H1 and more)
# New 11.5-11.3 H1 + 11.4-10.9 Ungava Bay
#                   'Hudson Strait':              np.array([[25.7, 24.0],
#                                                           [18.2, 17.7],
#                                                           [12.8, 12.3],
#                                                           [11.4, 10.9],
#                                                           [10.0,  9.4],
#                                                           [8.65, 8.55],
#                                                           [8.25, 8.15]]),
#
#                   'Mackenzie River':            np.array([[13.6,  11.5],
#                                                           [10.8,  10.7]]),
#
#                   'Mississippi River':          np.array([[26.5, 19.9],
#                                                           [19.3, 18.8],
#                                                           [18.4, 15.9],
#                                                           [13.6, 12.8],
#                                                           [11.7, 11.6],
#                                                           [10.8, 10.7]]),
#

# Generalized enhanced discharge: smoothes over general times when there should
# have been more discharge than at present
Enh_Q_gen = {      'Susquehanna River':          np.array([[26.5, 19.9],
                                                           [18.4, 15.9]]),

                   # Guessing based on ice extent that Q would be higher
                   'Hudson River':               np.array([[26.5, 18.4],
                                                           [15.9, 13.0]]),

                   'Saint Lawrence River':       np.array([[12.8,  8.8]]),

                   'Mackenzie River':            np.array([[15.6,  11.0]]),

                   'Mississippi River':          np.array([[26.5, 12.8],
                                                           [11.8, 11.3],
                                                           [11.0, 10.6]]),

                   'Hudson Strait':              np.array([[16.8, 16.3],
                                                           [15.9, 15.7],
                                                           [11.5, 11.3],
                                                           [11.0, 10.6],
                                                           [9.7,  8.0],
                                                           [7.55,  7.35]]),

                   'Columbia River':             np.array([[36.1, 10.6]]),

                   'Colorado River':             np.array([[23.0, 15.0],
                                                           [12.9, 11.7]]),

                   'Rio Grande':                 np.array([[21.0, 13.5],
                                                           [12.9, 11.7]])
}

# OLD H0 AGES (and H1 and more)
#                   'Hudson Strait':              np.array([[25.7, 24.0],
#                                                           [18.2, 17.7],
#                                                           [12.8, 8.2]]),
#
#                   'Mackenzie River':            np.array([[15.6,  9.3]]),
#
#                   'Hudson Strait':              np.array([[16.8, 16.3],
#                                                           [15.9, 15.7],
#                                                           [11.5, 8.15]]),
#

Qmodern = dict(Mississippi=18430., Mackenzie=9910., Colorado=665., Rio_Grande=150., Hudson=620., Columbia=7500., Susquehanna=1082., Saint_Lawrence=12000., Hudson_Strait=30900., Yukon=6428.)

model_runs = ['ICE3G', 'ICE5G', 'ICE6G', 'ANU', 'G12']
model_runs_longnames = ['ICE-3G/VM1', 'ICE-5G/VM2', 'ICE-6G/VM5a', 'ANU', 'G12/VM2']
colors = ['green', 'blue', 'red', 'purple', 'orange']

Q = []
for mrun in model_runs:
  Q.append( np.load(mrun+'.npy').item() )

shortRivers = np.array(['Mississippi River', 'Susquehanna River', 'Hudson River', 'Saint Lawrence River', 'Hudson Strait', 'Mackenzie River', 'Columbia River', 'Colorado River', 'Rio Grande'])
srnu = np.array(['Mississippi', 'Susquehanna', 'Hudson', 'Saint_Lawrence', 'Hudson_Strait', 'Mackenzie', 'Columbia', 'Colorado', 'Rio_Grande'])
srnu_spaces = np.array(['Mississippi', 'Susquehanna', 'Hudson', 'Saint Lawrence', 'Hudson Strait', 'Mackenzie', 'Columbia', 'Colorado', 'Rio Grande'])
plotRivers = np.array(['Mississippi River', 'Susquehanna/Chesapeake', 'Hudson River', 'Saint Lawrence River', 'Hudson Strait', 'Mackenzie River', 'Columbia River', 'Colorado River', 'Rio Grande'])


# First, convert to floating-point
for i in range(len(Q)):
  for Qtype in ['Q_i', 'Q_m', 'Q_t']:
    for rnu in srnu:
      Q[i][Qtype][rnu] = np.array(Q[i][Qtype][rnu]).astype(float)

# interpolate throgh the NaN's:
# This is a bit of a stopgap while I figure out why these are not exporting
# correctly, though on the scale of this paper, should be fine -- and I should
# publish the raw values (including NaN's) as well, anyway
for i in range(len(Q)):
  for Qtype in ['Q_i', 'Q_m', 'Q_t']:
    for rnu in srnu:
      not_nan = np.logical_not(np.isnan(Q[i][Qtype][rnu]))
      try:
        #Q[i][Qtype][rnu] = f(Q[i]['ages_numeric'])
        f = interp1d(Q[i]['ages_numeric'][not_nan][::-1], Q[i][Qtype][rnu][not_nan][::-1])
        Q[i][Qtype][rnu][Q[i]['ages_numeric'] < 20000] = f(Q[i]['ages_numeric'][Q[i]['ages_numeric'] < 20000])
      except:
        pass

# Should be using the Q_meteoric_uncorrected column in the database tables
# instead of this -- works, but redundant and something that could break
Q_meteoric_simulated_now = np.genfromtxt('Q_meteoric_simulated_now.txt')

j = 0
for Qi in Q:
  print model_runs[j]
  fig = plt.figure(figsize = (10, 14))
  fig.text((1+.17-.1)/2., 0.04, 'Age [ka]', ha='center', va='center', fontsize=24, fontweight='bold')
  fig.text(0.06, 0.5, 'Discharge [m$^3$s$^{-1}$]', ha='center', va='center', rotation='vertical', fontsize=24, fontweight='bold')

  iof = 1
  for i in range(len(shortRivers)):
    print shortRivers[i]
    rnu = srnu[i].replace(' ', '_') # Underscore for dict  
    ax = fig.add_subplot(5, 2, iof)
    for row in Enh_Q_gen[shortRivers[i]]:
      ax.axvspan(row[0], row[1], facecolor='0.85', linewidth=0)
    for row in EnhancedQ[shortRivers[i]]:
      ax.axvspan(row[0], row[1], facecolor='0.65', linewidth=0)
    ax.plot(Qi['ages_numeric']/1000., Qi['Q_t'][rnu],'k-', linewidth=6, label='Total discharge')
    ax.plot(Qi['ages_numeric']/1000., Qi['Q_i'][rnu], 'b-', linewidth=3, label='Meltwater discharge')
    ax.plot(Qi['ages_numeric']/1000., Qi['Q_m'][rnu],'g-', linewidth=2, label='Meteoric water\n(P-ET) discharge')
    ax.set_title(plotRivers[i], fontsize=16, fontweight='bold')
    ax.set_xlim((0, 20))
    ax.set_ylim((0, ax.get_ylim()[-1]))
    plt.plot(Qi['ages_numeric'][-1]/1000., Q_meteoric_simulated_now[(np.array(Qmodern.keys()) == rnu).nonzero()[0][0]], 'go', markersize=10, zorder=-100, label='Uncorrected modern\nmodeled discharge')
    # set y-lim based on model outputs
    if iof == 9:
      ax.legend(loc = 'center', bbox_to_anchor = (1.71, 0.5), fontsize=14, numpoints=1)
    iof += 1
    if rnu == 'Mississippi':
      ax.set_ylim((0, 160000))
    elif rnu == 'Susquehanna':
      ax.set_ylim((0, 40000))
    elif rnu == 'Hudson':
      ax.set_ylim((0, 160000)) # Could be 140000, but this is better for comparison
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
      ax.set_ylim((0, 1200))
  fig.tight_layout()
  fig.subplots_adjust(left=.17, bottom=.08, right=None, top=None, wspace=None, hspace=.4)
  plt.savefig('all_rivers_'+model_runs[j]+'.png')
  plt.savefig('../figures/all_rivers_'+model_runs[j]+'.pdf')
  plt.close()
  j += 1
  


# Now, to plot all total Q's together
# And add Licciardi et al. (1999)
L1999 = np.genfromtxt('/home/awickert/Dropbox/Papers/InProgress/DrainageMethods/Licciardi1999_tables/TotalDischarge.tsv', skip_header=1)
j = 0
fig = plt.figure(figsize = (10, 14))
fig.text((1+.17-.1)/2., 0.04, 'Age [ka]', ha='center', va='center', fontsize=24, fontweight='bold')
fig.text(0.06, 0.5, 'Total Discharge (meltwater + meteoric water) [m$^3$s$^{-1}$]', ha='center', va='center', rotation='vertical', fontsize=24, fontweight='bold')
iof = 1
# Licciardi 1999
for rnu in srnu:
  ax = fig.add_subplot(5, 2, iof)
  if rnu == 'Mississippi':
    ax.plot(L1999[:,1]/1000., L1999[:,2], color='black', alpha=0.55, linewidth=2, label='Licciardi et al. (1999)', zorder=1.5)
  elif rnu == 'Hudson':
    ax.plot(L1999[:,1]/1000., L1999[:,4], color='black', alpha=0.55, linewidth=2, label='Licciardi et al. (1999)', zorder=1.5)
  elif rnu == 'Saint_Lawrence':
    ax.plot(L1999[:,1]/1000., L1999[:,3], color='black', alpha=0.55, linewidth=2, label='Licciardi et al. (1999)', zorder=1.5)
  elif rnu == 'Hudson_Strait':
    ax.plot(L1999[:,1]/1000., L1999[:,5], color='black', alpha=0.55, linewidth=2, label='Licciardi et al. (1999)', zorder=1.5)
  elif rnu == 'Mackenzie':
    pass # L1999 is for whole Arctic
  iof += 1
for Qi in Q:
  print model_runs[j]
  iof = 1
  for i in range(len(shortRivers)):
    #print shortRivers[i]
    rnu = srnu[i].replace(' ', '_') # Underscore for dict  
    ax = fig.add_subplot(5, 2, iof)
    for row in Enh_Q_gen[shortRivers[i]]:
      ax.axvspan(row[0], row[1], facecolor='0.85', linewidth=0)
    for row in EnhancedQ[shortRivers[i]]:
      ax.axvspan(row[0], row[1], facecolor='0.65', linewidth=0)
    ax.plot(Qi['ages_numeric']/1000., Qi['Q_t'][rnu], color=colors[j], alpha=0.6, linewidth=3, label=model_runs_longnames[j])
    ax.set_title(plotRivers[i], fontsize=16, fontweight='bold')
    ax.set_xlim((0, 20))
    ax.set_ylim((0, ax.get_ylim()[-1]))
    # set y-lim based on model outputs
    iof += 1
    if rnu == 'Mississippi':
      ax.set_ylim((0, 160000))
    elif rnu == 'Susquehanna':
      ax.set_ylim((0, 40000))
    elif rnu == 'Hudson':
      ax.set_ylim((0, 160000))
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
      ax.set_ylim((0, 1200))
  j += 1

# L1999
#  elif rnu == 'Rio_Grande':
#    # Legend based off of this one
ax.plot([1000, 1005], [-1500, -1500], color='black', alpha=0.55, linewidth=2, label='Licciardi et al. (1999)', zorder=1.5)

ax.legend(loc = 'center', bbox_to_anchor = (1.71, 0.5), fontsize=14)
fig.tight_layout()
fig.subplots_adjust(left=.17, bottom=.08, right=None, top=None, wspace=None, hspace=.4)

plt.savefig('all_rivers.png')
plt.savefig('../figures/all_rivers.pdf')
plt.close()

# Now the distributions
from scipy import stats
p_values = []
chisquare = []
D_values = [] # KS statistic
MeanQ_fieldQbig = []
MeanQ_fieldQsmall = []

for j in range(len(Q)):
  fig = plt.figure(figsize = (10, 14))
  fig.text((1+.17-.1)/2., 0.04, 'Age [ka]', ha='center', va='center', fontsize=24, fontweight='bold')
  fig.text(0.06, 0.5, 'Total Discharge (meltwater + meteoric water) [m$^3$s$^{-1}$]', ha='center', va='center', rotation='vertical', fontsize=24, fontweight='bold')
  Qi = Q[j]
  ages_ka = Qi['ages_numeric']/1000.
  Q_fieldQbig = []
  Q_fieldQsmall = []
  print model_runs[j]
  iof = 1
  for i in range(len(shortRivers)):
    #print shortRivers[i]
    rnu = srnu[i].replace(' ', '_') # Underscore for dict  
    fieldQbig = np.zeros(ages_ka.shape, dtype=bool)
    for row in Enh_Q_gen[shortRivers[i]]:
      fieldQbig += ( (ages_ka <= row[0]) * (ages_ka > row[1]) )
    Q_fieldQbig.append(Qi['Q_t'][rnu][fieldQbig])
    Q_fieldQsmall.append(Qi['Q_t'][rnu][fieldQbig == False])
    # rmv nan
    Q_fieldQbig[i] = Q_fieldQbig[i][np.isnan(Q_fieldQbig[i]) == False]
    Q_fieldQsmall[i] = Q_fieldQsmall[i][np.isnan(Q_fieldQsmall[i]) == False]
    
    try:
      big = np.histogram(Q_fieldQbig[i], bins=np.ceil(len(Q_fieldQbig[i])/2.), density=True)
      small = np.histogram(Q_fieldQsmall[i], bins=np.ceil(len(Q_fieldQsmall[i])/2.), density=True)

      ax = fig.add_subplot(5, 2, iof)
      plt.bar((big[1][:-1] + big[1][1:])/2., big[0], width=np.mean(np.diff(big[1])), color=colors[j], alpha=0.6, label=model_runs_longnames[j])
      plt.bar((small[1][:-1] + small[1][1:])/2., small[0], width=np.mean(np.diff(small[1])), edgecolor=colors[j], alpha=0.6, linewidth=3, label=model_runs_longnames[j], facecolor='none')
      ax.set_title(plotRivers[i], fontsize=16, fontweight='bold')
    except:
      pass
    iof += 1
    
    #ax.plot((big[1][:-1] + big[1][1:])/2., big[0], '-', color=colors[j], alpha=0.6, linewidth=3, label=model_runs_longnames[j])
    #ax.plot((small[1][:-1] + small[1][1:])/2., small[0], '-', color='r', alpha=0.6, linewidth=3, label=model_runs_longnames[j])
    #ax.plot(small[:,0], (small[:-1,1] + small[1:,-1])/2., '--', color=colors[j], alpha=0.6, linewidth=3, label=model_runs_longnames[j])

  fig.tight_layout()
  fig.subplots_adjust(left=.17, bottom=.08, right=None, top=None, wspace=None, hspace=.4)
  plt.savefig('field_model_Q_comparison_'+model_runs[j]+'.png')
  plt.savefig('../figures/field_model_Q_comparison_'+model_runs[j]+'.pdf')
  plt.close()
  
  # Mean discharge during high and low periods
  tmp_meanQ = []
  for i in range(len(Q_fieldQbig)):
    tmp_meanQ.append(np.nanmean(Q_fieldQbig[i]))
  MeanQ_fieldQbig.append(tmp_meanQ)
  tmp_meanQ = []
  for i in range(len(Q_fieldQsmall)):
    tmp_meanQ.append(np.nanmean(Q_fieldQsmall[i]))
  MeanQ_fieldQsmall.append(tmp_meanQ)

  # Now the stats
  D_values_i = []
  for ii in range(len(shortRivers)):
    D_values_i.append(stats.ks_2samp(Q_fieldQbig[ii], Q_fieldQsmall[ii])[0])
  D_values.append(D_values_i)
  """
  p_values_i = []
  for ii in range(len(shortRivers)):
    p_values_i.append(stats.ks_2samp(Q_fieldQbig[ii], Q_fieldQsmall[ii])[1])
  p_values.append(p_values_i)

  #chisquare_i = []
    #chisquare_i.append(stats.chisquare(Q_fieldQbig[ii], Q_fieldQsmall[ii])[1]) 
    #p_values_i.append(stats.chisquare(Q_fieldQbig[ii], Q_fieldQsmall[ii])[0])
  """
  
MeanQ_fieldQbig = np.array(MeanQ_fieldQbig)
MeanQ_fieldQsmall = np.array(MeanQ_fieldQsmall)

#p_values = np.array(p_values)
#np.sum(p_values > 0.05, axis=1)
#np.sum(p_values, axis=1)
# Maybe plot all rivers (x) against all models (line colors) on semilogy axis

D_values = np.array(D_values)
Dvm = np.mean(D_values[:,:-2], axis=1)

plt.figure(figsize=(8,12))
plt.subplot(211)
for j in range(len(Q)):
  diffnorm = (MeanQ_fieldQbig[j] - MeanQ_fieldQsmall[j])/(((MeanQ_fieldQsmall[j] + MeanQ_fieldQsmall[j]))/2.)
  #diffnorm[diffnorm<0] = 1E-2
  plt.semilogy(diffnorm, '-o', markersize=10, color=colors[j], alpha=0.6, linewidth=3, label=model_runs_longnames[j])
  plt.plot([8.5,9], [np.mean(diffnorm),np.mean(diffnorm)], '-', color=colors[j], alpha=0.6, linewidth=6) # mean plt.ylim((1E-1,1E2))
#plt.ylabel('Normalized difference of large- and small-flow modeled $Q$\n$(Q_l - Q_s) / (Q_l + Q_s)$', fontsize=20)
plt.ylabel('$(Q_h - Q_l) / (Q_h + Q_l)$', fontsize=20)
#plt.ylim((0.1, 10))
plt.xticks(range((len(srnu_spaces)+1)), list(srnu_spaces) + ['MEAN OF\nGLACIATED\nDRAINAGES'], rotation=45, ha='right', fontweight='bold')
#######
plt.subplot(212)
for j in range(len(Q)):
  plt.plot(D_values[j], '-o', markersize=10, color=colors[j], alpha=0.6, linewidth=3, label=model_runs_longnames[j])
  plt.plot([8.5,9], [Dvm[j],Dvm[j]], '-', color=colors[j], alpha=0.6, linewidth=6) # mean of glacial
plt.legend(loc='lower right')
plt.ylim((0, 1))
plt.xticks(range((len(srnu_spaces)+1)), list(srnu_spaces) + ['MEAN OF\nGLACIATED\nDRAINAGES'], rotation=45, ha='right', fontweight='bold')
plt.ylabel('Kolmogorov-Smirnov distance ($D$)', fontsize=20)
#plt.text(0.1, 0.02, 'High values indicate that modeled discharges\ndiffer significantly during data-derived times\nofhigh and low flow', horizontalalignment='left', verticalalignment='bottom', fontsize=14)
plt.tight_layout()
plt.savefig('../figures/inprogress/modelTiming_difference_and_goodness_of_fit.svg')
plt.show()

"""
#  plt.plot(MeanQ_fieldQbig[j]/np.mean(, '-v', markersize=10, color=colors[j], alpha=0.6, linewidth=3, label=model_runs_longnames[j])
#  plt.plot(MeanQ_fieldQsmall[j], '-^', markersize=10, color=colors[j], alpha=0.6, linewidth=3, label=model_runs_longnames[j])
"""

"""
Dvm = np.mean(D_values[:,:-2], axis=1)
plt.figure()
for j in range(len(Q)):
  #plt.semilogy(p_values[j], 'o', markersize=18, color=colors[j], label=model_runs_longnames[j])
  plt.plot(D_values[j], '-o', markersize=10, color=colors[j], alpha=0.6, linewidth=3, label=model_runs_longnames[j])
  plt.plot([8.5,9], [Dvm[j],Dvm[j]], '-', color=colors[j], alpha=0.6, linewidth=6) # mean of glacial
plt.legend(loc='lower right')
plt.ylim((0, 1))
plt.xticks(range((len(srnu_spaces)+1)), list(srnu_spaces) + ['MEAN OF\nGLACIATED\nDRAINAGES'], rotation=45, ha='right', fontweight='bold')
#plt.xticks([9], ['MEAN'], rotation=45, ha='right', fontweight='bold', fontsize=16)
plt.ylabel('Kolmogorov-Smirnov distance ($D$)', fontsize=16)
plt.text(0.1, 0.02, 'High values indicate that modeled discharges\ndiffer significantly during data-derived times\nofhigh and low flow', horizontalalignment='left', verticalalignment='bottom', fontsize=14)
plt.tight_layout()
#plt.savefig('../figures/inprogress/modelTiming_goodness_of_fit.svg')
plt.show()
"""

"""
plt.figure()
for j in range(len(Q)):
  #plt.semilogy(p_values[j], 'o', markersize=18, color=colors[j], label=model_runs_longnames[j])
  plt.semilogy(p_values[j], '-o', markersize=10, color=colors[j], alpha=0.6, linewidth=3, label=model_runs_longnames[j])
plt.legend(loc='upper right')
plt.ylim((1E-6, 1))
plt.xticks(range(9), srnu_spaces, rotation=45, ha='right', fontweight='bold')
plt.ylabel('p-value', fontsize=16)
plt.tight_layout()
plt.show()

plt.figure()
for j in range(len(Q)):
  plt.plot(p_values[j], '-', color=colors[j], alpha=0.6, linewidth=3, label=model_runs_longnames[j])
plt.legend(loc='lower left')
plt.show()  
"""

"""
    ax.plot(Qi['ages_numeric']/1000., Qi['Q_t'][rnu], color=colors[j], alpha=0.6, linewidth=3, label=model_runs_longnames[j])
    ax.set_title(plotRivers[i], fontsize=16, fontweight='bold')
    ax.set_xlim((0, 20))
    ax.set_ylim((0, ax.get_ylim()[-1]))
    # set y-lim based on model outputs
    if iof == 9:
      ax.legend(loc = 'center', bbox_to_anchor = (1.71, 0.5), fontsize=14)
    iof += 1
    if rnu == 'Mississippi':
      ax.set_ylim((0, 160000))
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
      ax.set_ylim((0, 300))
  fig.tight_layout()
  fig.subplots_adjust(left=.17, bottom=.08, right=None, top=None, wspace=None, hspace=.4)
  j += 1

plt.savefig('all_rivers.png')
plt.savefig('../figures/all_rivers.pdf')
plt.close()
"""
