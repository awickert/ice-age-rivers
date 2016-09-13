#! /usr/bin/python

from glob import glob
import os
import numpy as np

# Rename the files to be 0-padded
# run in proper directory
#os.system("rename -v 's/(.*)(\d.*)/${1}0${2}/' *.?.mn") # -n 

# ICE-6G
indir = '/home/awickert/RUNS_ICE6G/60VM5/'
outdir = '/home/awickert/RUNS_ICE6G/60VM5out/'
ages = np.array( [25, 24, 23, 22, 21] + list(np.arange(0,21,.5))[::-1] ) * 1000.

# ICE-5G
#indir = '/home/awickert/RUNS_ICE5G/outaw090CVM2VM2/'
#outdir = '/home/awickert/RUNS_ICE5G/outaw090CVM2VM2out/'
#ages = np.genfromtxt('/home/awickert/Dropbox/all_codes_glacial/ContinentalDrainage/input/step_time')
#ages = ages[:,1]
#ages = ages[ages <= 26] # let's just take 26 ka to present
#ages *= 1000

# ICE-3G
#indir = '/home/awickert/RUNS_ICE3G/outaw12512/'
#outdir = '/home/awickert/RUNS_ICE3G/outaw12512out/'
#ages = np.genfromtxt('/home/awickert/IceModelAges/ages_lgm_to_present_3G_reversed')[::-1]

# ANU
#indir = '/home/awickert/RUNS_ANU/outaw72p35/'
#outdir = '/home/awickert/RUNS_ANU/outaw72p35out/'
#ages = np.genfromtxt('/home/awickert/Dropbox/all_codes_glacial/ContinentalDrainage/input/ages_lgm_to_present_ANU')
#ages *= 1000

# G12
#indir = '/home/awickert/RUNS_G12/outaw090CVM2VM2_gregoire/'
#outdir = '/home/awickert/RUNS_G12/outaw090CVM2VM2_gregoireout/'
#ages = np.arange(20900, 499, -100) # Model starts at 21,000, but JXM always rmv's ts #1 Hence 205 ts here but 206 in LG's output
#os.system("rename -v 's/(.*)(\d{2}.*)/${1}0${2}/' *.??.mn") # -n
#files = sorted(glob(indir+'topo_3.???.mn'))

files = sorted(glob(indir+'topo_3.??.mn'))

# STRING FORMATTING HAPPENS IN BASH
#ages = []
#for i in range(len(ages_numeric)):
#  ages.append('%06d' %ages_numeric[i])

# Modify files to fit ages
files=files[-len(ages):]

for i in range(len(files)):
    print ages[i], os.path.split(files[i])[-1]
    os.system('./grid_mn_ns30as.x << EOF\n' \
              +files[i]+'\n' \
              +outdir+'topo'+'%06d' %ages[i]+'.txt\n' \
              +'EOF')
    print ""

