# emf_l4_sort.py
#
# the goal here is to get the fraction of left hand waves between L = 2 and L = 3
# right now, I'm going to just the EMFISIS team L4 files and not bin by MLT
# this should just be a table.
#
# LKS, March 2016 - Tromso, Norway.
#
# imports
import numpy as np
from spacepy import pycdf
import glob
import pandas as pd
import os
import datetime
#
# counters
LeftHand=[ 0 for i in range(65)]
RightHand=[0 for i in range(65)]
Linear=[0 for i in range(65)]
DayCount=0
#
sats=['A', 'B']
for isat in range(2):
    # change to where EMFISIS data is
    os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EMFISIS_L4_'+sats[isat])
    globs=glob.glob('*.cdf')
    #
    # just get all the files and go from there
    for ifile in globs:
      try:
        pyf=pycdf.CDF(ifile)
        L=np.array(pyf['L'][...])
        Ell=np.array(pyf['ellsvd'][...])
        # nan the silly things and then remove
        nEll= Ell[np.where((L <= 3) & (L>= 2))[0]]
        # for each frequency
        for ifreq in range(65):
            fnEll=nEll[:,ifreq]
            LeftHand[ifreq]+=len(np.where(fnEll<-0.2)[0])
            RightHand[ifreq]+=len(np.where(fnEll>0.2)[0])
            Linear[ifreq]+=len(np.where((fnEll>-0.2) & (fnEll<0.2))[0])
            if ifreq==30:
                print(LeftHand[ifreq])
                print(RightHand[ifreq])
                print(Linear[ifreq])                
        DayCount+=1
    
      except:
          print('bad day')
#
# pickle
import pickle
os.chdir('/Users/loisks/Desktop/ResearchProjects/PlasmaWaves/')
pickle.dump(LeftHand, open('leftHandPol.p', 'wb'),protocol=2)
pickle.dump(RightHand, open('rightHandPol.p', 'wb'),protocol=2)
pickle.dump(Linear, open('linearPol.p', 'wb'),protocol=2)
