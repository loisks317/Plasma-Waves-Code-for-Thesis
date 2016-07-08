# SimpleModel2.0.py
#
# Schmitt 1977 model on cyclotron resonant heating
# Also use culmulative heating from broadband frequencies between 45 Hz and 100 Hz. 
#
# LKS, November 2015, updated December 2015 to focus on MLT instead of MLAT
# updated February 2016 in Svalbard to account for lower frequency
# broadband effects and linear polarization (1/2) 
#
# imports
import numpy as np
from spacepy import pycdf
import glob
from numpy import ma
import os
import pandas as pd
import datetime
import spacepy.pybats.kyoto as spk
import itertools as itert
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter 
import matplotlib.dates as dates
import pickle 
from dateutil.relativedelta import relativedelta 
import h5py
import spacepy.datamodel as dm
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
import numpy as np
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import (FixedLocator, MaxNLocator,
                                                 DictFormatter)
import matplotlib.pyplot as plt
os.chdir('/Users/loisks/Desktop/Functions/')
import pickling as pickling
import plots
os.chdir('/Users/loisks/Documents/ResearchProjects/PlasmaWaves/')
#
FREQUENCIES=[ 2.13599992e+00,   4.27199984e+00,   6.40899992e+00,
          8.54500008e+00,   1.06800003e+01,   1.28199997e+01,
          1.49499998e+01,   1.70900002e+01,   1.92299995e+01,
          2.13600006e+01,   2.35000000e+01,   2.56299992e+01,
          2.77700005e+01,   3.09799995e+01,   3.52500000e+01,
          3.95200005e+01,   4.37900009e+01,   4.91300011e+01,
          5.55400009e+01,   6.30200005e+01,   7.15599976e+01,
          8.01100006e+01,   8.97200012e+01,   1.00400002e+02,
          1.12199997e+02,   1.26000000e+02,   1.42100006e+02,
          1.59199997e+02,   1.78399994e+02,   1.99699997e+02,
          2.24300003e+02,   2.52100006e+02,   2.82000000e+02,
          3.16200012e+02,   3.54600006e+02,   3.98399994e+02,
          4.47500000e+02,   5.02100006e+02,   5.62900024e+02,
          6.31299988e+02,   7.09200012e+02,   7.95799988e+02,
          8.91900024e+02,   1.00100000e+03,   1.12400000e+03,
          1.26100000e+03,   1.41500000e+03,   1.58700000e+03,
          1.78100000e+03,   1.99800000e+03,   2.24300000e+03,
          2.51600000e+03,   2.82300000e+03,   3.16800000e+03,
          3.55500000e+03,   3.98800000e+03,   4.47400000e+03,
          5.02000000e+03,   5.63300000e+03,   6.32000000e+03,
          7.09100000e+03,   7.95600000e+03,   8.92700000e+03,
          1.00200000e+04,   1.12400000e+04]
def factorial(n):
    if n == 0:
        return 1
    else:
        return n * factorial(n-1)
# to get the cool floating axes
def setup_axes2(fig, rect):
          """
          With custom locator and formatter.
          Note that the extreme values are swapped.
          """
          tr = PolarAxes.PolarTransform()
      
          pi = np.pi
          angle_ticks = [((-30/180.)*np.pi, r"-30$^{o}$"),
                         ((-20/180.)*np.pi, r"-20$^{o}$"),
                         (-1*(10/180.)*np.pi, r"-10$^{o}$"),
                         (0, r"0$^{0}$"),
                         ((10/180.)*np.pi, r"10$^{o}$"),
                         ((20/180.)*np.pi, r"20$^{o}$"),
                         ]
          grid_locator1 = FixedLocator([v for v, s in angle_ticks])
          tick_formatter1 = DictFormatter(dict(angle_ticks))
      
          grid_locator2 = MaxNLocator(3)
      
          grid_helper = floating_axes.GridHelperCurveLinear(
              #tr, extremes=(.5*np.pi, -0.25*np.pi, 4, 1.5),
              tr,extremes=((1/6.)*np.pi, -(1/6.)*np.pi, 4,1.5),
              grid_locator1=grid_locator1,
              grid_locator2=grid_locator2,
              tick_formatter1=tick_formatter1,
              tick_formatter2=None)
      
          ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
          fig.add_subplot(ax1)
      
          # create a parasite axes whose transData in RA, cz
          aux_ax = ax1.get_aux_axes(tr)
      
          aux_ax.patch = ax1.patch  # for aux_ax to have a clip path as in ax
          ax1.patch.zorder = 0.9  # but this has a side effect that the patch is
          # drawn twice, and possibly over some other
          # artists. So, we decrease the zorder a bit to
          # prevent this.
      
          return ax1, aux_ax

nLbins=11
nmltbins=48
nMLATbins=31
LBins=np.linspace(1.5, 4, nLbins)+0.125
MLATbins=np.linspace(-20, 10, nMLATbins)+0.5
MLATbinsALL=np.linspace(-30,30,61)
HzIndex=range(15,24) # 39 Hz to 100 Hz, still N > 3
#HzAll=[100,145,200,250, 300, 600,1000]
# get the right directory
for iHz in HzIndex:
      CurHz=FREQUENCIES[iHz]
      os.chdir('/Users/loisks/Documents/ResearchProjects/PlasmaWaves/Sorted_WFR_Polarization')
      totalE=[[[] for i in range(nLbins)] for j in range(nmltbins)]
      totalPol=[[[] for i in range(nLbins)] for j in range(nmltbins)]
      medE=[[ 0 for i in range(nLbins)] for j in range(nmltbins)]
      medPol=[[ 0 for i in range(nLbins)] for j in range(nmltbins)]
      #
      # need to get the files and then sort by polarization
      EFiles=glob.glob('*_polarization_E_WFR_Spectra_freq='+str(FREQUENCIES[iHz])+'_rbsp'+'*')
      Polfiles=glob.glob('*_polarization_freq='+str(FREQUENCIES[iHz])+'_*')
      for iFile in range(len(EFiles)):
              data=pickling.hdf5_data_open(EFiles[iFile],nLbins, nmltbins)
              dP=pickling.hdf5_data_open(Polfiles[iFile],nLbins, nmltbins)
              print(str(iFile) + ' / ' +str(len(EFiles)))
              for iMLT in range(nmltbins):
                  for iL in range(nLbins):
                      try:
                          totalE[iMLT][iL]+=list(np.array(data[iMLT][iL]))
                          totalPol[iMLT][iL]+=list(np.array(dP[iMLT][iL]))
                      except(OSError):
                          totalE[iMLT][iL]+=[np.nan]
                          totalPol[iMLT][iL]+=[np.nan]
      # now get the median
      for iMLT in range(nmltbins):
              for iL in range(nLbins):
                  tempE=np.array(totalE[iMLT][iL])[np.array(totalPol[iMLT][iL])<=0.2]
                  medE[iMLT][iL]=np.nanmedian(tempE)/2.0 # the 2.0 accounts for
                # the linear polarization
      os.chdir('..')
      os.chdir('Sorted_BMag')
      Bmed=[[ 0 for i in range(nLbins)] for j in range(nmltbins)]
      files=glob.glob('*_MLT_*')
      tempB1=np.array(pickling.hdf5_data_open(files[0], nLbins, nmltbins))
      tempB2=np.array(pickling.hdf5_data_open(files[1], nLbins, nmltbins))
      #
      for iMLT in range(nmltbins):
          for iL in range(nLbins):
              try:
                  temp=np.array(list(tempB1[iMLT][iL])+list(tempB2[iMLT][iL]))*1e-9
                  temp[temp<0]=np.nan
                  temp[temp>1e-4]=np.nan
                  Bmed[iMLT][iL]=np.nanmedian(temp)
              except(OSError):
                  Bmed[iMLT][iL]=np.nan
      os.chdir('..')
      # now plot this
      modMedE=[[ np.nan for i in range(nLbins)] for j in range(nmltbins)]
      for ii in range(nmltbins):
         for jj in range(nLbins):
             modMedE[ii][jj]=medE[ii][jj]
            

      # use the EMFISIS B data 
      ER=6378000 # in meters
      mi=1.67*1e-27
      q=1.6*1e-19
      kperpArr=[1e-2,5*1e-3,1e-3, 1e-4]
      DegTrap=[89.9,89,88,80,70,60,50,40,30,20,10]
      for iK in range(len(kperpArr)):
       for iDeg in range(len(DegTrap)):
          kperp=kperpArr[iK]
          PA=DegTrap[iDeg]
          vTot=9788.18*np.sqrt(2) # for a 1 eV particle
          vperp=np.sqrt((np.sin(PA*np.pi/180.)**2)*(vTot**2))
          vpar=np.sqrt((vTot**2) - (vperp**2))
          sqrtE=np.sqrt(np.array(medE)*CurHz)# need to multiply by  Hz in here
          # LKS unit workout, Nov 2015
          LBinsMod=np.array(np.linspace(1.5, 4, nLbins))+0.125
          MLTbinsArr=np.linspace(0.25, 24.25,nmltbins+1)
          n=np.ones((nmltbins,nLbins))*np.nan
          t=np.ones((nmltbins,nLbins))*np.nan
          BounceTime=np.ones((nmltbins,nLbins))*np.nan
          tB=np.ones((nmltbins,nLbins))*np.nan
          os.chdir('/Users/loisks/Desktop/Functions/')
          # this assumes a centered dipole where mag lat = latitude. This isn't
          # quite right
          #
          # check out subharmonic resonances in AM
          for iL in range(nLbins):
              for iMLT in range(nmltbins):
                  nT=int(2.0*np.pi*CurHz*mi/(q*Bmed[iMLT][iL]))
                  if nT < 1:
                      # pass
                      print(q*Bmed[iMLT][iL]/mi)
                  else:
                     n[iMLT][iL]=int(nT)
                     gamma=kperp*sqrtE[iMLT][iL]/Bmed[iMLT][iL]
                     Ro=kperp*vperp*mi/(q*Bmed[iMLT][iL])
                     t[iMLT][iL]=(1.0/gamma)*(np.power(2,(n[iMLT][iL]-1)))*factorial(n[iMLT][iL]-1)/((n[iMLT][iL]-2)*(np.power(Ro,(n[iMLT][iL]-2))))
                     BounceTime[iMLT][iL]=(LBins[iL]*ER/vpar)*(3.7-1.6*np.sin(np.pi*PA/180.)) # used PA = 30
                     tB[iMLT][iL]=t[iMLT][iL]/BounceTime[iMLT][iL]
          # need to check those values
          os.chdir('/Users/loisks/Documents/ResearchProjects/PlasmaWaves/')
          #
          # now plot log of time for these
          modMedT=[[ np.nan for i in range(nLbins)] for j in range(nmltbins)]
          modMedTb=[[ np.nan for i in range(nLbins)] for j in range(nmltbins)]
          for ii in range(nmltbins):
             for jj in range(nLbins):
                 modMedT[ii][jj]=t[ii][jj]
                 modMedTb[ii][jj]=tB[ii][jj]

          plots.fancy_plot(modMedTb, nmltbins, nLbins, -1, 2, 'Resonance Time/Bounce Period', 'SimpleModel2.0_Results', 'low', str(CurHz)+'Hz_En=1_K='+str(iK)+'_PA='+str(PA), 'SimpleModel2.0_Results', 'log', 'PiYG')
          plots.fancy_plot(Bmed, nmltbins, nLbins, -4, -7, '|B| [T]', 'SimpleModel_Results', 'low', 'Bmag', 'SimpleModel_Results', 'log', 'viridis') 
          


      
