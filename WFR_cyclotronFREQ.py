# WFR_cyclotronFREQ.py
#
# determine qb/m find the wave amplitude
# sort by L, MLT
#
# LKS January 2015
#
# imports
import numpy as np
from spacepy import pycdf
import glob
import os
import pandas as pd
import datetime
import spacepy.pybats.kyoto as spk
import itertools as itert
import matplotlib.dates as dates
from scipy.interpolate import interp1d
import pickle 
from dateutil.relativedelta import relativedelta 
import h5py
import spacepy.datamodel as dm
os.chdir('/Users/loisks/Documents/Functions')
import pickling as pickling
os.chdir('/Users/loisks/Desktop/ResearchProjects/PlasmaWaves/')
#
# parameters and starting conditions
dateStart='20130201'
dateEnd='20150401'
date1=dateStart
date2=dateEnd
dt0=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')
def diff_month(d1, d2):
    return (d1.year - d2.year)*12 + d1.month - d2.month
total_month=diff_month(dt2, dt0)
nLbins=11
nmlt_bins=48
harmonic=16
sat=['A','B']
satellite=['a','b']
#dir=['EMFISIS_WFR_A', 'EMFISIS_WFR_B']
#dir2=['EMPHEMERIS_A', 'EMPHEMERIS_B']

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
name=['rbsp-a', 'rbsp-b']
name2=['rbspa', 'rbspb']
# pickle function


date=dateStart
endDt=datetime.datetime.strptime(dateEnd,'%Y%m%d')
DT=datetime.datetime.strptime(date, '%Y%m%d')


q=1.6*1e-19
m=1.67*1e-27
for iSat in range(1,len(sat)):
  E_Spectra_sorted=[[[] for i in range(nLbins)] for j in range(nmlt_bins)]
  B_Spectra_sorted=[[[] for i in range(nLbins)] for j in range(nmlt_bins)]
  DT=datetime.datetime.strptime(dateStart, '%Y%m%d')
  monthCur=DT.month
  while DT != endDt:
    try:
     date=datetime.datetime.strftime(DT, '%Y%m%d')
     DT=datetime.datetime.strptime(date, '%Y%m%d')
     yrCur=DT.year
     # first get the magnetic field first
     os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EMFISIS_'+sat[iSat])
     g=glob.glob('rbsp-'+satellite[iSat]+'_magnetometer_4sec-gse_emfisis-L3_'+date+'*')
     pyfem=pycdf.CDF(g[0])
     Bepoch=pd.DatetimeIndex(pyfem['Epoch'][...])
     # just make it cyclotron frequency to make it easier 
     Bmag=harmonic*(q*np.array(pyfem['Magnitude'][...])*1e-9)/(2.0*np.pi*m) 
     # resample both to 1 minute periods
     rt = pd.period_range(date,periods=1440, freq='T').to_timestamp()
     rng=rt[::1]
     BmagDF=pd.DataFrame({'Bmag':Bmag}, index=Bepoch)
     BmagDF=BmagDF['Bmag'].resample('1min', how='median').reindex(index=rng,fill_value=np.nan)
     # find the nearest frequency in the EMFISIS data
     nearestFreq=np.zeros(len(BmagDF))
     for iCyc in range(len(BmagDF)):
         # this gives the index
       nearestFreq[iCyc]=min(range(len(FREQUENCIES)), key=lambda i: abs(FREQUENCIES[i]-BmagDF[iCyc]))  
     
     #
     # got that, now onto wave amplitudes
     os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EMFISIS_WFR_'+sat[iSat])
     g=glob.glob('rbsp-'+satellite[iSat]+'_WFR-spectral-matrix-diagonal_emfisis-L2_'+date+'*')
     pyfW=pycdf.CDF(g[0])
     Wepoch=pd.DatetimeIndex(pyfW['Epoch'][...])
     #BSpectra=np.zeros(len(BmagDF))
     ESpectra=np.zeros(len(BmagDF))
     for iTime in range(len(BmagDF)):
        Freq=nearestFreq[iTime]
        if Freq >1: # these are at perigee
          #BuBu=np.swapaxes(pyfW['BuBu'][...],1,0)[Freq]
          #BvBv=np.swapaxes(pyfW['BvBv'][...],1,0)[Freq]
          #BwBw=np.swapaxes(pyfW['BwBw'][...],1,0)[Freq]
          EuEu=np.swapaxes(pyfW['EuEu'][...],1,0)[Freq]
          EvEv=np.swapaxes(pyfW['EvEv'][...],1,0)[Freq]
          EwEw=np.swapaxes(pyfW['EwEw'][...],1,0)[Freq]
          #BT=BuBu+BvBv+BwBw
          ET=EuEu+EvEv+EwEw
          # now index and resample this shit
          #BTDT=pd.DataFrame({'B':BT}, index=Wepoch)
          ETDT=pd.DataFrame({'E':ET},index=Wepoch)
          #BTDT=BTDT['B'].resample('1min', how='median').reindex(index=rng,fill_value=np.nan)
          ETDT=ETDT['E'].resample('1min', how='median').reindex(index=rng,fill_value=np.nan)
          # now take the appropriate value
          #BSpectra[iTime]=BTDT[iTime]
          ESpectra[iTime]=ETDT[iTime]
        else: # just nan the apogee times
          #BSpectra[iTime]=np.nan
          ESpectra[iTime]=np.nan
      # now we have amplitudes at corresponding frequencies
      # now get HOPE ephemeris
     os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EMPHEMERIS_'+sat[iSat])
     hope_emphem=glob.glob('rbsp'+satellite[iSat]+'_def_MagEphem_OP77Q_'+date+'*.txt')[0]
     pyf2=dm.readJSONheadedASCII(hope_emphem)
     L_emphem=np.nanmedian(pyf2['L'], axis=1) # L from emphemeris
     L_emphem[L_emphem<0]=np.nan
     
     MLT_emphem=pyf2['CDMAG_MLT'] # MLT from ephemeris
     MLT_emphem[MLT_emphem<0]=np.nan
     Kp_emphem=pyf2['Kp'] # Kp index from ephemeris file
     epoch_ephem=pd.DatetimeIndex(pyf2['DateTime'])
     #
     # great, now resample this stuff
     LEMP=pd.DataFrame({'L':L_emphem},index=epoch_ephem)
     MLTEMP=pd.DataFrame({'MLT':MLT_emphem},index=epoch_ephem)
     KPEMP=pd.DataFrame({'KP':Kp_emphem},index=epoch_ephem)
     L=LEMP['L'].resample('1min',how='median').reindex(index=rng,fill_value=np.nan)
     MLT=MLTEMP['MLT'].resample('1min',how='median').reindex(index=rng,fill_value=np.nan)
     KP=KPEMP['KP'].resample('1min',how='median').reindex(index=rng,fill_value=np.nan)
     #
     # get rid of KP < 3 data
     KPnan=np.where(np.array(KP) > 3)[0]
     #BSpectra[KPnan]=np.nan
     ESpectra[KPnan]=np.nan
     
     L=np.array(L); MLT=np.array(MLT)
     #
     # Alright now, let's sort everything
     for imlt in range(nmlt_bins):
            for iL in range(nLbins):
                L_indexes=np.where((L >= 1.5+0.25*iL) & (L < 1.5+0.25*(iL+1)))[0]
                MLT_indexes=np.where((MLT >= 0.5*imlt) & (MLT < 0.5*(imlt+1)))[0]
                overlap=list(set(L_indexes) & set(MLT_indexes))
                E_Spectra_sorted[imlt][iL]=E_Spectra_sorted[imlt][iL]+list(ESpectra[overlap])
                #B_Spectra_sorted[imlt][iL]=B_Spectra_sorted[imlt][iL]+list(BSpectra[overlap])
     # save by month
     print(date)
    except:
      print('ERROR ' + str(date))
    yrCur=DT.year
    DT=DT+datetime.timedelta(days=1)
    if DT.month != monthCur:
      # save this file
        os.chdir('/Users/loisks/Documents/Functions/')
        import pickling
        os.chdir('/Users/loisks/Desktop/ResearchProjects/PlasmaWaves/')
        pickling.hdf5_data_save(E_Spectra_sorted, 'ESorted_Cyc_'+str(monthCur)+'_'+str(yrCur)+'_'+ name2[iSat], 'WFR_FREQ_240',  nLbins, nmlt_bins)
        monthCur=DT.month
        print(str(yrCur)+'_'+str(monthCur))
        E_Spectra_sorted=[[[] for i in range(nLbins)] for j in range(nmlt_bins)]
       
