# waveAmplitudesMLAT.py
#
# get the wave amplitudes at different magnetic latitudes
# to see how equatorial noise attenuates off the equator
#
# LKS, November 2015
#
# 100 Hz = 23, 316 Hz = 33, 631 Hz = 39,  1000 Hz = 43
# 150 Hz = 26, 200 Hz = 29, 250 Hz = 31
# 50 Hz = 17
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
# parameters and starting conditions
dateStart='20130201'
dateEnd='20140901'
date1=dateStart
date2=dateEnd
dt0=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')
def diff_month(d1, d2):
    return (d1.year - d2.year)*12 + d1.month - d2.month
total_month=diff_month(dt2, dt0)
nLbins=11
nMLATbins=31
Lbins=np.linspace(1.5, 4, nLbins)
MLATbins=np.linspace(-20, 10, nMLATbins)
sat=['A','B']
satellite=['a','b']
date=dateStart
endDt=datetime.datetime.strptime(dateEnd,'%Y%m%d')
DT=datetime.datetime.strptime(date, '%Y%m%d')
monthCur=DT.month
freqIndex=[17]
freqLabel=[50]
#freqIndex=[23, 26,29,31, 33]
#freqLabel=[100, 150, 200, 250, 300]

# keep it simple, only get the E Spectra
#
#
for iFreq in range(len(freqIndex)):
  iFr=freqIndex[iFreq]  
  for iSat in range(len(sat)):
    E_Spectra_sorted=[[[] for i in range(nLbins)] for j in range(nMLATbins)]
    DT=datetime.datetime.strptime(dateStart, '%Y%m%d') 
    while DT != endDt:
        try:
            date=datetime.datetime.strftime(DT, '%Y%m%d')
            DT=datetime.datetime.strptime(date, '%Y%m%d')
            rt = pd.period_range(date,periods=1440, freq='T').to_timestamp()
            # got that, now onto wave amplitudes
            os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EMFISIS_WFR_'+sat[iSat])
            g=glob.glob('rbsp-'+satellite[iSat]+'_WFR-spectral-matrix-diagonal_emfisis-L2_'+date+'*')
            pyfW=pycdf.CDF(g[0])
            Wepoch=pd.DatetimeIndex(pyfW['Epoch'][...])
            EUEU=pd.DataFrame({'EU':np.swapaxes(pyfW['EuEu'][...],1,0)[iFr]}, index=Wepoch)
            EVEV=pd.DataFrame({'EV':np.swapaxes(pyfW['EvEv'][...],1,0)[iFr]},index=Wepoch)
            EWEW=pd.DataFrame({'EW':np.swapaxes(pyfW['EwEw'][...],1,0)[iFr]},index=Wepoch)
            EU=np.array(EUEU['EU'].resample('1min', how='median').reindex(index=rt,fill_value=np.nan))
            EV=np.array(EVEV['EV'].resample('1min',how='median').reindex(index=rt,fill_value=np.nan))
            EW=np.array(EWEW['EW'].resample('1min',how='median').reindex(index=rt,fill_value=np.nan))
            ESpectra=np.zeros(len(rt))
            for iTime in range(len(rt)):
               ESpectra[iTime]=EU[iTime]+EV[iTime]+EW[iTime]
        
            os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EMPHEMERIS_'+sat[iSat])
            hope_emphem=glob.glob('rbsp'+satellite[iSat]+'_def_MagEphem_OP77Q_'+date+'*.txt')[0]
            pyf2=dm.readJSONheadedASCII(hope_emphem)
            L_emphem=np.nanmedian(pyf2['L'], axis=1) # L from emphemeris
            L_emphem[L_emphem<0]=np.nan
            MLT_emphem=pyf2['CDMAG_MLT']
            MLT_emphem[MLT_emphem<0]=np.nan
            MLAT_emphem=pyf2['CDMAG_MLAT'] # MLT from ephemeris
            MLAT_emphem[MLAT_emphem<0]=np.nan
            Kp_emphem=pyf2['Kp'] # Kp index from ephemeris file
            epoch_ephem=pd.DatetimeIndex(pyf2['DateTime'])
            #
            # great, now resample this stuff
            LEMP=pd.DataFrame({'L':L_emphem},index=epoch_ephem)
            MLATEMP=pd.DataFrame({'MLAT':MLAT_emphem},index=epoch_ephem)
            MLTEMP=pd.DataFrame({'MLT':MLT_emphem}, index=epoch_ephem)
            KPEMP=pd.DataFrame({'KP':Kp_emphem},index=epoch_ephem)
            L=LEMP['L'].resample('1min',how='median').reindex(index=rt,fill_value=np.nan)
            MLAT=MLATEMP['MLAT'].resample('1min',how='median').reindex(index=rt,fill_value=np.nan)
            MLT=MLTEMP['MLT'].resample('1min',how='median').reindex(index=rt, fill_value=np.nan)
            KP=KPEMP['KP'].resample('1min',how='median').reindex(index=rt,fill_value=np.nan)
            #
            # get rid of KP < 3 data
            KPnan=np.where(np.array(KP) > 3)[0]
            ESpectra[KPnan]=np.nan
     
            L=np.array(L); MLAT=np.array(MLAT); MLT=np.array(MLT)
            #
            # Alright now, let's sort everything
            for imlat in range(nMLATbins):
                for iL in range(nLbins):
                        L_indexes=np.where((L >= 1.5+0.25*iL) & (L < 1.5+0.25*(iL+1)))[0]
                        MLAT_indexes=np.where((MLAT >= 0.5*imlat) & (MLAT < 0.5*(imlat+1)))[0]
                        # only take dayside
                        MLT_indexes=np.where((MLT >= 6) & (MLT <=18))[0]
                        overlap=list(set(L_indexes) & set(MLAT_indexes) & set(MLT_indexes))
                        E_Spectra_sorted[imlat][iL]=E_Spectra_sorted[imlat][iL]+list(ESpectra[overlap])
               
        except:
            print 'ERROR ' + str(date)
        yrCur=DT.year
        DT=DT+datetime.timedelta(days=1)
        if DT.month != monthCur:
            # save this file
            os.chdir('/Users/loisks/Documents/Functions/')
            import pickling
            os.chdir('/Users/loisks/Desktop/ResearchProjects/PlasmaWaves/')
            pickling.hdf5_data_save(E_Spectra_sorted,'ESpectra_sat=rbsp'+str(satellite[iSat])+'_'+str(yrCur)+'_'+str(monthCur) , 'WaveAmplitudes_'+str(freqLabel[iFreq])+'_Hz', nMLATbins, nLbins)
    
            monthCur=DT.month
            print monthCur
            E_Spectra_sorted=[[[] for i in range(nLbins)] for j in range(nMLATbins)]
        
