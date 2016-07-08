# broadbandsortFreqBand.py
#
# this script specifically looks at broadband cyclotron resonant frequencies
# for H+ vs He+ 
# with polarization of < 0.3
# keep MLT and L shell data
# resample to 1 second?
#
# LKS, SANSA 2016 March
#
#
#
# imports
import numpy as np
from spacepy import pycdf
import glob
import pandas as pd
import os
import datetime
import spacepy.pybats.kyoto as spk
import itertools as itert
import matplotlib.dates as dates
from scipy.interpolate import interp1d
import pickle 
from dateutil.relativedelta import relativedelta 
import h5py
import spacepy.datamodel as dm
import generateMatrices as gM
os.chdir('/Users/loisks/Documents/Functions')
import pickling as pickling
os.chdir('/Users/loisks/Desktop/ResearchProjects/PlasmaWaves/')
#
# parameters and starting conditions
# need to check to make sure we have enough EMFISIS data tomorrow and
# start running this script
date1='20130201'
date2='20150401'
dt0=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')
def diff_month(d1, d2):
    return (d1.year - d2.year)*12 + d1.month - d2.month
total_month=diff_month(dt2, dt0)
Lbins=11
mltbins=48
dir=['EMFISIS_WFR_A', 'EMFISIS_WFR_B']
sat=['A','B']
lsat=['a','b']
#
# need frequency indexes for 10 Hz and 100 Hz
index10=1
index100=24

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

#
# load in the Kp
KpArr=pickle.load(open('KpFeb2013_Apr2015.p', 'rb'))
#
name=['rbsp-a', 'rbsp-b']
name2=['rbspa', 'rbspb']
# pickle function
#
def save_data_pickle(data, EN_LABEL, data_label) :
#
    import pickle
    import os
    with open(data_label + '_' +EN_LABEL, 'w') as f:
        pickle.dump(data, f)
    os.chdir('..')
#
# begin loop 
#for ifreq in range(len(FREQUENCIES)):
# frequencies from 70 to 1000 Hz
for idir in range(2):# CHANGE FOR BOTH A AND B
   # start overarching loop through days
    dt1=dt0
    E_Spectra_sorted=[[ [] for x in range(Lbins)] for x in range(mltbins)]
    Polarization_sorted=[[ [[] for z in range(index100-index10)] for x in range(Lbins)] for x in range(mltbins)]
    month_cur=dt1.month
    while dt1 < dt2:

      curDate=str(dt1.month)+'_'+str(dt1.year)
      try:        
        date1=datetime.datetime.strftime(dt1, '%Y%m%d')
        print(date1)
        # to wake up the nfs server I hope
        os.chdir('/Users/loisks/Desktop/liemohn10/loisks/'+dir[idir])
        emf_file=name[idir]+'_WFR-spectral-matrix_emfisis-L2_'+date1+'*.cdf'
        gemf=glob.glob(emf_file)
        pyf=pycdf.CDF(gemf[0])
        os.chdir('..')
        # get the polarization for each frequency in broad band
        polarizationH= [ [] for i in range(index100-index10)]
        E_SpectraH= [[] for i in range(index100-index10) ]
        polarizationHe= [ [] for i in range(index100-index10)]
        E_SpectraHe= [[] for i in range(index100-index10) ]
        for ifreq in range(index100-index10):
            ii = ifreq+index10
            polarization[ifreq]=gM.mag_SVD(pyf, ii)
            Eu_Spectra=np.array(np.swapaxes(pyf['EuEu'][...],1 , 0)[ii]) # spectra, 65 frequncies 
            Ev_Spectra=np.array(np.swapaxes(pyf['EvEv'][...],1 , 0)[ii])
            Ew_Spectra=np.array(np.swapaxes(pyf['EwEw'][...],1 , 0)[ii])
            E_Spectra[ifreq]=0.5*np.sqrt(np.array((Eu_Spectra+Ev_Spectra+Ew_Spectra)))#*FREQUENCIES[ii]))
          
        #
        # get the emphermeris file
        os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EMPHEMERIS_'+sat[idir])
        hope_emphem=glob.glob('rbsp'+lsat[idir]+'_def_MagEphem_OP77Q_'+date1+'*.txt')[0]
        pyf2=dm.readJSONheadedASCII(hope_emphem)
        L_emphem=np.nanmedian(pyf2['L'], axis=1) # L from emphemeris
        L_emphem[L_emphem<0]=np.nan
     
        MLT_emphem=pyf2['CDMAG_MLT'] # MLT from ephemeris
        MLT_emphem[MLT_emphem<0]=np.nan     
        epoch_ephem=epochL4=pd.DatetimeIndex(pyf2['DateTime'])
        #
        # great, now resample this stuff
        LEMP=pd.DataFrame({'L':L_emphem},index=epoch_ephem)
        MLTEMP=pd.DataFrame({'MLT':MLT_emphem},index=epoch_ephem)
          
        # get the coordinates from the EMFISIS L4 files
        # time array is 14400, reduce by factor of 10
        Kp=np.array(KpArr[date1])

# Spectra as time x frequency array
        #rng= pd.date_range(date1, periods=28, freq='30S')
        rng2=pd.date_range(date1, periods=1440, freq='1T')
        Kp=np.array(pd.DataFrame({'Kp':Kp}, index=rng2))[:,0]
       
        L=np.array(LEMP['L'].resample('1T',how='median').reindex(index=rng2,fill_value=np.nan))
        MLT=np.array(MLTEMP['MLT'].resample('1T',how='median').reindex(index=rng2,fill_value=np.nan))
        #Kp=np.array(pdKp['Kp'].resample('30S',how='median').reindex(index=rng,fill_value=np.nan))

        #
        # now get the EMFISIS spectral matrices
        epoch=pyf['Epoch'][...] # time
        # interpolate the spectra to  grid like the emphemeris data
        epochD=pd.DatetimeIndex(epoch)
        #
        # set up temp arrays for resampling
        polarizationT= [ [] for i in range(index100-index10)]
        E_SpectraT= [[] for i in range(index100-index10) ]
        ESP=np.zeros(len(rng2))
        for ii in range(index100-index10):
            ESpec=pd.DataFrame({'ESpec':E_Spectra[ii]},index=epochD)
            E_SpectraT[ii]=np.array(ESpec['ESpec'].resample('1T',how='median').reindex(index=rng2,fill_value=np.nan))
            pdPOL=pd.DataFrame({'Pol':polarization[ii]}, index=epochD)
            polarizationT[ii]=np.array(pdPOL['Pol'].resample('1T',how='median').reindex(index=rng2,fill_value=np.nan))
            E_SpectraT[ii][np.array(polarizationT[ii]) > 0.3] = 0
            E_SpectraT[ii][np.isnan(E_SpectraT[ii])]= 0 # just make nans zeros
            ESP+=E_SpectraT[ii] # ESP is the broadband
        
        # now in the loop 
        os.chdir('/Users/loisks/Desktop/ResearchProjects/PlasmaWaves/')
        ESP[Kp>3]=0
        #
        # screen for ellipticity and polarization
        # linearly polarized waves = thsvd -> 90 and ellsvd -> 0
# now go through the fluxes
        for imlt in range(mltbins):
            for iL in range(Lbins):
                    L_indexes=np.where((L >= 1.5+0.25*iL) & (L < 1.5+0.25*(iL+1)))[0]
                    MLT_indexes=np.where((MLT >= 0.5*imlt) & (MLT < 0.5*(imlt+1)))[0]
                    overlap=list(set(L_indexes) & set(MLT_indexes))
                    E_Spectra_sorted[imlt][iL]=E_Spectra_sorted[imlt][iL]+list(ESP[overlap][ESP[overlap]!=0])# reduce number of nans
                    for ifreq in range(index100-index10):
                        # make sure this is doing what it is supposed to 
                        Polarization_sorted[imlt][iL][ifreq]=Polarization_sorted[imlt][iL][ifreq]+list(polarizationT[ifreq][overlap][ESP[overlap]!=0])
      except(IndexError):
                print('ERROR!')
                print(dt1)
      dt1=dt1+datetime.timedelta(days=1)

      if dt1.month != month_cur:
          subdir_name='Sorted_Broadband_WFR'
          #print(curDate + 'freq='+str(FREQUENCIES[ifreq]))
          os.chdir('/Users/loisks/Desktop/ResearchProjects/PlasmaWaves/')
          pickling.hdf5_data_save(E_Spectra_sorted, curDate+'_polarization_PSD_broadband_'+ name[idir], subdir_name,Lbins, mltbins)
          pickling.hdf5_data_save(Polarization_sorted, curDate+'_PSD_polarization_broadband_'+ name[idir], subdir_name,Lbins, mltbins)
          #pickling.hdf5_data_save(E_Spectra_sorted, curDate+'_polarization_broadband_'+ name[idir], subdir_name,Lbins, mltbins)
          #pickling.hdf5_data_save(Polarization_sorted, curDate+'_polarization_broadband_'+ name[idir], subdir_name,Lbins, mltbins)
          os.chdir('..')
          month_cur=dt1.month
          E_Spectra_sorted=[[ [] for x in range(Lbins)] for x in range(mltbins)]
          Polarization_sorted=[[ [[] for z in range(index100-index10)] for x in range(Lbins)] for x in range(mltbins)]

          
