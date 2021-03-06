# broadband_10Hz_100Hz_sort.py
#
# looking at total PSD between 10 and 100 Hz measurements in EMFISIS
# sum the amplitudes of the linearly polarized portion
#
# LKS, SANSA 2016 March
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
indexH2=24
indexH1=0
#for ifreq in range(len(FREQUENCIES)):
# frequencies from 70 to 1000 Hz
for idir in range(2):# CHANGE FOR BOTH A AND B
   # start overarching loop through days
    dt1=dt0
    EH_Spectra_sorted=[[ [] for x in range(Lbins)] for x in range(mltbins)]
    PolarizationH_sorted=[[[ [] for y in range(indexH2-indexH1)] for x in range(Lbins)] for x in range(mltbins)]
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
        polarizationH = [ [] for i in range(indexH2-indexH1)]
        MLTH = [ [] for i in range(indexH2-indexH1)]
        LH = [ [] for i in range(indexH2-indexH1)]
        weightsH = [ [] for i in range(indexH2-indexH1)]
        E_SpectraH = [[] for i in range(indexH2-indexH1) ]

        # H loop
        for ii in range(24):
            #ii = ifreq1+indexH1
            polarizationH[ii]=gM.mag_SVD(pyf, ii)
            Eu_SpectraH=np.array(np.swapaxes(pyf['EuEu'][...],1 , 0)[ii]) # spectra, 65 frequncies 
            Ev_SpectraH=np.array(np.swapaxes(pyf['EvEv'][...],1 , 0)[ii])
            Ew_SpectraH=np.array(np.swapaxes(pyf['EwEw'][...],1 , 0)[ii])
            #
            #
            Eu_SpectraH[Eu_SpectraH<0]=np.nan
            Ev_SpectraH[Ev_SpectraH<0]=np.nan
            Ew_SpectraH[Ew_SpectraH<0]=np.nan
            
            #
            # change nans to zero
            Eu_SpectraH[np.isnan(Eu_SpectraH)]=0
            Ev_SpectraH[np.isnan(Eu_SpectraH)]=0
            Ew_SpectraH[np.isnan(Eu_SpectraH)]=0
            #
            # change places where polarization > .2 to 0
            #polarizationH[ii][polarizationH[ii] > .2]=0
            E_SpectraH[ii]=0.5*np.sqrt(np.array((Eu_SpectraH+Ev_SpectraH+Ew_SpectraH)))
        EHintegrated=np.zeros(len(E_SpectraH[0]))
        for ifreq1 in range(23):
            #ii = ifreq1+indexH1
            #
            # Riemann Sum
            EHintegrated+=np.array( ((E_SpectraH[ifreq1]+E_SpectraH[ifreq1+1])/2.0))*(FREQUENCIES[ifreq1+1]-FREQUENCIES[ifreq1])

        #
        # now get polarization
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
        #
        L = np.array(L_emphem)
        MLT = np.array(MLT_emphem)
        #
        #
        Len= np.array(list(filter(lambda EHintegrated: EHintegrated != 0, L)))
        MLTen= np.array(list(filter(lambda EHintegrated: EHintegrated != 0, MLT)))
        for ij in range(24):
            temp=np.array(polarizationH[ij])
            MLTH[ij]=np.array(list(filter(lambda temp: temp <= 0.2, MLT)))
            LH[ij]=np.array(list(filter(lambda temp: temp <= 0.2, L)))
            polarizationH[ij]=np.array(list(filter(lambda temp : temp <= 0.2, temp)))           
        EHintegrated = np.array(list(filter(lambda EHintegrated: EHintegrated != 0, EHintegrated)))
        
        # now in the loop 
        os.chdir('/Users/loisks/Desktop/ResearchProjects/PlasmaWaves/')
        #
        # screen for ellipticity and polarization
        # linearly polarized waves = thsvd -> 90 and ellsvd -> 0
# now go through the fluxes
        for imlt in range(mltbins):
            for iL in range(Lbins):
                    L_indexes=np.where((Len >= 1.5+0.25*iL) & (Len < 1.5+0.25*(iL+1)))[0]
                    MLT_indexes=np.where((MLTen >= 0.5*imlt) & (MLTen < 0.5*(imlt+1)))[0]
                    overlapEn=list(set(L_indexes) & set(MLT_indexes))
                    EH_Spectra_sorted[imlt][iL]=EH_Spectra_sorted[imlt][iL]+list(EHintegrated[overlapEn])
                    for ifreq in range(24):
                        LinPol=np.where((LH[ij] >= 1.5+0.25*iL) & (LH[ij] < 1.5+0.25*(iL+1)))[0]
                        MLTPol=np.where((MLTH[ij] >= 0.5*imlt) & (MLTH[ij] < 0.5*(imlt+1)))[0]
                        overlapP = list(set(LinPol) & set(MLTPol))
                        temp=np.swapaxes(PolarizationH_sorted,0,2)[ifreq]
                        temp=np.swapaxes(temp, 0,1)                        
                        PolarizationH_sorted[imlt][iL][ifreq]=PolarizationH_sorted[imlt][iL][ifreq]+list(temp[overlapP])
                   
      except(IndexError):
      #except(OSError):
          
                print('ERROR!')
                print(dt1)
      dt1=dt1+datetime.timedelta(days=1)
      if dt1.month != month_cur:
          subdir_name='Sorted_Broadband_WFR_Species'
          os.chdir('/Users/loisks/Desktop/ResearchProjects/PlasmaWaves/')
          pickling.hdf5_data_save(EH_Spectra_sorted, curDate+'_integrated_2_100_broadband_'+ name[idir], subdir_name,Lbins, mltbins)

          for ifreq in range(24):
              temp=np.swapaxes(PolarizationH_sorted,0,2)[ifreq]
              temp=np.swapaxes(temp, 0,1)
              
              pickling.hdf5_data_save(temp, curDate+'_2_100_polarization_broadband_'+ name[idir]+'_freq='+str(FREQUENCIES[ifreq]), subdir_name,Lbins, mltbins)
          
          #pickling.hdf5_data_save(E_Spectra_sorted, curDate+'_polarization_broadband_'+ name[idir], subdir_name,Lbins, mltbins)
          #pickling.hdf5_data_save(Polarization_sorted, curDate+'_polarization_broadband_'+ name[idir], subdir_name,Lbins, mltbins)
          os.chdir('..')
          month_cur=dt1.month
          EH_Spectra_sorted=[[ [] for x in range(Lbins)] for x in range(mltbins)]
          PolarizationH_sorted=[[[ [] for y in range(24)  ] for x in range(Lbins)] for x in range(mltbins)]
          
