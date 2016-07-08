# emfisis_b_sort.py
#
# get the magnetic field magnitude for quiet time
# sort by L and MLAT
#
# LKS January 2015
#
# imports
import numpy as np
from spacepy import pycdf
import glob
import os
import datetime
import spacepy.pybats.kyoto as spk
import itertools as itert
import matplotlib.dates as dates
from scipy.interpolate import interp1d
import pickle 
import h5py
import spacepy.datamodel as dm
import pandas as pd
os.chdir('/Users/loisks/Documents/Functions')
import pickling as pickling
os.chdir('/Users/loisks/Desktop/')
#
# parameters and starting conditions
date1='20130101'
date2='20150401'
dt0=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')
nLbins=11
nMLATbins=31
nmltbins=48
MLTbinsArr=np.linspace(0.25, 24.25,nmltbins+1)
LBins=np.linspace(1.5, 4, nLbins)+0.125
MLATbins=np.linspace(-20, 10, nMLATbins)+0.5
dir=['A', 'B']
dir2=['EMPHEMERIS_A', 'EMPHEMERIS_B']

#MLAT_bins=[-20, -15, -10, -5, 0, 5, 10, 15, 20]
name=['rbsp-a', 'rbsp-b']
name2=['rbspa', 'rbspb']
# pickle function
#
#
# begin loop 
for idir in range(len(dir)):
    dt1=dt0
    Bsorted=[[[] for i in range(nLbins)] for j in range(nmltbins)]
    while dt1 < dt2:
        try:
                date=datetime.datetime.strftime(dt1, '%Y%m%d')
                os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EMFISIS_'+dir[idir])
                emf_file='*'+date+'*.cdf'
                gemf=glob.glob(emf_file)
                pyf=pycdf.CDF(gemf[0])
                bTime=pd.DatetimeIndex(pyf['Epoch'][...])
                bMag=pyf['Magnitude'][...]
                dB=pd.DataFrame({'Bmag':bMag}, index=bTime)
                
                #
                # now have to get the emphermeris data
                os.chdir('..')
                os.chdir('EMPHEMERIS_'+dir[idir])
                emp_file='*'+date+'*.txt'
                gemp=glob.glob(emp_file)
                pyfemp=dm.readJSONheadedASCII(gemp[0])                
                L_emphem=pyfemp['L'] # L from emphemeris
                L_emphem=np.nanmedian(L_emphem, axis=1)
                L=np.array(L_emphem)
                Kp_emphem=pyfemp['Kp'] # Kp index from ephemeris file
                MLAT_emphem=pyfemp['CDMAG_MLAT'] # MLAT from ephemeris
                MLAT_emphem[MLAT_emphem<-100]=np.nan
                MLAT=np.array(MLAT_emphem)
                MLT_emphem=pyfemp['CDMAG_MLT'] # MLAT from ephemeris
                MLT_emphem[MLT_emphem<-100]=np.nan
                MLT=np.array(MLT_emphem)
                # emphemeris data is every minute
                rt = pd.period_range(date,periods=1441, freq='T').to_timestamp()
                #
                # resample the B data
                B=np.array(dB['Bmag'].resample('1min', how='mean').reindex(index=rt,fill_value=np.nan))

                #
                # now sort appropriately
                for imlt in range(nmltbins):
                    for iL in range(nLbins):
                        L_indexes=np.where((L >= 1.25+.125+0.25*iL) & (L < 1.25+0.125+0.25*(iL+1)))[0]
                        MLT_indexes=np.where((MLT >= 0.5*imlt) & (MLT < 0.5*(imlt+1)))[0]
                        matches=list(set(L_indexes) & set(MLT_indexes))
                        Bsorted[imlt][iL]+=list(B[matches])
                
        except(IndexError):
                print(dt1)
        dt1=dt1+datetime.timedelta(days=1)
    print('saving first file')
    os.chdir('/Users/loisks/Desktop/ResearchProjects/PlasmaWaves')
    subdir_name='Sorted_BMag'
    if not os.path.exists(subdir_name):
        os.umask(0) # unmask if necessary
        os.makedirs(subdir_name) 
    pickling.hdf5_data_save(Bsorted, 'BSorted_MLT_'+ dir[idir], subdir_name,  nLbins, nmltbins)
    os.chdir('..')


