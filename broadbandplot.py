# broadbandplot.py
#
# plot the 10-100 Hz Equatorial Noise
#
# LKS, SANSA, March 2016
#
# imports
import numpy as np
from spacepy import pycdf
import glob
import os
import pickle
import datetime
import matplotlib.pyplot as plt
from matplotlib.dates import HourLocator, DayLocator, DateFormatter, MinuteLocator
import matplotlib.collections as collections
import matplotlib.dates as dt
from dateutil.relativedelta import relativedelta
from matplotlib.colors import LogNorm
os.chdir('/Users/loisks/Desktop/Functions')
import plots
import pickling
os.chdir('/Users/loisks/Documents/ResearchProjects/PlasmaWaves/')

#
# globals
Lbins=11
mltbins=48
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
# Functions
# gets the number of months between two dates for reading in the files
date1='20130201' # start date
date2='20150401' # end date
index10=4
index100=24
dt0=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')
def diff_month(d1, d2):
        return (d1.year - d2.year)*12 + d1.month - d2.month
    #

total_month=diff_month(dt2, dt0)
n1=6
dt1=dt0 
#
# loop through each month
# find the bad data month and remove it - I think august 2014
cfilee=[[ [] for x in range(Lbins)] for x in range(mltbins)]
meanfe=[[ [] for x in range(Lbins)] for x in range(mltbins)]
medPol=[[ [] for x in range(Lbins)] for x in range(mltbins)]
for imonth in range(total_month):
    cur_date=str(dt1.month)+'_'+str(dt1.year)
    dt1=dt1+relativedelta(months=1)
    os.chdir('Sorted_Broadband_WFR')
    filename1=glob.glob(cur_date+'_polarization_E_WFR_Spectra_broadband'+'_rbsp-a.h5')
    filename2=glob.glob(cur_date+'_polarization_E_WFR_Spectra_broadband'+'_rbsp-b.h5')
    try:
        de=pickling.hdf5_data_open(filename1[0],Lbins,mltbins)
    except:
        de = [[ [] for x in range(Lbins)] for x in range(mltbins)]        
    try:
        de2=pickling.hdf5_data_open(filename2[0],Lbins,mltbins)        
    except:
        de2 = [[ [] for x in range(Lbins)] for x in range(mltbins)]        
    os.chdir('..')
    for imlt in range(mltbins):
        for iL in range(Lbins):
    # combine A and B
            cfilee[imlt][iL]=list(cfilee[imlt][iL])+list(de[imlt][iL])+list(de2[imlt][iL])
    #
    # okay now exclude anamolies 
for imlt in range(mltbins):
    for iL in range(Lbins):
        try:
            a=np.array(cfilee[imlt][iL])*1000
            #a[a>=.2]=np.nan
            meanfe[imlt][iL]=np.nanmedian(a) # to get mV/m
            #allFiles[ifreq][imlt][iL]=np.nanmedian(a)
                          #  except(TypeError, AttributeError):
        except(IndexError):
            meanfe[imlt][iL]=np.nan
            #allFiles[ifreq][imlt][iL]=np.nan
            
# Now plotstop
plots.fancy_plot(meanfe, mltbins, Lbins, -3, -1, 'mV/m', 'Wave_Intensity', 'low', 'E_Wave_EMFISIS_wave_intensity_broadband', 'BroadbandPlots', 'log', 'viridis')
 # polarization 

