# polarizationCompare.py
#
# compare the right hand polarization for different MLTs at one L-Shell
# position via line plot
# so frequencies vs polarization
# and frequencies vs wave amplitude for pol > 0.7
# sounds good to me
#
# LKS, March 2016, Tromso, Norway. I have such a headache, so we'll see
# how this code goes
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

def depickle(name):
    with open(name, 'rb') as f:
        u = pickle._Unpickler(f)
        u.encoding = 'latin1'
        outf = u.load()
        return outf
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
date1='20130201' # start date
date2='20150401' # end date
dt0=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')
def diff_month(d1, d2):
        return (d1.year - d2.year)*12 + d1.month - d2.month
    #

total_month=diff_month(dt2, dt0)
n1=6
allFiles=[[[ np.nan for x in range(Lbins)] for x in range(mltbins)] for j in range(len(FREQUENCIES))]
PolAll=np.ones((4,50))*np.nan
IntensityAll=np.ones((4,50))*np.nan
IntensityLow=np.ones((4,50))*np.nan
for ifreq in range(0,50):
 dt1=dt0 
    #
    # loop through each month
    # find the bad data month and remove it - I think august 2014
 cfilee=[[ [] for x in range(Lbins)] for x in range(mltbins)]
 cPol=[[ [] for x in range(Lbins)] for x in range(mltbins)]
 meanfe=[[ [] for x in range(Lbins)] for x in range(mltbins)]
 medPol=[[ [] for x in range(Lbins)] for x in range(mltbins)]
 perPol=[[ 0 for x in range(Lbins)] for x in range(mltbins)]
 for imonth in range(total_month):
    cur_date=str(dt1.month)+'_'+str(dt1.year)
    dt1=dt1+relativedelta(months=1)
    os.chdir('Sorted_WFR_Polarization')
    filename1=glob.glob(cur_date+'_polarization_E_WFR_Spectra_freq='+str(FREQUENCIES[ifreq])+'_rbsp-a.h5')
    filename2=glob.glob(cur_date+'_polarization_E_WFR_Spectra_freq='+str(FREQUENCIES[ifreq])+'_rbsp-b.h5')
    filePol1=glob.glob(cur_date+'_polarization_freq='+str(FREQUENCIES[ifreq])+'_rbsp-a.h5')
    filePol2=glob.glob(cur_date+'_polarization_freq='+str(FREQUENCIES[ifreq])+'_rbsp-b.h5')
    try:
        de=pickling.hdf5_data_open(filename1[0],Lbins,mltbins)
        deP1=pickling.hdf5_data_open(filePol1[0],Lbins,mltbins)
    except:
        de = [[ [] for x in range(Lbins)] for x in range(mltbins)]
        deP1 = [[ [] for x in range(Lbins)] for x in range(mltbins)]
    try:
        de2=pickling.hdf5_data_open(filename2[0],Lbins,mltbins)
        deP2=pickling.hdf5_data_open(filePol2[0],Lbins,mltbins)
    except:
        de2 = [[ [] for x in range(Lbins)] for x in range(mltbins)]
        deP2 = [[ [] for x in range(Lbins)] for x in range(mltbins)]
    os.chdir('..')
    for imlt in range(mltbins):
        for iL in range(Lbins):
    # combine A and B
                cfilee[imlt][iL]=list(cfilee[imlt][iL])+list(de[imlt][iL])+list(de2[imlt][iL])
                cPol[imlt][iL]=list(cPol[imlt][iL])+list(deP1[imlt][iL])+list(deP2[imlt][iL])
    
    PolAll[0,ifreq]=np.nanmedian(cPol[5][3])
    PolAll[1,ifreq]=np.nanmedian(cPol[17][3])
    PolAll[2,ifreq]=np.nanmedian(cPol[29][3])
    PolAll[3,ifreq]=np.nanmedian(cPol[41][3])
    #
    IntensityAll[0,ifreq]=np.nanmedian(np.array(cfilee[5][3])[np.array(cPol[5][3])>0.7])
    IntensityAll[1,ifreq]=np.nanmedian(np.array(cfilee[17][3])[np.array(cPol[17][3])>0.7])
    IntensityAll[2,ifreq]=np.nanmedian(np.array(cfilee[29][3])[np.array(cPol[29][3])>0.7])
    IntensityAll[3,ifreq]=np.nanmedian(np.array(cfilee[41][3])[np.array(cPol[41][3])>0.7])
    #
    IntensityLow[0,ifreq]=np.nanmedian(np.array(cfilee[5][3])[np.array(cPol[5][3])<0.2])
    IntensityLow[1,ifreq]=np.nanmedian(np.array(cfilee[17][3])[np.array(cPol[17][3])<0.2])
    IntensityLow[2,ifreq]=np.nanmedian(np.array(cfilee[29][3])[np.array(cPol[29][3])<0.2])
    IntensityLow[3,ifreq]=np.nanmedian(np.array(cfilee[41][3])[np.array(cPol[41][3])<0.2])
    # 
#
# now go into line plot
xfreqs=FREQUENCIES[0:50]
f, ax=plt.subplots(3,sharex=True)
plt.subplots_adjust(right=0.7, top=0.9, bottom=0.15)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
ax[0].set_title('Median Polarization', fontsize=18, fontweight='bold')
ax[0].plot(xfreqs,PolAll[0], lw=2, c='DarkOrange')
ax[0].plot(xfreqs,PolAll[1], lw=2, c='blue')
ax[0].plot(xfreqs,PolAll[2], lw=2, c='turquoise')
ax[0].plot(xfreqs,PolAll[3], lw=2, c='DarkPink')
ax[0].set_xscale('log')
ax[0].set_xlim(50,1000)
ax[0].axvline(x=150, c='k', lw=2, ls='--')
ax[0].axvline(x=600, c='k', lw=2, ls='--')
ax[0].set_ylabel('Polarization')

ax[1].set_title('Power Spectral Density of Right Hand Waves', fontsize=18, fontweight='bold')
ax[1].set_xlim(50,1000)
ax[1].plot(xfreqs,IntensityAll[0], lw=2, c='DarkOrange')
ax[1].plot(xfreqs,IntensityAll[1], lw=2, c='b')
ax[1].plot(xfreqs,IntensityAll[2], lw=2, c='turquoise')
ax[1].plot(xfreqs,IntensityAll[3], lw=2, c='DarkPink')
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].axvline(x=150, c='k', lw=2, ls='--')
ax[1].axvline(x=600, c='k', lw=2, ls='--')
ax[1].set_ylim(1e-13, 1e-10)
ax[1].set_ylabel('V$^{2}$/m$^{2}$/Hz')
#ax[1].set_xlabel('Frequency [Hz]')
#
ax[2].set_title('Power Spectral Density of Linearly Polarized Waves', fontsize=18, fontweight='bold')
ax[2].set_xlim(50,1000)
ax[2].plot(xfreqs,IntensityLow[0], lw=2, c='DarkOrange')
ax[2].plot(xfreqs,IntensityLow[1], lw=2, c='b')
ax[2].plot(xfreqs,IntensityLow[2], lw=2, c='turquoise')
ax[2].plot(xfreqs,IntensityLow[3], lw=2, c='DarkPink')
ax[2].set_xscale('log')
ax[2].set_yscale('log')
ax[2].axvline(x=150, c='k', lw=2, ls='--')
ax[2].axvline(x=600, c='k', lw=2, ls='--')
ax[2].set_ylim(1e-13, 1e-9)
ax[2].set_ylabel('V$^{2}$/m$^{2}$/Hz')
ax[2].set_xlabel('Frequency [Hz]')

plt.legend(['MLT = 3', 'MLT = 9', 'MLT = 15', 'MLT = 21'], bbox_to_anchor=[1.5, 2.3])
subdir_name='Comparison_Pol'
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name) 
os.chdir(subdir_name)#
f.set_size_inches(13,9)
plt.savefig('polarization_compare.pdf')
plt.close(f)
os.chdir('..')
