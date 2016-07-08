# broadbandintegrateplot.py
#
# plot the two species integrated frequency
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

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
     r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
     The Savitzky-Golay filter removes high frequency noise from data.
     It has the advantage of preserving the original shape and
     features of the signal better than other types of filtering
     approaches, such as moving averages techniques.
     Parameters
     ----------
     y : array_like, shape (N,)
         the values of the time history of the signal.
     window_size : int
         the length of the window. Must be an odd integer number.
     order : int
         the order of the polynomial used in the filtering.
         Must be less then `window_size` - 1.
     deriv: int
         the order of the derivative to compute (default = 0 means only smoothing)
     Returns
     -------
     ys : ndarray, shape (N)
         the smoothed signal (or it's n-th derivative).
     Notes
     -----
     The Savitzky-Golay is a type of low-pass filter, particularly
     suited for smoothing noisy data. The main idea behind this
     approach is to make for each point a least-square fit with a
     polynomial of high order over a odd-sized window centered at
     the point.
     Examples
     --------
     t = np.linspace(-4, 4, 500)
     y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
     ysg = savitzky_golay(y, window_size=31, order=4)
     import matplotlib.pyplot as plt
     plt.plot(t, y, label='Noisy signal')
     plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
     plt.plot(t, ysg, 'r', label='Filtered signal')
     plt.legend()
     plt.show()
     References
     ----------
     .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
        Data by Simplified Least Squares Procedures. Analytical
        Chemistry, 1964, 36 (8), pp 1627-1639.
     .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
        W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
        Cambridge University Press ISBN-13: 9780521880688
     """
     import numpy as np
     from math import factorial
     
     try:
         window_size = np.abs(np.int(window_size))
         order = np.abs(np.int(order))
     except ValueError, msg:
         raise ValueError("window_size and order have to be of type int")
     if window_size % 2 != 1 or window_size < 1:
         raise TypeError("window_size size must be a positive odd number")
     if window_size < order + 2:
         raise TypeError("window_size is too small for the polynomials order")
     order_range = range(order+1)
     half_window = (window_size -1) // 2
     # precompute coefficients
     b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
     m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
     # pad the signal at the extremes with
     # values taken from the signal itself
     firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
     lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
     y = np.concatenate((firstvals, y, lastvals))
     return np.convolve( m[::-1], y, mode='valid')
    #
# Functions
# gets the number of months between two dates for reading in the files
date1='20130201' # start date
date2='20150401' # end date

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
cfileeH=[[ [] for x in range(Lbins)] for x in range(mltbins)]
meanfeH=[[ [] for x in range(Lbins)] for x in range(mltbins)]
medPolH=[[ [] for x in range(Lbins)] for x in range(mltbins)]
cfileeHe=[[ [] for x in range(Lbins)] for x in range(mltbins)]
meanfeHe=[[ [] for x in range(Lbins)] for x in range(mltbins)]
medPolHe=[[ [] for x in range(Lbins)] for x in range(mltbins)]
for imonth in range(total_month):
    cur_date=str(dt1.month)+'_'+str(dt1.year)
    dt1=dt1+relativedelta(months=1)
    os.chdir('Sorted_Broadband_WFR_Species')
    filename1H=glob.glob(cur_date+'_H_polarization_integrated_broadband'+'_rbsp-a.h5')
    filename2H=glob.glob(cur_date+'_H_polarization_integrated_broadband'+'_rbsp-b.h5')
    filename1He=glob.glob(cur_date+'_He_polarization_integrated_broadband'+'_rbsp-a.h5')
    filename2He=glob.glob(cur_date+'_He_polarization_integrated_broadband'+'_rbsp-b.h5')
    try:
        deH=pickling.hdf5_data_open(filename1H[0],Lbins,mltbins)
    except:
        deH = [[ [] for x in range(Lbins)] for x in range(mltbins)]        
    try:
        de2H=pickling.hdf5_data_open(filename2H[0],Lbins,mltbins)        
    except:
        de2H = [[ [] for x in range(Lbins)] for x in range(mltbins)]
    try:
        deHe=pickling.hdf5_data_open(filename1He[0],Lbins,mltbins)
    except:
        deHe = [[ [] for x in range(Lbins)] for x in range(mltbins)]        
    try:
        de2He=pickling.hdf5_data_open(filename2He[0],Lbins,mltbins)        
    except:
        de2He = [[ [] for x in range(Lbins)] for x in range(mltbins)]    
    os.chdir('..')
    for imlt in range(mltbins):
        for iL in range(Lbins):
    # combine A and B
            cfileeH[imlt][iL]=list(cfileeH[imlt][iL])+list(deH[imlt][iL])+list(de2H[imlt][iL])
            cfileeHe[imlt][iL]=list(cfileeHe[imlt][iL])+list(deHe[imlt][iL])+list(de2He[imlt][iL])
    #
    # okay now exclude anamolies 
for imlt in range(mltbins):
    for iL in range(Lbins):
        try:
            aH=np.array(cfileeH[imlt][iL])*1000
            aHe=np.array(cfileeHe[imlt][iL])*1000
            #a[a>=.2]=np.nan
            meanfeH[imlt][iL]=np.nanmedian(aH) # to get mV/m
            meanfeHe[imlt][iL]=np.nanmedian(aHe) # to get mV/m
            #allFiles[ifreq][imlt][iL]=np.nanmedian(a)
                          #  except(TypeError, AttributeError):
        except(IndexError):
            meanfeH[imlt][iL]=np.nan
            meanfeHe[imlt][iL]=np.nan
            #allFiles[ifreq][imlt][iL]=np.nan
            
# Now plotstop
plots.fancy_plot(meanfeH, mltbins, Lbins, -3, -1, 'mV/m', 'Wave_Intensity', 'low', 'H_EMFISIS_wave_intensity_broadband', 'BroadbandPlots', 'log', 'viridis')
plots.fancy_plot(meanfeHe, mltbins, Lbins, -3, -1, 'mV/m', 'Wave_Intensity', 'low', 'He_EMFISIS_wave_intensity_broadband', 'BroadbandPlots', 'log', 'viridis')
 # polarization
 #
 #make a quick line plot at L = 2
 #
spe=['P', 'He']
totalDen=[[[ [] for j in range(len(spe))] for x in range(Lbins)] for x in range(mltbins)]
medDen=[[[ 0 for j in range(len(spe))] for x in range(Lbins)] for x in range(mltbins)]
os.chdir('PlasmaDensity10eV')
for iSpe in range(len(spe)):
    files=glob.glob('*'+'_species='+spe[iSpe])
    for iFile in range(len(files)):
        dataA=pickle.load(open(files[0], "rb"))
        dataB=pickle.load(open(files[0], "rb"))
        for iMLT in range(mltbins):
            for iL in range(Lbins):
                try:
                    totalDen[iMLT][iL][iSpe]+=list(np.array(dataA[iMLT][iL]))+list(np.array(dataB[iMLT][iL]))
                except(OSError):
                    totalDen[iMLT][iL][iSpe]=np.nan
    # now get the median
    for iMLT in range(mltbins):
        for iL in range(Lbins):
            medDen[iMLT][iL][iSpe]=np.nanmedian(totalDen[iMLT][iL][iSpe])
os.chdir('..')
for iL in range(2, 7):
    fig=plt.figure()
    ax=fig.add_subplot(111)
    plt.subplots_adjust(right=0.6, top=0.65, bottom=0.15)
    ax.set_xlabel('MLT', fontsize=25, fontweight='bold')
    ax.set_ylabel('mV/m', fontsize=25, fontweight='bold')
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 22}
    plt.rc('font', **font)
    y1=np.swapaxes(meanfeH,1,0)[iL]
    y2=np.swapaxes(meanfeHe,1,0)[iL]
    yHden=np.swapaxes(medDen,2,0)[0][iL]
    yHeden=np.swapaxes(medDen,2,0)[1][iL]
    #y1=savitzky_golay(np.swapaxes(meanfeH,1,0)[iL], 5, 1)
    #y2=savitzky_golay(np.swapaxes(meanfeHe,1,0)[iL], 5, 1)
    #yHden=savitzky_golay(np.swapaxes(medDen,2,0)[0][iL], 5, 1)
    #yHeden=savitzky_golay(np.swapaxes(medDen,2,0)[1][iL], 5, 1)
    #
    # shift the data by x hours
    nmltbins=mltbins/2
    ty1=np.zeros(nmltbins)
    ty2=np.zeros(nmltbins)
    #
    # first average by hour bins
    for iN in range(nmltbins):
      ty1[iN]=(y1[iN]*2 + (y1[iN]*2)+1)/2.
      ty2[iN]=(y2[iN]*2 + (y2[iN]*2)+1 )/2.
    newy1=np.zeros(nmltbins)
    newy2=np.zeros(nmltbins)
    offset=5 # number of hours
    for imlt in range(nmltbins-offset):
        newy1[imlt]=ty1[imlt+offset]#[imlt+3]
        newy2[imlt]=ty2[imlt+offset]
    for imlt in range(offset):
        newy1[nmltbins-offset + imlt]=ty1[imlt]
        newy2[nmltbins-offset + imlt]=ty2[imlt]

     
    p1,=ax.plot(np.linspace(0,24,nmltbins),newy1, c='b', lw=2, ls='--')
    p2,=ax.plot(np.linspace(0,24,nmltbins),newy2, c='r', lw=2, ls='--')
    print newy1
    print newy2
    
    ax.set_xlim(0,24)
    ax1=ax.twinx()
    ax1.set_yscale('log')
    ax1.set_xlim(0,24)
    ax1.set_ylabel('1-10 eV Density [cm$^{-3}$]')
    p3,=ax1.plot(np.linspace(0,24,mltbins-3),yHden[:-3], c='b', lw=2)
    p4,=ax1.plot(np.linspace(0,24,mltbins-3),yHeden[:-3], c='r', lw=2)
   # ax.set_ylim(0.02, .09)
    ax.set_xlim(0,24)
    plt.legend([p1,p2,p3,p4],['H$^{+}$ PSD', 'He$^{+}$ PSD', 'H$^{+}$ Density',  'He$^{+}$ Density'], bbox_to_anchor=[1.75, 0.75], fontsize=20)
    subdir_name='BroadBand_LinePlots'
    if not os.path.exists(subdir_name):
        os.umask(0) # unmask if necessary
        os.makedirs(subdir_name) 
    os.chdir(subdir_name)#
    fig.set_size_inches(13,9)
    plt.savefig('L='+str(iL*0.25 + 1.5)+'_offset='+str(offset/2.0)+'.pdf')
    plt.close(fig)
    os.chdir('..')
