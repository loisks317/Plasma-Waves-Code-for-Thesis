# WFR_spectra.py
#
# low frequency wave spectra from EMFISIS for case studies
#
# LKS MAY 2015 OMG FIVE DAYS TILL QUALS.
#
import numpy as np 
import os
os.chdir('/Users/loisks/Documents/Functions')
import pickling as pickling
os.chdir('/Users/loisks/Desktop/PA_Distributions/stat_corr/')
import glob
import datetime
from dateutil.relativedelta import relativedelta 
import scipy
import pickle
import matplotlib as mpl
from spacepy import pycdf
import scipy.stats
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from matplotlib.dates import HourLocator, DayLocator, DateFormatter, MinuteLocator
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter 
import matplotlib.collections as collections
import matplotlib.dates as dt
days = DayLocator(interval=1) 
#hours = HourLocator(interval=1) 
#hours2 = HourLocator(interval=3)
hours=MinuteLocator(interval=30)
hours2=MinuteLocator(interval=15)
daysFmt = DateFormatter('%H:%M')
#
# dates and stuff
date='20130702'
hour1='18'
hour2='20'
FREQUENCIES=[  2.13599992e+00,   4.27199984e+00,   6.40899992e+00,
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
ENERGIES=[0.99, 1.19, 1.33, 1.55, 1.82, 2.18, 2.53, 2.95, 3.38, 3.94, 4.64, 5.34, 6.26, 7.31, 8.51, 9.91]
echan=12 # 8 = 2.5 eV
#fchan=17
#frequency=str(FREQUENCIES[fchan]) + ' Hz'
energy=str(ENERGIES[echan]) + ' eV'

os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EMFISIS_WFR_A')
emf_file=glob.glob('rbsp-a_WFR-spectral-matrix-diagonal_emfisis-L2_'+date+'*.cdf')
emf_pyf=pycdf.CDF(emf_file[0])
emf_epoch=emf_pyf['Epoch'][...]
BuBu=emf_pyf['BuBu'][...]
BvBv=emf_pyf['BvBv'][...]
BwBw=emf_pyf['BwBw'][...]
EuEu=emf_pyf['EuEu'][...]
EvEv=emf_pyf['EvEv'][...]
EwEw=emf_pyf['EwEw'][...]
Emag=[[ 0 for x in xrange(len(FREQUENCIES))] for x in xrange(len(EuEu))]
#BuBu=np.swapaxes(BuBu, 1, 0) 
#BvBv=np.swapaxes(BvBv, 1, 0)
#BwBw=np.swapaxes(BwBw, 1, 0)
#EuEu=np.swapaxes(EuEu, 1, 0)
#EvEv=np.swapaxes(EvEv, 1, 0)
#EwEw=np.swapaxes(EwEw, 1, 0)
for itime in range(len(emf_epoch)):
    for ifreq in range(len(FREQUENCIES)):
        #Bmag=np.sqrt(np.array(BuBu**2) + np.array(BvBv**2) + np.array(BwBw**2))
        Emag[itime][ifreq]=np.sqrt(EuEu[itime][ifreq]**2 + EvEv[itime][ifreq]**2 + EwEw[itime][ifreq]**2)
#
os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EMFISIS_A')                
emf_file= 'rbsp-a_magnetometer_4sec-gse_emfisis-L3_'+date+'*'
emf_pyf2=pycdf.CDF(glob.glob(emf_file)[0])
emf_B=emf_pyf2['Magnitude'][...]
gyro=(1.6*1e-19)*(np.array(emf_B)*1e-9)/(1.67*1e-27)
magneto_time=emf_pyf2['Epoch'][...]
os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_A')
hope_file=glob.glob('rbspa_rel02_ect-hope-PA-L3_'+date+'*.cdf')
hope_pyf=pycdf.CDF(hope_file[0])
hope_epoch=hope_pyf['Epoch_Ion'][...]
hope_Hflux=hope_pyf['FPDO'][...]
L=hope_pyf['L_Ion'][...]
MLT=hope_pyf['MLT_Ion'][...]
os.chdir('/Users/loisks/Desktop/PA_Distributions/wave_spectra/')
flux={}
n1=2
n2=10
#
# choose specifice frequencies and energies to plot
for iFlux in range(n1, n2):
    flux[iFlux]=np.swapaxes(hope_Hflux, 1, 0)[iFlux] # E=2.9 eV
#
# modify for hours instead of a full day of data
datetime1=datetime.datetime(int(date[0:4]), int(date[4:6]), int(date[6:8]), int(hour1), 0, 0)
datetime2=datetime.datetime(int(date[0:4]), int(date[4:6]), int(date[6:8]), int(hour2), 0, 0)
try:
        start_hope=np.where(hope_epoch>datetime1)[0][0]
        start_emf=np.where(emf_epoch>datetime1)[0][0]
except(IndexError):
        start_hope=0; start_emf=0
try:
        end_hope=np.where(hope_epoch>datetime2)[0][0]
        end_emf=np.where(emf_epoch>datetime2)[0][0]
except(IndexError):
        end_hope=len(hope_epoch)-1; end_emf=len(emf_epoch)-1

vmin=-13
vmax=-10
vmin1=1e-13
vmax1=1e-10

# let's plot this
fig=plt.figure()
fig.gca().xaxis.set_major_locator(hours)
fig.gca().xaxis.set_major_formatter(daysFmt)
fig.gca().xaxis.set_minor_locator(hours2)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
#ax1=fig.add_subplot(111)
#
# define axes
left, width = 0.1, 0.73
bottom, height = 0.3, 0.45
left_h= left+width+0.03
bottom_h= bottom+height+0.03

rect_scatter= [left, bottom, width, height]
rect_top = [left, bottom_h, width, 0.18]
fig=plt.figure(1, figsize=(8,8))
ax1 = plt.axes(rect_scatter)
axtop = plt.axes(rect_top)
ax1.set_yscale('log')
ax1.set_ylim(50, 1000)
#ax1.set_xlim(0,24)
ax1.set_xlabel('Time')
ax1.set_ylabel('Frequency (Hz)', fontweight='bold')
plt.subplots_adjust(left=0.11, right=0.7, top=0.92, bottom=0.28)
timeArr=emf_epoch[start_emf:end_emf]
Emag_t=Emag[start_emf:end_emf]
allFreq=FREQUENCIES
pshe=ax1.pcolormesh(timeArr,np.array(allFreq),np.swapaxes(np.array(Emag_t), 1, 0), cmap='jet',
                    norm=LogNorm(vmin=vmin1, vmax=vmax1)
                        )
ax1.plot(magneto_time,6*np.array(gyro)/(2*np.pi), lw=3, color='silver') 
#ax1.plot(MLT_bins, np.array(np.ones(len(MLT_bins)))*1e4, ls='--', lw=3, color='k')
cbaxes = fig.add_axes([0.84, 0.3, 0.03, 0.45])  
#cb = plt.colorbar(pshe, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(vmin, vmax+1)) 
cb = plt.colorbar(pshe, cax = cbaxes)
cb.set_label('V$^{2}$/m$^{2}$/Hz', fontsize=25, fontweight='bold')
ax1.tick_params(axis='both', which='major', labelsize=22)
cb.ax.tick_params(labelsize=25)

ax2 = ax1.twinx()
ax2.set_yscale('log')
ax2.set_ylim(1e5, 1e8)
#ax2.set_ylabel('Differential # Flux', fontsize=22, fontweight='bold')
p2=ax2.scatter(hope_epoch[start_hope:end_hope][flux[n1][start_hope:end_hope]>0], flux[n1][start_hope:end_hope][flux[n1][start_hope:end_hope]>0], s=20, color='k', label= 'Flux at 3 eV', alpha=0)

#plt.setp(ax2.get_yticklabels(), visible=False)
#plt.setp(ax2.get_xticklabels(), visible=False)
#p2=ax2.scatter(hope_epoch[start_hope:end_hope][flux[start_hope:end_hope]>0], flux[start_hope:end_hope][flux[start_hope:end_hope]>0], s=10, color='greenyellow', label= 'Flux at 3 eV')
#ax2.set_ylabel
#ax2.yaxis.label.set_color('blue')
a= L[start_hope:end_hope] <= 3
b= L[start_hope:end_hope] >= 2
c= MLT[start_hope:end_hope] <= 4
d= MLT[start_hope:end_hope] >= 1
overlap=np.zeros(len(L[start_hope:end_hope]))
for io in range(len(overlap)):
    if ((a[io]==True) & (b[io]==True) & (c[io]==True) & (d[io]==True)):
        overlap[io]=True
    else:
        overlap[io]=False
#
collection_L = collections.BrokenBarHCollection.span_where(
           dt.date2num(hope_epoch[start_hope:end_hope]), ymin=0.0, ymax=1e12, where=overlap, facecolor='black', alpha=0.25)

ax2.add_collection(collection_L)
ax2.set_xlim(hope_epoch[start_hope], hope_epoch[end_hope] )
plt.setp(ax2.get_yticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)

cmap = [ 'lightseagreen', 'lime', 'blue', 'magenta', 'darkorange']
nsteps=(n2-n1)/2
labels=[]
for iPlot in range(n1/2, n2/2):
    axtop.scatter(hope_epoch[start_hope:end_hope][flux[iPlot*2][start_hope:end_hope]>0], flux[iPlot*2][start_hope:end_hope][flux[iPlot*2][start_hope:end_hope]>0], s=10, color=cmap[iPlot])
    labels.append(str(ENERGIES[iPlot*2]) + ' eV')
    
axtop.set_xlim(ax1.get_xlim() )
axtop.set_ylim(1e5,1e10)
axtop.set_yscale('log')
axtop.set_ylabel('Diff # Flux', fontsize='20', fontweight='bold')
axtop.legend(labels,bbox_to_anchor=(1.01, 1.05), loc=2, fontsize='15')
collection_L = collections.BrokenBarHCollection.span_where(
           dt.date2num(hope_epoch[start_hope:end_hope]), ymin=0.0, ymax=1e12, where=overlap, facecolor='black', alpha=0.25)

axtop.add_collection(collection_L)
axtop.tick_params(axis='x', bottom='off', top='off', labelbottom='off')

fig.gca().xaxis.set_major_locator(hours)
fig.gca().xaxis.set_major_formatter(daysFmt)
fig.gca().xaxis.set_minor_locator(hours2)


#
# Add L LABEL
ax3 = fig.add_axes(ax2.get_position())
ax3.patch.set_visible(False)
ax3.yaxis.set_visible(False)
ax3.set_yscale('log')
#ax3.set_ylim(1e3, 1e9)
for spinename in ax3.spines:
    if spinename != 'bottom':
        ax3.spines[str(spinename)].set_visible(False)
ax3.spines['bottom'].set_position(('outward', 65))
L_conc=[round(x, 2) for x in L[start_hope:end_hope]]
decoy=np.swapaxes(flux[n1], 1, 0)
p3=ax3.plot(hope_epoch[start_hope:end_hope],  decoy[start_hope:end_hope], 'g', alpha=0)
# I hate doing this but
spacing = len(L_conc)/7
ax3.set_xticks(hope_epoch[start_hope:end_hope][::spacing])
ax3.set_xticklabels(L_conc[::spacing])

ax3.set_xlabel('L-Shell', fontsize=22, fontweight='bold')
#
# Add MLT Label
ax5 = fig.add_axes(ax2.get_position())
ax5.patch.set_visible(False)
ax5.yaxis.set_visible(False)
ax5.set_yscale('log')
for spinename in ax3.spines:
    if spinename != 'bottom':
        ax3.spines[str(spinename)].set_visible(False)
ax5.spines['bottom'].set_position(('outward', 125))
MLT_conc=[round(x, 1) for x in MLT[start_hope:end_hope]]
p4=ax5.plot(hope_epoch[start_hope:end_hope],  decoy[start_hope:end_hope], 'purple', alpha=0)
spacing = len(MLT_conc)/8
ax5.set_xticks(hope_epoch[start_hope:end_hope][::spacing])
ax5.set_xticklabels(MLT_conc[::spacing])
ax5.set_xlabel('MLT', fontsize=22, fontweight='bold')
ax5.set_xlim(hope_epoch[start_hope], hope_epoch[end_hope] )
    


font = {'family' : 'normal',
                  'weight' : 'bold',
                  'size'   : 22}
plt.rc('font', **font)
fig.set_size_inches(13,9)
plotdir = "WFR_CaseStudies"
if not os.path.exists(plotdir):
      os.umask(0) # unmask if necessary
      os.makedirs(plotdir, 0777) 
os.chdir(plotdir)
plt.savefig(date+'_'+hour1+'_'+hour2+'_WFR_CaseStudy.png')
plt.close()
os.chdir('..')
