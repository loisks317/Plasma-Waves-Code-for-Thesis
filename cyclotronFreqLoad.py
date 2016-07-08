# cyclotronFreqLoad.py
#
# load in the cyclotronFrequency wave amplitudes
#
# LKS, September 2015
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
import matplotlib.pyplot as plt
import spacepy.datamodel as dm
os.chdir('/Users/loisks/Desktop/Functions')
import pickling as pickling
import plots

def depickle(name):
    with open(name, 'rb') as f:
        u = pickle._Unpickler(f)
        u.encoding = 'latin1'
        outf = u.load()
        return outf
# need to modify this once we save the first harmonic
os.chdir('/Users/loisks/Documents/ResearchProjects/PlasmaWaves/WFR_FREQ_240/')
#
# Globals
nLbins=11
nmlt_bins=48

MLT=np.linspace(0.25, 23.75, 48)
totalE=[[[] for i in range(nLbins)] for j in range(nmlt_bins)]

medE=[[ 0 for i in range(nLbins)] for j in range(nmlt_bins)]

#
# load the files and combine
EFiles=glob.glob('ESorted_Cyc_'+'*')
for iFile in range(len(EFiles)):
        dataB=pickling.hdf5_data_open(EFiles[iFile],nLbins, nmlt_bins)
        print(str(iFile) + ' / ' +str(len(EFiles)))
        for iMLT in range(nmlt_bins):
            for iL in range(nLbins):
                try:
                   
                    totalE[iMLT][iL]+=list(np.array(dataB[iMLT][iL]))
                except(OSError):
                    
                    totalE[iMLT][iL]=np.nan
# now get the median
for iMLT in range(nmlt_bins):
        for iL in range(nLbins):
            medE[iMLT][iL]=np.nanmedian(totalE[iMLT][iL])
            #medB[iMLT][iL]=np.nanmedian(totalB[iMLT][iL])


os.chdir('/Users/loisks/Documents/ResearchProjects/PlasmaWaves/CyclotronAmplitudes/')
# save the median
pickle.dump( medE,open('medE_16.p', 'wb'))

# then redo for 6
#medE=pickle.load(open('medE_1.p', 'rb'))
medE6=depickle('medE_6.p')
medE10=depickle('medE_10.p')
            
##
## plot as a fancy plot
#os.chdir('..')
#os.chdir('..')
#plots.fancy_plot(medE,nmlt_bins,nLbins, -14, -10, 'V$^{2}$/m$^{2}$/Hz', 'E', 'low', 'ESpectra_Amplitude_'+resonance, 'CyclotronAmplitudes', 'log', '')
#plots.fancy_plot(medB,nmlt_bins,nLbins, -10, -6, 'nT$^{2}$/m$^{2}$', 'B', 'low', 'BSpectra_Amplitude_'+resonance, 'CyclotronAmplitudes', 'log', '')
#
spe=['P']
n1=0
n2=16
nLbins=11
nmlt_bins=48
lines=[[] for x in range(3)]
MLT=np.linspace(0.25, 23.75, 48)
#
totalDen=[[[] for i in range(nLbins)] for j in range(nmlt_bins)]
medDen=[[ 0 for i in range(nLbins)] for j in range(nmlt_bins)]
#

for iSpe in range(len(spe)):
    os.chdir('/Users/loisks/Documents/ResearchProjects/PlasmaWaves/PlasmaDensity10eV')
    aFiles=glob.glob(spe[iSpe]+'_sat=rbspA_'+'*')
    bFiles=glob.glob(spe[iSpe]+'_sat=rbspB_'+'*')
    for iFile in range(len(aFiles)):
        dataA=pickling.hdf5_data_open(aFiles[iFile],nLbins, nmlt_bins)
        dataB=pickling.hdf5_data_open(bFiles[iFile],nLbins, nmlt_bins)
        for iMLT in range(nmlt_bins):
            for iL in range(nLbins):
                try:
                    totalDen[iMLT][iL]+=list(np.array(dataA[iMLT][iL]))+list(np.array(dataB[iMLT][iL]))
                except:
                    totalDen[iMLT][iL]=np.nan
    # now get the median
    for iMLT in range(nmlt_bins):
        for iL in range(nLbins):
            medDen[iMLT][iL]=np.nanmedian(totalDen[iMLT][iL])
os.chdir('..')

    #
    # now calculate median
    # add the L=2.5 line to each one
#    os.chdir('..')
    #mdensity[[ [] for x in xrange(nLbins)] for x in xrange(nmlt_bins)]
#    plots.fancy_plot(medDen,nmlt_bins,nLbins, -1, 2, 'Density [cm$^{-3}$]', spe[iSpe], 'low', '6_density_'+spe[iSpe], 'DensityPlots10eV', 'log', '')

## plot the lines
#fig=plt.figure()
#ax=fig.add_subplot(111)
#plt.subplots_adjust(right=0.75, top=0.85, bottom=0.15, left=0.17)
Lindex=6
temp=np.swapaxes(medDen,1,0)[Lindex]
tempE=np.swapaxes(medE,1,0)[Lindex]
tempE6=np.swapaxes(medE6,1,0)[Lindex]
tempE10=np.swapaxes(medE10,1,0)[Lindex]
newData=[ 0 for i in range(len(MLT))]
newEData=[ 0 for i in range(len(MLT))]
newEData6=[ 0 for i in range(len(MLT))]
newEData10=[ 0 for i in range(len(MLT))]
for iMLT in range(len(MLT)):
    if MLT[iMLT] < 12:
        #MLT[iMLT]+=12
        newData[iMLT+24]=temp[iMLT]
        newEData[iMLT+24]=tempE[iMLT]
        newEData6[iMLT+24]=tempE6[iMLT]
        newEData10[iMLT+24]=tempE10[iMLT]
    else:
        #MLT[iMLT]-=12
        newData[-24+iMLT]=temp[iMLT]
        newEData[-24+iMLT]=tempE[iMLT]
        newEData6[-24+iMLT]=tempE6[iMLT]
        newEData10[-24+iMLT]=tempE10[iMLT]
#ax.plot(np.array(MLT)[~np.isnan(newData)],np.array( newData)[~np.isnan(newData)], c='r', label='P', lw=2)
#
#
#labels = [item.get_text() for item in ax.get_xticklabels()]
#
#ax.set_yscale('log')
#ax.set_ylabel('Density [cm$^{-3}$]')
#ax.set_xlabel('MLT')
#ax.set_xlim(0,24)
#ax.set_xticklabels(['12', '17', '22', '02', '07'])
#ax2=ax.twinx()
#ax2.plot(np.array(MLT)[~np.isnan(newBData)],np.array( newBData)[~np.isnan(newBData)], c='b', label='Cyclotron', lw=2)
#ax2.set_yscale('log')
#ax2.set_ylabel('nT$^{2}$/m$^{2}$')
#plt.legend()
#font = {'family' : 'normal',
#        'weight' : 'bold',
#        'size'   : 22}
#plt.rc('font', **font)
#os.chdir('DensityPlots10eV')
#plt.savefig('6_linesB_'+resonance+'.pdf')
#plt.close()
#os.chdir('..')
#
#
#
fig=plt.figure()
ax=fig.add_subplot(111)
plt.subplots_adjust(left=0.15, right=0.8, top=0.85, bottom=0.15)
ax.plot(np.array(MLT)[~np.isnan(newData)],np.array( newData)[~np.isnan(newData)], c='Black', ls='--', lw=2)
ax.set_yscale('log')
ax.set_ylabel('Density [cm$^{-3}$]', color='Black')
ax.set_xlabel('MLT')
ax.set_xlim(0,24)
ax.set_xticklabels(['12', '17', '22', '02', '07'])
ax2=ax.twinx()

ax2.plot(np.array(MLT)[~np.isnan(newEData6)],np.array( newEData6)[~np.isnan(newEData6)], c='Blue',  lw=2)
ax2.plot(np.array(MLT)[~np.isnan(newEData10)],np.array( newEData10)[~np.isnan(newEData10)], c='Green',  lw=2)
ax2.plot(np.array(MLT)[~np.isnan(newEData)],np.array( newEData)[~np.isnan(newEData)],  lw=2, c='gold')
plt.legend([ '6 $\Omega_{H^{+}}$', '10 $\Omega_{H^{+}}$', '16 $\Omega_{H^{+}}$'], bbox_to_anchor=[1.3, 1.18], fontsize=20)
ax2.set_yscale('log')
ax2.set_ylabel('V$^{2}$/m$^{2}$/Hz', color='Blue')
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
#plt.legend()
plt.rc('font', **font)
os.chdir('DensityPlots10eV')
plt.savefig('16_linesE_3_Lindex='+str(Lindex)+'.pdf')
plt.close()
os.chdir('..')
