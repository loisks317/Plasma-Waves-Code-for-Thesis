# loadPlasmaDensity.py
#
# load in the combined plasma density
# and plot everything
#
# LKS August 2015 # Last day of August :)
# editted September for new ion moments
#
# import list 
import os
import glob 
import numpy as np
from matplotlib import pyplot as plt
import datetime
import pickle
import matplotlib.dates as dates
os.chdir('/Users/loisks/Desktop/Functions/')
import pickling
import plots
os.chdir('/Users/loisks/Documents/ResearchProjects/PlasmaWaves/')
#
#
spe=['P']
n1=0
n2=16
nLbins=11
nmlt_bins=48
lines=[[] for x in range(3)]
MLT=np.linspace(0.25, 23.75, 48)
FREQUENCIES_WFR=[  2.13599992e+00,   4.27199984e+00,   6.40899992e+00,
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

totalDen=[[[] for i in range(nLbins)] for j in range(nmlt_bins)]
medDen=[[ 0 for i in range(nLbins)] for j in range(nmlt_bins)]
os.chdir('/Users/loisks/Documents/ResearchProjects/Misc/stat_corr/GyroFreqEmfisis/')
with open('meanAll.pickle', 'rb') as handle:
           meanMedian=pickle.load(handle)
os.chdir('/Users/loisks/Documents/ResearchProjects/PlasmaWaves/')

for iSpe in range(len(spe)):
    os.chdir('PlasmaDensity10eV')
    aFiles=glob.glob(spe[iSpe]+'_sat=rbspA_'+'*')
    bFiles=glob.glob(spe[iSpe]+'_sat=rbspB_'+'*')
    for iFile in range(len(aFiles)):
        dataA=pickling.hdf5_data_open(aFiles[iFile],nLbins, nmlt_bins)
        dataB=pickling.hdf5_data_open(bFiles[iFile],nLbins, nmlt_bins)
        for iMLT in range(nmlt_bins):
            for iL in range(nLbins):
                try:
                    totalDen[iMLT][iL]+=list(np.array(dataA[iMLT][iL]))+list(np.array(dataB[iMLT][iL]))
                except(OSError):
                    totalDen[iMLT][iL]=np.nan
    # now get the median
    for iMLT in range(nmlt_bins):
        for iL in range(nLbins):
            medDen[iMLT][iL]=np.nanmedian(totalDen[iMLT][iL])


    #
    # now calculate median
    # add the L=2.5 line to each one
    os.chdir('..')
    #mdensity[[ [] for x in xrange(nLbins)] for x in xrange(nmlt_bins)]
    plots.fancy_plot(medDen,nmlt_bins,nLbins, -1, 2, 'Density [cm$^{-3}$]', spe[iSpe], 'low', 'density_'+spe[iSpe], 'DensityPlots10eV', 'log', '')

# plot the lines
fig=plt.figure()
ax=fig.add_subplot(111)
plt.subplots_adjust(left=0.15, right=0.8, top=0.85, bottom=0.15)
temp=np.swapaxes(medDen,1,0)[3] # L of 2
freq1=np.where(np.array(FREQUENCIES_WFR)>300)[0][0]
freq2=np.where(np.array(FREQUENCIES_WFR)>50)[0][0]
freq3=np.where(np.array(FREQUENCIES_WFR)>100)[0][0]



tempW1=np.swapaxes(meanMedian, 1, 0)[freq1]
tempW2=np.swapaxes(meanMedian,1, 0)[freq2]
tempWAll=np.array(tempW2)
for iFreq in range(freq2+1, freq3+1):
    tempWAll+=np.array(np.swapaxes(meanMedian,1, 0)[iFreq])
tempW3=np.swapaxes(meanMedian,1,0)[freq3]
medWaves1=[ 0 for j in range(len(MLT))]
medWaves2=[ 0 for j in range(len(MLT))]
medWaves3=[ 0 for j in range(len(MLT))]
medWavesAll=[ 0 for j in range(len(MLT))]
newData=[ 0 for i in range(len(MLT))]
for iMLT in range(len(MLT)):
    if MLT[iMLT] < 12:
        #MLT[iMLT]+=12
        newData[iMLT+24]=temp[iMLT]
        medWaves1[iMLT+24]=tempW1[iMLT]
        medWaves2[iMLT+24]=tempW2[iMLT]
        medWaves3[iMLT+24]=tempW3[iMLT]
        medWavesAll[iMLT+24]=tempWAll[iMLT]
    else:
        #MLT[iMLT]-=12
        newData[-24+iMLT]=temp[iMLT]
        medWaves1[iMLT-24]=tempW1[iMLT]
        medWaves2[iMLT-24]=tempW2[iMLT]
        medWaves3[iMLT-24]=tempW3[iMLT]
        medWavesAll[iMLT-24]=tempWAll[iMLT]
ax.plot(np.array(MLT)[~np.isnan(newData)],np.array( newData)[~np.isnan(newData)], c='k', label='P', lw=3, ls='--')
labels = [item.get_text() for item in ax.get_xticklabels()]

#ax.plot(MLT, lines[1], c='b', label='He', lw=2)
#ax.plot(MLT, lines[2], c='g', label='O', lw=2)
#plt.legend(bbox_to_anchor=[1.35, 0.5])
ax.set_yscale('log')
ax.set_ylabel('Density [cm$^{-3}$]')
ax.set_xlabel('MLT')
ax.set_xlim(0,24)
ax.set_xticklabels(['12', '17', '22', '02', '07'])

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
axl2=ax.twinx()
axl2.set_yscale('log')
#axl2.set_ylim(1e-12, 1e-7)
axl2.set_xlim(0,24)
axl2.set_ylabel('V$^{2}$ m$^{-2}$ Hz$^{-1}$', fontweight='bold', color='r')
axl2.plot(MLT, medWavesAll,lw=3, color='r') 
#freq1=np.where(np.array(FREQUENCIES_WFR)>300)[0][0]
#freq2=np.where(np.array(FREQUENCIES_WFR)>50)[0][0]
#freq3=np.where(np.array(FREQUENCIES_WFR)>100)[0][0]
#axl2.plot(MLT,medWaves2 , lw=3, color='r')
#axl2.plot(MLT, medWaves3, lw=3, color='b')
#axl2.plot(MLT, medWaves1, lw=3, color='g')
#plt.legend(['50 Hz', '100 Hz', '300 Hz'], bbox_to_anchor=[1.3, 1.18], fontsize=20)
font = {'family' : 'normal',
                  'weight' : 'bold',
                  'size'   : 22}
plt.rc('font', **font)





os.chdir('DensityPlots10eV')
plt.savefig('lines.pdf')
plt.close()
os.chdir('..')
