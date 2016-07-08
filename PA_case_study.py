# PA_case_study.py
#
# PA spectogram case study
#
# LKS, March, 2015
#
# 
import glob
import matplotlib.pyplot as plt
from spacepy import pycdf
import numpy as np
import os
import datetime
from matplotlib.dates import HourLocator, DayLocator, DateFormatter, MinuteLocator
import matplotlib.collections as collections
import matplotlib.dates as dt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter 
from numpy import ma
import pandas as pd
from collections import deque,Counter
from bisect import insort, bisect_left
import generateMatrices as gM
os.chdir('/Users/loisks/Documents/Functions/')
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
os.chdir('/Users/loisks/Desktop/ResearchProjects/PlasmaWaves/')
from itertools import islice
days = DayLocator(interval=1) 
hours = MinuteLocator(interval=15) 
#hours2 = HourLocator(interval=3) 
hours2=MinuteLocator(interval=30)
daysFmt = DateFormatter('%H:%M')
#
# define an energy threshold
energyThres=2
# 2 ev +/- .2 
#
#
def RunningMedian(seq, M):
    """
     Purpose: Find the median for the points in a sliding window (odd number in size) 
              as it is moved from left to right by one point at a time.
      Inputs:
            seq -- list containing items for which a running median (in a sliding window) 
                   is to be calculated
              M -- number of items in window (window size) -- must be an integer > 1
      Otputs:
         medians -- list of medians with size N - M + 1
       Note:
         1. The median of a finite list of numbers is the "center" value when this list
            is sorted in ascending order. 
         2. If M is an even number the two elements in the window that
            are close to the center are averaged to give the median (this
            is not by definition)
    """   
    seq = iter(seq)
    s = []   
    m = M // 2

    # Set up list s (to be sorted) and load deque with first window of seq
    s = [item for item in islice(seq,M)]    
    d = deque(s)

    # Simple lambda function to handle even/odd window sizes    
    median = lambda : s[m] if bool(M&1) else (s[m-1]+s[m])*0.5

    # Sort it in increasing order and extract the median ("center" of the sorted window)
    s.sort()    
    medians = [median()]   

    # Now slide the window by one point to the right for each new position (each pass through 
    # the loop). Stop when the item in the right end of the deque contains the last item in seq
    for item in seq:
        old = d.popleft()          # pop oldest from left
        d.append(item)             # push newest in from right
        del s[bisect_left(s, old)] # locate insertion point and then remove old 
        insort(s, item)            # insert newest such that new sort is not required        
        medians.append(median())  
    return medians
#
#
#  Pitch Angles
PA_labels=[4.5, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0, 175.5]


PA_labels=[4.5, 18.0, 36.0, 54.0, 72.0, 90.0, 108.0, 126.0, 144.0, 162.0, 175.5]
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
#
# do loop
#date1='20130605'
date1='20130702'
date2='20130703'
dt1=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')
fchan=31 # like 250 Hz
n1=0
n2=10
while dt1 != dt2:
        #
        date=datetime.datetime.strftime(dt1,'%Y%m%d')
        frequency=str(FREQUENCIES[fchan]) + ' Hz'
        os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EMFISIS_WFR_A')
        emf_file=glob.glob('rbsp-a_WFR-spectral-matrix_emfisis-L2_'+date+'*.cdf')
        emf_pyf=pycdf.CDF(emf_file[0])
        emf_epoch=emf_pyf['Epoch'][...]
        eeePOCH=pd.DatetimeIndex(emf_epoch)
        EuEu=np.array(emf_pyf['EuEu'][...])
        EvEv=np.array(emf_pyf['EvEv'][...])
        EwEw=np.array(emf_pyf['EwEw'][...])
        EPSD=EuEu+EvEv+EwEw
        #
        # need to pandas the HOPE and EFW data
        rt = pd.period_range(date,periods=1440, freq='T').to_timestamp()
        rng=rt[::1]
        #
        #
        # get the ellipticity
        os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EMFISIS_L4_A')
        file=glob.glob('rbsp-a_wna-survey_emfisis-L4_20130702_v1.5.3.cdf')
        pyfL4=pycdf.CDF(file[0])
        ellsvd=pyfL4['ellsvd'][...]
        Lsvd=pyfL4['L'][...]
        MLTsvd=pyfL4['MLT'][...]
        emfL4epoch=pyfL4['Epoch'][...]
        #
        os.chdir('/Users/loisks/Desktop/liemohn10/loisks/HOPE_A')
        hope_file=glob.glob('rbspa_rel02_ect-hope-PA-L3_'+date+'*.cdf')
        hope_pyf=pycdf.CDF(hope_file[0])
        hope_epoch=hope_pyf['Epoch_Ion'][...]
        hope_Hflux=hope_pyf['FPDU'][...]
        hope_energy=hope_pyf['HOPE_ENERGY_Ion'][...]
        L=hope_pyf['L_Ion'][...]        
        MLT=hope_pyf['MLT_Ion'][...]

        # get the spacecraft potential
        os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EFW_L3_A')
        f=glob.glob('*'+date+'*')
        pyfe=pycdf.CDF(f[0])
        potential=np.abs(pyfe['Vavg'][...])
        potential[potential>5]=np.nan

        #
        pEpoch=pd.DatetimeIndex(pyfe['epoch'][...])
        df=pd.DataFrame(potential, index=pEpoch, columns=['potential'])
        Phi=np.array(df['potential'].resample('1min', how='median').reindex(rng, fill_value=np.nan))

        # HOPE data
        hEpoch=pd.DatetimeIndex(hope_epoch)
        PCfluxes=[ np.zeros(len(rng)) for j in range(11) ]
        LDF=pd.DataFrame(L, index=hEpoch, columns=['L'])
        MLTDF=pd.DataFrame(MLT,index=hEpoch,columns=['MLT'])
        Lpd=np.array(LDF['L'].resample('1min', how='median').reindex(rng, fill_value=np.nan))
        MLTpd=np.array(MLTDF['MLT'].resample('1min', how='median').reindex(rng, fill_value=np.nan))
        
        # should have several itmes
        for iEn in range(n1,n2):
          hopeEnergyPD=pd.DataFrame(np.swapaxes(hope_energy,1,0)[iEn], index=hEpoch,columns=['energy'])
          energyH=np.array(hopeEnergyPD['energy'].resample('1min', how='median').reindex(rng, fill_value=np.nan))
          energyH[np.isnan(energyH)]=0
          correctedEN=energyH+Phi
          # keep only energies within +/- .2 of the energy threshold
          #correctedEN[correctedEN< energyThresh-.2]=np.nan
          #correctedEN[correctedEN> energyThresh+.2]=np.nan
                    
          for iPA in range(11):
            screen=np.swapaxes(hope_Hflux,2,0)[iEn][iPA]
            screen[screen<=0]=np.nan
            fluxDF=pd.DataFrame(screen, index=hEpoch, columns=['flux']) # this leave time vs PA
            fluxPD=np.array(fluxDF['flux'].resample('1min',how='mean').reindex(rng,fill_value=0))
            fluxPD[np.isnan(fluxPD)]=0
            #
            # now we have to spacecraft potential correct
            # sum at the end across all energy channels
            fluxPD[correctedEN<(energyThres-0.3)]=0
            fluxPD[correctedEN>(energyThres+0.3)]=0
            #
            # check
            PCfluxes[iPA][np.isnan(PCfluxes[iPA])]=0
            PCfluxes[iPA]+=np.array(fluxPD)
            # to get the total fluxes within this range
        
        #
        # change Epochs 
    #
    # change Epochs
        fig=plt.figure(figsize=(20,20))
        #fig.gca().xaxis.set_major_locator(hours)
        #fig.gca().xaxis.set_major_formatter(daysFmt)
        #fig.gca().xaxis.set_minor_locator(hours2)
        font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 22}
        plt.rc('font', **font)
        #
        # first do the magnetic field
        ax0=fig.add_subplot(311)
       # iEMFSpectra=np.array(iEMFSpectra)
        plt.subplots_adjust(right=0.1+0.65, top=0.95, bottom=0.5,left=0.1)
        ax0.set_ylabel('Frequency [Hz]', fontsize=22, fontweight='bold')
        time=dt.date2num(emf_epoch)
        X,Y=np.meshgrid(time, FREQUENCIES)
        #a[a<=0]=np.nan
        data2=ma.masked_invalid(np.log10(EPSD))
        vmin=-13
        vmax=-10
        ax0.set_xlim(dt.date2num(datetime.datetime(2013,7,2,9,0,0)),dt.date2num(datetime.datetime(2013,7,2,11,30,0)) )
        ax0.set_yscale('log')
        ax0.set_ylim(1e2,800)
        #ax0.set_xlabel("Time", fontsize=20, fontweight='bold')
        col=ax0.pcolormesh(X,Y,data2.transpose(), cmap='viridis', vmin=vmin, vmax=vmax)
        ax0.set_xlim(dt.date2num(datetime.datetime(2013,7,2,9,0,0)),dt.date2num(datetime.datetime(2013,7,2,11,30,0)) )
        font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 22}
        plt.rc('font', **font)
        cbaxes = fig.add_axes([0.78, 0.82, 0.03, 0.13]) 
        cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(vmin, vmax+1))
        fig.gca().xaxis.set_major_locator(hours)
        fig.gca().xaxis.set_major_formatter(daysFmt)
        fig.gca().xaxis.set_minor_locator(hours2)
        ax0.set_xlim(dt.date2num(datetime.datetime(2013,7,2,9,0,0)),dt.date2num(datetime.datetime(2013,7,2,11,30,0)) )        
        cb.set_label('V$^{2}$/m$^{2}$/Hz', fontsize=22, fontweight='bold')
        ax0.tick_params(axis='both', which='major', labelsize=22)
        ax0.set_xticklabels([])
        cb.ax.tick_params(labelsize=30)        
        #a= Lpd <= 3
        #b= Lpd >= 2
        #mlta=MLTpd >= 1
        #mltb=MLTpd <= 4
        #overlap=np.zeros(len(Lpd))
        #for io in range(len(overlap)):
        #    if ((a[io]==True) & (b[io]==True) & (mlta[io]==True) & (mltb[io]==True)):
        #        overlap[io]=True
        #    else:
        #        overlap[io]=False
        #collection_Lpd = collections.BrokenBarHCollection.span_where(
        #       dt.date2num(rng.to_pydatetime()), ymin=0.0, ymax=1e12, where=overlap, facecolor='black', alpha=0.25)
        #ax0.add_collection(collection_Lpd)
##################
        axpol=fig.add_subplot(312)
        axpol.set_xlim(dt.date2num(datetime.datetime(2013,7,2,9,0,0)),dt.date2num(datetime.datetime(2013,7,2,11,30,0)) )
        axpol.plot(RunningMedian(dt.date2num(emfL4epoch), 100),RunningMedian(ellsvd[:,23],100), lw=2, color='DarkOrange')        
        axpol.plot(RunningMedian(dt.date2num(emfL4epoch), 100),RunningMedian(ellsvd[:,27],100), lw=2, color='DeepPink')
        axpol.plot(RunningMedian(dt.date2num(emfL4epoch), 100),RunningMedian(ellsvd[:,31],100), lw=2, color='Blue')
        axpol.plot(RunningMedian(dt.date2num(emfL4epoch), 100),RunningMedian(ellsvd[:,35],100), lw=2, color='Turquoise')
        axpol.legend(['f = 100 Hz', 'f = 160 Hz', 'f = 250 Hz', 'f = 400 Hz'], bbox_to_anchor=[1.32, 1])
        axpol.set_ylim(-1,1.05)
        axpol.plot(RunningMedian(dt.date2num(emfL4epoch), 100), np.zeros(len(RunningMedian(dt.date2num(emfL4epoch), 100))), lw=2, ls='--', c='k')
        axpol.set_ylabel('Ellipticity', fontweight='bold')
        axpol.set_xticklabels([])
        axpol.set_xlim(dt.date2num(datetime.datetime(2013,7,2,9,0,0)),dt.date2num(datetime.datetime(2013,7,2,11,30,0)) )
        #a= np.array(RunningMedian(np.array(Lsvd),100)) <= 3
        #b= np.array(RunningMedian(np.array(Lsvd),100)) >= 2
        #mlta=np.array(RunningMedian(np.array(MLTsvd),100)) >= 1
        #mltb=np.array(RunningMedian(np.array(MLTsvd),100)) <= 4
        #overlapc=np.zeros(len(RunningMedian(Lsvd,100)))
        #for io in range(len(overlapc)):
        #    if ((a[io]==True) & (b[io]==True) & (mlta[io]==True) & (mltb[io]==True)):
        #        overlapc[io]=True
        #    else:
        #        overlapc[io]=False
        #collection_Lsvd = collections.BrokenBarHCollection.span_where(
        #    RunningMedian(dt.date2num(emfL4epoch),100), ymin=-100, ymax=100, where=overlapc, facecolor='black', alpha=0.25)
        #axpol.add_collection(collection_Lsvd)
        axpol.set_xlim(dt.date2num(datetime.datetime(2013,7,2,9,0,0)),dt.date2num(datetime.datetime(2013,7,2,11,30,0)) )
##################

        ax1=fig.add_subplot(313)
        #plt.subplots_adjust(right=0.1+0.65, top=0.95, bottom=0.5,left=0.1)
        ax1.set_ylabel('PA', fontsize=22, fontweight='bold')

        time=dt.date2num(rng.to_pydatetime())
        X,Y=np.meshgrid(time, PA_labels)
        a=np.array(PCfluxes)
        a[a<=0]=np.nan
        data=ma.masked_invalid(np.log10(a)).transpose()
        vmin=6
        vmax=10
        ax1.set_xlim(dt.date2num(datetime.datetime(2013,7,2,9,0,0)),dt.date2num(datetime.datetime(2013,7,2,11,30,0)) )
        ax1.set_xlabel("Time", fontsize=20, fontweight='bold')
        col=ax1.pcolormesh(X,Y,data.transpose(), cmap='viridis', vmin=vmin, vmax=vmax)
        ax1.set_xlim(dt.date2num(datetime.datetime(2013,7,2,9,0,0)),dt.date2num(datetime.datetime(2013,7,2,11,30,0)) )
        font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 22}
        plt.rc('font', **font)
        cbaxes = fig.add_axes([0.84, 0.5, 0.03, 0.14]) 
        cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(vmin, vmax+1)) 
        cb.set_label('cm$^{-2}$ s$^{-1}$ sr$^{-1}$ keV$^{-1}$', fontsize=22, fontweight='bold')
        ax1.tick_params(axis='both', which='major', labelsize=22)
        cb.ax.tick_params(labelsize=30) 
        # plot EMFISIS data
        ax2 = ax1.twinx()
        ax2.set_yscale('log')
        ax2.set_ylim(9*1e-14, 1e-10)
        ax2.set_xlim(dt.date2num(datetime.datetime(2013,7,2,9,0,0)),dt.date2num(datetime.datetime(2013,7,2,11,0,0)) )
        #ax2.set_xlim(hope_epoch[start_hope], hope_epoch[end_hope] )
        ax2.set_ylabel('V$^{2}$/m$^{2}$/Hz', fontsize=22, fontweight='bold')
        ax2.plot(RunningMedian(dt.date2num(emf_epoch), 100), RunningMedian(EPSD[:,fchan], 100),'-', lw=3, c='k')

        fig.gca().xaxis.set_major_locator(hours)
        fig.gca().xaxis.set_major_formatter(daysFmt)
        fig.gca().xaxis.set_minor_locator(hours2)
        ax2.set_xlim(dt.date2num(datetime.datetime(2013,7,2,9,0,0)),dt.date2num(datetime.datetime(2013,7,2,11,30,0)) )
    #
    #
        #a= L <= 3
        #b= L >= 2
        ##mlta=MLT >= 1
        ##mltb=MLT <= 4
        #overlap=np.zeros(len(L))
        #for io in range(len(overlap)):
        #    if ((a[io]==True) & (b[io]==True) ):
        #        overlap[io]=True
        #    else:
        #        overlap[io]=False
        #collection_L = collections.BrokenBarHCollection.span_where(
        #       dt.date2num(hope_epoch), ymin=0.0, ymax=1e12, where=overlap, facecolor='black', alpha=0.25)
        #ax2.add_collection(collection_L)
        
        # Add L LABEL
        ax3 = fig.add_axes(ax2.get_position())
        ax3.patch.set_visible(False)
        ax3.yaxis.set_visible(False)
        ax3.set_yscale('log')
        ax3.set_xlim(dt.date2num(datetime.datetime(2013,7,2,9,0,0)),dt.date2num(datetime.datetime(2013,7,2,11,30,0)) )
        for spinename in ax3.spines:
            if spinename != 'bottom':
                ax3.spines[str(spinename)].set_visible(False)
        ax3.spines['bottom'].set_position(('outward', 65))

        BegP=np.where(hope_epoch > datetime.datetime(2013,7,2,9,0,0))[0][0]
        EndP=np.where(hope_epoch > datetime.datetime(2013,7,2,11,30,0))[0][0]
        L_conc=[round(x, 2) for x in L[BegP:EndP]]
        decoy=data[:,5]
        p3=ax3.plot(time,  decoy, 'g', alpha=0)
        ax3.set_xlim(dt.date2num(datetime.datetime(2013,7,2,9,0,0)),dt.date2num(datetime.datetime(2013,7,2,11,30,0)) )
        # I hate doing this but
        spacing = int(len(L_conc)/6)
        ax3.set_xticks(hope_epoch[BegP:EndP][::spacing])
        ax3.set_xticklabels(L_conc[::spacing])    
        ax3.set_xlabel('L-Shell', fontsize=22, fontweight='bold')
        
        # Add MLT Label
        ax5 = fig.add_axes(ax2.get_position())
        ax5.patch.set_visible(False)
        ax5.yaxis.set_visible(False)
        ax5.set_yscale('log')
        ax5.set_xlim(dt.date2num(datetime.datetime(2013,7,2,9,0,0)),dt.date2num(datetime.datetime(2013,7,2,11,30,0)) )
        for spinename in ax3.spines:
            if spinename != 'bottom':
                ax3.spines[str(spinename)].set_visible(False)
        ax5.spines['bottom'].set_position(('outward', 125))
        MLT_conc=[round(x, 1) for x in MLT[BegP:EndP]]
        p4=ax5.plot(time,  decoy, 'purple', alpha=0)
        ax5.set_xlim(dt.date2num(datetime.datetime(2013,7,2,9,0,0)),dt.date2num(datetime.datetime(2013,7,2,11,30,0)) )
        spacing = int(len(MLT_conc)/6)
        
        ax5.set_xticks(hope_epoch[BegP:EndP][::spacing])
        ax5.set_xticklabels(MLT_conc[::spacing])
        ax5.set_xlabel('MLT', fontsize=22, fontweight='bold')
        font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 22}
        plt.rc('font', **font)
        #fig.set_size_inches(13,9)
        os.chdir('/Users/loisks/Desktop/ResearchProjects/PlasmaWaves/')
        subdir_name='Polarized_Partial_Case_Studies_WFR'
        if not os.path.exists(subdir_name):
            os.umask(0) # unmask if necessary
            os.makedirs(subdir_name) 
        os.chdir(subdir_name)#
        plt.savefig(date+'_Polarized_PA_case_FullDay_study.png')
        os.chdir('..')
        plt.close()
        date=datetime.datetime.strftime((datetime.datetime.strptime(date,'%Y%m%d')+datetime.timedelta(days=1)), '%Y%m%d')
        print(dt1)
        dt1+=datetime.timedelta(days=1)
