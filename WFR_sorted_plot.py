# WFR_sorted_plot.py
#
# plot the output of the EMFISIS HFR data in a fancy plot by each frequency
#
# LKS, January 2015
#
#
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
def factorial(n):
    if n == 0:
        return 1
    else:
        return n * factorial(n-1)
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
dt0=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')
def diff_month(d1, d2):
        return (d1.year - d2.year)*12 + d1.month - d2.month
    #

total_month=diff_month(dt2, dt0)
n1=6
q = 1.6*1e-19
m = 1.67*1e-27
allFiles=[[[ np.nan for x in range(Lbins)] for x in range(mltbins)] for j in range(len(FREQUENCIES))]
for ifreq in range(27,43):
 dt1=dt0 
    #
    # loop through each month
    # find the bad data month and remove it - I think august 2014
 cfilee=[[ [] for x in range(Lbins)] for x in range(mltbins)]
 cPol=[[ [] for x in range(Lbins)] for x in range(mltbins)]
 meanfe=[[ [] for x in range(Lbins)] for x in range(mltbins)]
 medPol=[[ [] for x in range(Lbins)] for x in range(mltbins)]
 polDiff=[[ [] for x in range(Lbins)] for x in range(mltbins)]
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

    
    #
    # okay now exclude anamolies 
 for imlt in range(mltbins):
    for iL in range(Lbins):
        try:
            a=np.array(cfilee[imlt][iL])
            b=np.array(cPol[imlt][iL])
            # screen for polarization
            #noGo=np.where(np.array(cPol[imlt][iL] > 0.2))[0]
            good=len(b[b<=0.2])
            if len(b) != 0:
                perPol[imlt][iL]=good/(1.0*len(b[~np.isnan(b)]))
            else:
                perPol[imlt][iL]=np.nan
            #a[b>0.2]=np.nan # nan the ones with the wrong polarization
            medPol[imlt][iL]=np.nanmedian(b)
            meanfe[imlt][iL]=np.nanmedian(a*(1-b)) # multiply by the polarization to get the right fraction
            polDiff[imlt][iL]=np.nanmedian(b/(1-b)) # right / left 
            allFiles[ifreq][imlt][iL]=np.nanmedian(a)
                          #  except(TypeError, AttributeError):
        except(IndexError):
            perPol[imlt][iL]=np.nan
            medPol[imlt][iL]=np.nan
            meanfe[imlt][iL]=np.nan
            polDiff[imlt][iL]=np.nan
            allFiles[ifreq][imlt][iL]=np.nan
    
#
    # Now plotstop
 plots.fancy_plot(meanfe, mltbins, Lbins, -14, -10, 'V$^{2}$/m$^{2}$ /Hz', 'Wave_Intensity', 'low', 'E_Wave_EMFISIS_wave_intensit_freq='+str(FREQUENCIES[ifreq]), 'Wave_Percent_Polarized', 'log', 'viridis')
 # polarization 
 plots.fancy_plot(np.array(medPol)*100.0, mltbins, Lbins, 0,90, 'Polarization', 'Wave_Intensity', 'low', 'E_Wave_EMFISIS_Polarization_freq='+str(FREQUENCIES[ifreq]), 'Wave_Percent_Polarized', '', 'viridis')
 plots.fancy_plot(np.array(perPol)*100.0, mltbins, Lbins, 0,100, 'Percentage [$\%$]', 'Wave_Intensity', 'low', 'E_Wave_EMFISIS_Polarization_Perc_freq='+str(FREQUENCIES[ifreq]), 'Wave_Percent_Polarized', '', 'viridis')
 #
 # make a line plot at L = 2.25 or so
 fig=plt.figure()
 ax=fig.add_subplot(111)
 plt.subplots_adjust(right=0.75, top=0.9, bottom=0.15, left=.2)
 ax.set_yscale('log')
 ax.set_xlim(0,24)
 font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
 plt.rc('font', **font)
 ax.set_xlabel('MLT')
 ax.set_ylabel('V$^{2}$/m$^{2}$/Hz')
 colors=['r','b', 'g', 'Orange']
 lshells=[0,1,2,3]
 mlt=np.linspace(0,24,mltbins)
 for item in range(len(lshells)):
     temp=np.swapaxes(meanfe,1,0)[lshells[item]]
     newData=np.zeros(mltbins)
     for imlt in range(mltbins):
         if imlt+24 < mltbins:
             newData[imlt+24]=temp[imlt]
         else:
             newData[imlt-24]=temp[imlt]
            
     ax.plot(mlt,newData, lw=2, c=colors[item]) 
 a=ax.get_xticks().tolist()
 for ilist in range(len(a)):
     a[ilist]+=12
     if a[ilist] >= 24:
         a[ilist]=int(a[ilist]-24)
     else:
         a[ilist]=int(a[ilist])
 ax.set_xticklabels(a)
 ax.axvspan(13,16, alpha=0.1, color='black')
 ax.axvline(13,lw=2, ls='--', c='k')
 ax.axvline(16,lw=2, ls='--', c='k')
 plt.legend(['L=1.5', "L=2.0", "L=2.5", "L=3.0"],bbox_to_anchor=[1.43, 0.75], fontsize=20)
 os.chdir('Wave_Percent_Polarized')
 plt.savefig('linePlot_freq='+str(FREQUENCIES[ifreq])+'.png')
 os.chdir('..')


 # now use analytic model to show bounce time and heating rate for amplitudes
 LBins=np.linspace(1.5, 4, Lbins)
 ER=6378000 # in meters
 mi=1.67*1e-27
 q=1.6*1e-19
 #kperpArr=[1e-2,5*1e-3,1e-3, 1e-4]
 kperpArr=[4*1e-3,3*1e-3, 2*1e-3, 1e-3, 1e-4]
 DegTrap=[89.9,89,88,80,70,60,50,40,30,20,10]
 CurHz=FREQUENCIES[ifreq]
 #
 # get the median magnetic field data
 os.chdir('Sorted_BMag')
 Bmed=[[ 0 for i in range(Lbins)] for j in range(mltbins)]
 files=glob.glob('*_MLT_*')
 tempB1=np.array(pickling.hdf5_data_open(files[0], Lbins, mltbins))
 tempB2=np.array(pickling.hdf5_data_open(files[1], Lbins, mltbins))
 #
 for iMLT in range(mltbins):
  for iL in range(Lbins):
   try:
     temp=np.array(list(tempB1[iMLT][iL])+list(tempB2[iMLT][iL]))*1e-9
     temp[temp<0]=np.nan
     temp[temp>1e-4]=np.nan
     Bmed[iMLT][iL]=np.nanmedian(temp)
   except(OSError):
     Bmed[iMLT][iL]=np.nan
 os.chdir('/Users/loisks/Documents/ResearchProjects/PlasmaWaves/')
 for iK in range(len(kperpArr)):
     for iDeg in range(len(DegTrap)):
          kperp=kperpArr[iK]
          PA=DegTrap[iDeg]
          vTot=np.sqrt((0.5*q)/mi) # for a 1 eV particle
          vperp=np.sqrt((np.sin(PA*np.pi/180.)**2)*(vTot**2))
          print vperp
          vpar=np.sqrt((vTot**2) - (vperp**2))
          sqrtE=np.sqrt(np.array(meanfe)*CurHz)# need to multiply by  Hz in here
          # LKS unit workout, Nov 2015
          LBinsMod=np.array(np.linspace(1.5, 4, Lbins))#+0.125
          MLTbinsArr=np.linspace(0, 24,mltbins)
          n=np.zeros((mltbins,Lbins))
          t=np.zeros((mltbins,Lbins))
          BounceTime=np.zeros((mltbins,Lbins))
          tB=np.zeros((mltbins,Lbins))
          os.chdir('/Users/loisks/Desktop/Functions/')
          # this assumes a centered dipole where mag lat = latitude. This isn't
          # quite right
          for iL in range(Lbins):
              for iMLT in range(mltbins):
                  nT=int(2*np.pi*CurHz*mi/(q*Bmed[iMLT][iL]))
                  if nT < 3:
                        # just set to 1
                        n[iMLT][iL]=3
                  else:
                            n[iMLT][iL]=int(nT)
                  gamma=kperp*sqrtE[iMLT][iL]/Bmed[iMLT][iL]
                  larmor=mi*vperp/(q*Bmed[imlt][iL])
                  #
                  # check and make sure we can have resonance
                  condition=(kperp*larmor)**(-2)
                  if condition <= polDiff[imlt][iL]:
                      print 'Violation!'
                      print FREQUENCIES[ifreq]
                      print kperp
                      print PA
                  Ro=kperp*vperp*mi/(q*Bmed[iMLT][iL])
                  t[iMLT][iL]=(1.0/gamma)*(np.power(2,(n[iMLT][iL]-1)))*factorial(n[iMLT][iL]-1)/((n[iMLT][iL]-2)*(np.power(Ro,(n[iMLT][iL]-2))))
                  BounceTime[iMLT][iL]=(LBins[iL]*ER/vpar)*(3.7-1.6*np.sin(np.pi*PA/180.)) # used PA = 30
                  tB[iMLT][iL]=t[iMLT][iL]/BounceTime[iMLT][iL]
          # need to check those values
          os.chdir('/Users/loisks/Documents/ResearchProjects/PlasmaWaves/')
          #
          # now plot log of time for these
          modMedT=[[ np.nan for i in range(Lbins)] for j in range(mltbins)]
          modMedTb=[[ np.nan for i in range(Lbins)] for j in range(mltbins)]
          for ii in range(mltbins):
             for jj in range(Lbins):
                 modMedT[ii][jj]=t[ii][jj]
                 modMedTb[ii][jj]=tB[ii][jj]
          os.chdir('/Users/loisks/Documents/ResearchProjects/PlasmaWaves/')
          plots.fancy_plot(modMedTb, mltbins, Lbins, -1, 2, 'Resonance Time/Bounce Period', 'SimpleModel_Results', 'low', str(CurHz)+'Hz_En=1_K='+str(iK)+'_PA='+str(PA), 'SimpleModel_Results', 'log', 'PiYG')
          plots.fancy_plot(Bmed, mltbins, Lbins, -4, -7, '|B| [T]', 'SimpleModel_Results', 'low', 'Bmag', 'SimpleModel_Results', 'log', 'viridis') 
          # make a line plot at L = 2.25 or so
          fig=plt.figure()
          ax=fig.add_subplot(211)
          #plt.subplots_adjust(right=0.75, top=0.9, bottom=0.15, left=.2)
          ax.set_yscale('log')
          ax.set_xlim(0,24)
          font = {'family' : 'normal',
                  'weight' : 'bold',
                  'size'   : 22}
          plt.rc('font', **font)
          
          ax.set_ylabel('Resonance Time/'+'\n'+'Bounce Period')
          colors=['r','b', 'g', 'Orange']
          lshells=[0,2,4,6]
          #
          # set up the heating rate
                    # plot the figure 
          ax2=fig.add_subplot(212)
          plt.subplots_adjust(right=0.75, top=0.9, bottom=0.15, left=.2)
          ax2.set_yscale('log')
          ax2.set_xlim(0,24)
          font = {'family' : 'normal',
                  'weight' : 'bold',
                  'size'   : 20}
          plt.rc('font', **font)
          ax2.set_xlabel('MLT')
          ax2.set_ylabel('dW [eV]')          
          mlt=np.linspace(0,24,mltbins)
          for item in range(len(lshells)):
              temp=np.swapaxes(modMedTb,1,0)[lshells[item]]
              BounceTimes=np.zeros(mltbins)
              for imlt in range(mltbins):
                if imlt+24 < mltbins:
                  BounceTimes[imlt+24]=temp[imlt]
                  heatingRate=BounceTimes*(1/2.)*(q**2)*np.array(newData)/(m*q)
                else:
                  BounceTimes[imlt-24]=temp[imlt]
                  heatingRate=BounceTimes*(1/2.)*(q**2)*np.array(newData)/(m*q)

                  # need to make it plot within here on separate axes for bounce time
              ax.plot(mlt,BounceTimes, lw=2, c=colors[item])
              ax2.plot(mlt,heatingRate, lw=2, c=colors[item])
          a=ax.get_xticks().tolist()
          for ilist in range(len(a)):
              a[ilist]+=12
              if a[ilist] >= 24:
                  a[ilist]=int(a[ilist]-24)
              else:
                  a[ilist]=int(a[ilist])
          ax.set_xticklabels([])
          ax2.set_xticklabels(a)
          ax.axvspan(13,16, alpha=0.1, color='black')
          ax.axvline(13,lw=2, ls='--', c='k')
          ax.axvline(16,lw=2, ls='--', c='k')
          ax2.axvspan(13,16, alpha=0.1, color='black')
          ax2.axvline(13,lw=2, ls='--', c='k')
          ax2.axvline(16,lw=2, ls='--', c='k')
          plt.legend(['L=1.5', "L=1.75", "L=2.0", "L=2.25"],bbox_to_anchor=[1.45, 1.5], fontsize=20)
          os.chdir('BounceTimeLinePlots')
          subdir_name='Freq='+str(FREQUENCIES[ifreq])+'/k='+str(kperpArr[iK])
          if not os.path.exists(subdir_name):
              os.umask(0) # unmask if necessary
              os.makedirs(subdir_name) 
          os.chdir(subdir_name)#
          #plt.figure(figsize=(13,9))
          plt.savefig('linePlot_freq='+str(FREQUENCIES[ifreq])+'_K='+str(iK)+'_PA='+str(PA)+'.png')
          os.chdir('..')
          os.chdir('..')
          os.chdir('..')
          #
          ###

         
#
#dataF=np.swapaxes(allFiles, 0,2)[4] # L = 2, MLT x ENERGY
#subdir = "/Users/loisks/Documents/ResearchProjects/Misc/stat_corr/GyroFreqEmfisis"
## switch and get the hybrids
#os.chdir(subdir)
#gyroIon=depickle( 'gyroIon.pickle')
#gyroEle=depickle('gyroEle.pickle')
#
#os.chdir('/Users/loisks/Documents/ResearchProjects/PlasmaWaves/')
#gyroIon=(np.swapaxes(gyroIon, 1, 0)[10])/(2*np.pi)
#gyroEle=(np.swapaxes(gyroEle, 1, 0)[10])/(2*np.pi)
#lowerhybrid=np.sqrt(gyroIon*gyroEle)
##
## now go back
#MLT_bins=np.linspace(0,24,mltbins)
#xlabel='MLT'
#ylabel='Frequency [Hz]'
#from numpy import ma
## intitialize the figure
#fig=plt.figure()
#ax1=fig.add_subplot(111)
#ax1.set_xlabel(xlabel, fontsize=30, fontweight='bold')
#ax1.set_ylabel(ylabel, fontsize=30, fontweight='bold')
#plt.subplots_adjust(right=0.8, top=0.9, bottom=0.15)
## mask the nans
## log mods here 
#data=np.log10(np.array(dataF))
#datah_m=ma.masked_invalid(data.transpose())
#X,Y=np.meshgrid(np.linspace(0,24,mltbins), FREQUENCIES)
#colorbar_max=-11
#colorbar_min=-14
#col=ax1.pcolormesh(X,Y,datah_m,cmap='viridis', vmax=colorbar_max, vmin=colorbar_min)
## colorbar_min=int(np.log10(colorbar_min))
## colorbar_max=int(np.log10(colorbar_max))+1
#
## col=ax.pcolormesh(X,Y,datah_m,cmap='jet')
#xmin= 0
#xmax = 24
#ymin=10
#ymax=10000
#ax1.set_xlim(xmin,xmax)
#ax1.set_ylim(ymin,ymax)
#ax1.set_yscale('log')
#cbaxes = fig.add_axes([0.81, 0.15, 0.03, 0.75]) 
#cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(colorbar_min, colorbar_max+1) ) 
#ax1.plot(MLT_bins, gyroIon[0:-1]*6, color='k', lw=6)
#ax1.plot(MLT_bins, gyroIon[0:-1]*6, color='silver', lw=5)
#ax1.plot(MLT_bins, lowerhybrid[0:-1], color='k', lw=6)
#ax1.plot(MLT_bins, lowerhybrid[0:-1], color='white', lw=5)
#cb.set_label('V$^{2}$/ m$^{2}$ / Hz', fontsize=30, fontweight='bold')
#ax1.tick_params(axis='both', which='major', labelsize=30)
#cb.ax.tick_params(labelsize=30) 
#ax1.set_title('', fontsize=25)
#subdir_name='All_Freqs'
#plt.draw()
#if not os.path.exists(subdir_name):
#    os.umask(0) # unmask if necessary
#    os.makedirs(subdir_name) 
#os.chdir(subdir_name)#
#fig.set_size_inches(13,9)
#plt.savefig('All_Freqs.png')
#plt.close(fig)
#os.chdir('..')
#
#
#    
#
#
#
#
## normalize by mlt
#normE=[[ [] for x in range(len(FREQUENCIES))] for x in range(mltbins)]
#for ifreq in range(len(FREQUENCIES)):
#        d=np.swapaxes(dataF, 1,0)[ifreq]
#        max=np.nanmax(d)
#        for imlt in range(mltbins):
#           normE[imlt][ifreq]=dataF[imlt][ifreq]/(1.0*max)
#        
##plots.log_box_plot(normE, np.linspace(0,24,mltbins), FREQUENCIES, 0, 24,10, 10000, 'MLT', 'Frequency [Hz]',-2,0 ,'' ,'Normalized Power Spectral Density', 'E', 'All_Freqs', 'Normalized_All_Freqs')
#fig=plt.figure()
#ax1=fig.add_subplot(111)
#ax1.set_xlabel(xlabel, fontsize=30, fontweight='bold')
#ax1.set_ylabel(ylabel, fontsize=30, fontweight='bold')
#plt.subplots_adjust(right=0.8, top=0.9, bottom=0.15)
## mask the nans
## log mods here 
#data=np.log10(np.array(normE))
#datah_m=ma.masked_invalid(data.transpose())
#X,Y=np.meshgrid(np.linspace(0,24,mltbins), FREQUENCIES)
#colorbar_max=0
#colorbar_min=-2
#col=ax1.pcolormesh(X,Y,datah_m,cmap='viridis', vmax=colorbar_max, vmin=colorbar_min)
## colorbar_min=int(np.log10(colorbar_min))
## colorbar_max=int(np.log10(colorbar_max))+1
#
## col=ax.pcolormesh(X,Y,datah_m,cmap='jet')
#xmin= 0
#xmax = 24
#ymin=10
#ymax=10000
#ax1.set_xlim(xmin,xmax)
#ax1.set_ylim(ymin,ymax)
#ax1.set_yscale('log')
#cbaxes = fig.add_axes([0.81, 0.15, 0.03, 0.75]) 
#cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(colorbar_min, colorbar_max+1) ) 
#ax1.plot(MLT_bins, gyroIon[0:-1]*6, color='k', lw=6)
#ax1.plot(MLT_bins, gyroIon[0:-1]*6, color='silver', lw=5)
#ax1.plot(MLT_bins, lowerhybrid[0:-1], color='k', lw=6)
#ax1.plot(MLT_bins, lowerhybrid[0:-1], color='white', lw=5)
#cb.set_label('V$^{2}$/ m$^{2}$ / Hz', fontsize=30, fontweight='bold')
#ax1.tick_params(axis='both', which='major', labelsize=30)
#cb.ax.tick_params(labelsize=30) 
#ax1.set_title('', fontsize=25)
#subdir_name='All_Freqs'
#plt.draw()
#if not os.path.exists(subdir_name):
#    os.umask(0) # unmask if necessary
#    os.makedirs(subdir_name) 
#os.chdir(subdir_name)#
#fig.set_size_inches(13,9)
#plt.savefig('Normalized_All_Freqs.png')
#plt.close(fig)
#os.chdir('..')
   
