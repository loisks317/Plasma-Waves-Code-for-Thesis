# WFR_sorted_ksolve.py
#
# use bounce time = interaction time = resonance time to solve for k
# and then verify if k vector is reasonable
# still not separating by pitch angle. Figure this out!!
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

import os
import numpy as np
from matplotlib import pyplot as plt
import datetime
from matplotlib.colors import LogNorm
from numpy import ma
from mpl_toolkits.axisartist.grid_finder import FixedLocator, \
     MaxNLocator, DictFormatter
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter 
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
#
# create a half circle for density maps
def dual_half_circle(center, radius, angle=90, ax=None, colors=('w','k'),
                     **kwargs):
    from matplotlib.patches import Wedge
    """
    Add two half circles to the axes *ax* (or the current axes) with the 
    specified facecolors *colors* rotated at *angle* (in degrees).
    """
    if ax is None:
        ax = plt.gca()
    theta1, theta2 = angle, angle + 180
    w1 = Wedge(center, radius, theta1, theta2, fc=colors[0],transform=ax.transData._b, **kwargs)
    w2 = Wedge(center, radius, theta2, theta1, fc=colors[1], transform=ax.transData._b,**kwargs)
    for wedge in [w1, w2]:
        ax.add_artist(wedge)
    return [w1, w2]

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
#lowerB=31
#upperB=36 # to 560 Hz
lowerB=4
upperB=43 # 43 = 1 kHz
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
WIDTHS=[    2.13599992,     2.13599992,     2.13599992,     2.13599992,
            2.13599992,     2.13599992,     2.13599992,     2.13599992,
            2.13599992,     2.13599992,     2.13599992,     2.13599992,
            2.13599992,     4.27266645,     4.27266645,     4.27266645,
            4.27266645,     6.40866661,     6.40866661,     8.54466724,
            8.54466724,     8.54466724,    10.68133259,    10.68133259,
           12.81733322,    14.95333385,    17.09000015,    17.09000015,
           21.36199951,    21.36199951,    27.77066612,    27.77066612,
           32.04333115,    36.31599808,    40.58866501,    46.99733353,
           51.26933289,    57.67799759,    64.08666229,    72.63199615,
           83.31266785,    89.72199249,   102.53933716,   115.35666656,
          130.30999756,   145.26399231,   162.35333252,   181.57933044,
          205.07800293,   230.71266174,   258.48400879,   288.3913269 ,
          324.70733643,   365.2953186 ,   408.02001953,   459.28933716,
          512.69537354,   578.91864014,   647.27801514,   726.31866455,
          816.04003906,   914.30664062,  1027.5267334 ,  1151.42797852,
         1292.41931152]
    #
    # Functions
    # gets the number of months between two dates for reading in the files
date1='20130201' # start date
date2='20150401' # end date


energyRange=[0.5,0.75, 1, 1.5, 2]
paRange=[30,60, 80, 89, 89.9 ]
dt0=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')
def diff_month(d1, d2):
        return (d1.year - d2.year)*12 + d1.month - d2.month
    #
SummedHeating=[[[[ 0 for x in range(Lbins)] for x in range(mltbins)] for y in range(len(energyRange))] for z in range(len(paRange))]

total_month=diff_month(dt2, dt0)
n1=6
q = 1.6*1e-19
mi = 1.67*1e-27
kVecLim=1e-2
allowableTime=3600

allFiles=[[[ np.nan for x in range(Lbins)] for x in range(mltbins)] for j in range(len(FREQUENCIES))]
WaveAmplitudes=[[[ 0 for x in range(Lbins)] for x in range(mltbins)] for j in range(len(FREQUENCIES))]
FreqHeatingRate=[[[[[ 0 for x in range(Lbins)] for x in range(mltbins)] for j in range(len(FREQUENCIES))] for k in range(5)] for m in range(5)]
FreqN=[[[ 0 for x in range(Lbins)] for x in range(mltbins)] for j in range(len(FREQUENCIES))]
for ifreq in range(lowerB,upperB):
     heatingTimes=[[[[ 0 for x in range(Lbins)] for x in range(mltbins)] for y in range(len(energyRange))] for z in range(len(paRange))]
     cover=[[[[ np.nan for x in range(Lbins)] for x in range(mltbins)] for y in range(len(energyRange))] for z in range(len(paRange))]
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
     kVector=[[ 0 for x in range(Lbins)] for x in range(mltbins)]
     for imonth in range(total_month):
        cur_date=str(dt1.month)+'_'+str(dt1.year)
        dt1=dt1+relativedelta(months=1)
        os.chdir('/Users/loisks/Documents/ResearchProjects/PlasmaWaves/Sorted_WFR_Polarization')
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
                     meanfe[imlt][iL]=np.nanmedian(a*(1-b)/2.0) # multiply by the polarization to get the right fraction
                     polDiff[imlt][iL]=np.nanmedian(b/((1-b)/2.0)) # right / left 
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
     plots.fancy_plot(meanfe, mltbins, Lbins, -13, -11, 'V$^{2}$/m$^{2}$ /Hz', 'Wave_Intensity', 'low', 'E_Wave_EMFISIS_wave_intensit_freq='+str(FREQUENCIES[ifreq]), 'Wave_Percent_Polarized', 'log', 'autumn')
          # polarization 
    # plots.fancy_plot(np.array(medPol)*100.0, mltbins, Lbins, 0,90, 'Polarization', 'Wave_Intensity', 'low', 'E_Wave_EMFISIS_Polarization_freq='+str(FREQUENCIES[ifreq]), 'Wave_Percent_Polarized', '', 'viridis')
    # plots.fancy_plot(np.array(perPol)*100.0, mltbins, Lbins, 0,100, 'Percentage [$\%$]', 'Wave_Intensity', 'low', 'E_Wave_EMFISIS_Polarization_Perc_freq='+str(FREQUENCIES[ifreq]), 'Wave_Percent_Polarized', '', 'viridis')
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
#     plt.savefig('linePlot_freq='+str(FREQUENCIES[ifreq])+'.png')
     os.chdir('..')
     
     
     # now use analytic model to show bounce time and heating rate for amplitudes
     LBins=np.linspace(1.5, 4, Lbins)
     ER=6378000 # in meters
     os.chdir('Sorted_BMag')
     Bmed=[[ 0 for i in range(Lbins)] for j in range(mltbins)]
     files=glob.glob('*_MLT_*')
     tempB1=np.array(pickling.hdf5_data_open(files[0], Lbins, mltbins))
     tempB2=np.array(pickling.hdf5_data_open(files[1], Lbins, mltbins))
               #
     # Begin Model Loop
     for iEN in range(len(energyRange)):
           for iPA in range(len(paRange)):
               PA=paRange[iPA]
               particleEnergy=energyRange[iEN]
               CurHz=FREQUENCIES[ifreq]
               #
               # get the median magnetic field data
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
               # print the frequnecy status
               print('freq is: '+str(FREQUENCIES[ifreq]) + ' PA = ' + str(PA) + ' Energy = ' + str(particleEnergy))
               vTot=np.sqrt((2.0*particleEnergy*q)/mi) # for a 1 eV particle
               vperp=np.sqrt((np.sin(PA*np.pi/180.)**2)*(vTot**2))
               vpar=np.sqrt((vTot**2) - (vperp**2))
               
               # LKS unit workout, Nov 2015
               LBinsMod=np.array(np.linspace(1.5, 4, Lbins))#+0.125
               MLTbinsArr=np.linspace(0, 24,mltbins)
               n=np.zeros((mltbins,Lbins))
               t=np.zeros((mltbins,Lbins))
               BounceTime=np.zeros((mltbins,Lbins))
               tB=np.zeros((mltbins,Lbins))
               #heatingTimes=[[ np.nan for i in range(Lbins)] for j in range(mltbins)]
               amplitude=[[ np.nan for i in range(Lbins)] for j in range(mltbins)]
               
               os.chdir('/Users/loisks/Desktop/Functions/')
               
               # this assumes a centered dipole where mag lat = latitude. This isn't
               # quite right
           
               for iL in range(Lbins):
                    for iMLT in range(mltbins):
                        nT= FREQUENCIES[ifreq]/(q* Bmed[iMLT][iL]/(1.0*mi))
                        if nT < 3:
                              # just set to 1
                              n[iMLT][iL]=3
                        else:
                                  n[iMLT][iL]=int(nT)
                        
                        larmor=mi*vperp/(q*Bmed[iMLT][iL])
                        FreqN[ifreq][iMLT][iL]=n[iMLT][iL]
                        #

                        # OVER-RIDE bounce TIME
                        #
                        # I think it's better to account for resonance time instead of bounce time
                        # maybe just assume that k vector is 1e-2
                        # and then calculate a resonance time off of that
                       # BounceTime[iMLT][iL]=(LBins[iL]*ER/vTot)*(3.7-1.6*np.sin(np.pi*PA/180.)) # used PA = 30
                        # define the resonant time as the time it takes for the particles to interact
                        # with the eq noise waves as a fraction of their pitch angle
                        resonantTime=PA/90.
                        BounceTime[iMLT][iL]=allowableTime*PA/90. # 2 hour time period
                        ### TEST THIS WHEN YOU GET BACK!
                        ## AA 
                        AA=np.power(2,(n[iMLT][iL]-1))
                        BB=factorial(n[iMLT][iL]-1)
                        CC=(n[iMLT][iL]-2)
                        DD=(BounceTime[iMLT][iL])*np.power(larmor,n[iMLT][iL]-2)
                        EE= Bmed[iMLT][iL]/np.sqrt(meanfe[iMLT][iL]*WIDTHS[ifreq])
                        #EE= Bmed[iMLT][iL]/(np.sqrt(meanfe[iMLT][iL])*FREQUENCIES[ifreq])
                        #cc=3.0*1e8
                        #EE= Bmed[iMLT][iL]/(meanfe[iMLT][iL]*(WIDTHS[ifreq]* mi * cc**2)/(q*FREQUENCIES[ifreq]**2)) 
                        #EE= Bmed[iMLT][iL]/(meanfe[iMLT][iL]*(FREQUENCIES[ifreq]* mi * cc**2)/(q*WIDTHS[ifreq]**2))       
                        kVector[iMLT][iL]= np.power((AA*BB*EE)/(CC*DD), 1/(1.0*(n[iMLT][iL] - 1)))
                        #
                        # To put back in Kvector limits uncomment this! 
                        if kVector[iMLT][iL] < 0:
                            kVector[iMLT][iL] = np.nan
                            cover[iPA][iEN][iMLT][iL]=1
            
                        condition=(kVector[iMLT][iL]*larmor)**(-2)
                        if condition <= polDiff[iMLT][iL]:
                           cover[iPA][iEN][iMLT][iL]=1
                        
                        elif kVector[iMLT][iL] > kVecLim:
                            cover[iPA][iEN][iMLT][iL]=1

                        # note we divide by a q here to accunt for heating rate in eV!
                        heatingTimes[iPA][iEN][iMLT][iL]=3600.0*meanfe[iMLT][iL]*(q)/(2*mi) 
    
                        # need to check those values
                        if cover[iPA][iEN][iMLT][iL] != 1:
                            FreqHeatingRate[iPA][iEN][ifreq][iMLT][iL]=heatingTimes[iPA][iEN][iMLT][iL]
               print("energy is: " +str(particleEnergy) + " Kvector is: " + str(kVector[10][2]) )            
               os.chdir('/Users/loisks/Documents/ResearchProjects/PlasmaWaves/')
               ####
                #
      
                #       modMedTb[ii][jj]=tB[ii][jj]
               os.chdir('/Users/loisks/Documents/ResearchProjects/PlasmaWaves/')
#               plots.fancy_plot(BounceTime, mltbins, Lbins, -1, 2, 'Bounce Time', 'SimpleModel_Results', 'low', str(CurHz)+'Hz_En=1_PA='+str(PA), 'SimpleModel_Results/EN='+str(particleEnergy)+'_PA='+str(PA), 'log', 'PiYG')
#               os.chdir('..')
#               plots.fancy_plot(Bmed, mltbins, Lbins, -4, -7, '|B| [T]', 'SimpleModel_Results', 'low', 'Bmag', 'SimpleModel_Results/EN='+str(particleEnergy)+'_PA='+str(PA), 'log', 'autumn')
#               os.chdir('..')
#               plots.fancy_plot(kVector, mltbins, Lbins, -4, -1, '|k|', 'KVectorMaps', 'low', str(CurHz)+'Hz_En='+str(particleEnergy)+'_PA='+str(PA), 'KVectorMaps/EN='+str(particleEnergy)+'_PA='+str(PA), 'log', 'autumn')
#               os.chdir('..')
#               plots.fancy_plot(n, mltbins, Lbins, 0, 2, 'Harmonic Number', 'n', 'low', str(CurHz)+'Hz_En='+str(particleEnergy)+'_PA='+str(PA), 'SimpleModel_Results/EN='+str(particleEnergy)+'_PA='+str(PA), 'log', 'autumn')
#               os.chdir('..')
               #
               # HEATING TIMES PLOT
               #
               # WORK ON GETTING THIS TO GRAY SCALE ON RETURN
               #
               #

               directory='HeatingTimes'
               label=str(CurHz)+'Hz_En='+str(iEN)+'_PA='+str(PA)
               scale='log'
               vmin=-2; vmax=0; colorbar_Label='eV / Hour'; cmap='autumn'
               mlts=(np.linspace(0, 24, mltbins))*(15*np.pi/180.)
               L_shells=np.linspace(1.5,4 , Lbins)
               # intitialize the figure
               fig=plt.figure()
               ax=fig.add_subplot(111, polar=True)
               datah_m=ma.masked_invalid(np.array(heatingTimes[iPA][iEN]).transpose())
               data2=ma.masked_invalid(np.array(cover[iPA][iEN]).transpose())
               X,Y=np.meshgrid(mlts, L_shells)
               ax.set_ylim(0,4)
          
               dual_half_circle((0,0), 1.0, angle=90, ax=None, colors=('w','k'))
               cbaxes = fig.add_axes([0.75, 0.15, 0.03, 0.75])
               m2 = np.linspace(0,1,Lbins*mltbins).reshape(Lbins,mltbins) 
               if scale=='log':
                   col=ax.pcolormesh( X, Y, np.log10(datah_m), cmap=cmap, vmin=vmin, vmax=vmax )
                   p=ax.pcolormesh( X, Y, data2, cmap='Greys', alpha=0.6, vmin=0, vmax=2 )
                   #for ik,jk in zip(p.get_facecolors(),m2.flatten()):
                   #    ik[3] = 0.5 # Set the alpha value of the RGBA tuple using m2
                   cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(vmin, vmax+1))
               else:
                   col=ax.pcolormesh( X, Y,datah_m, vmin=vmin,cmap=cmap, vmax=vmax)
                   cb = plt.colorbar(col, cax = cbaxes,ticks=range(vmin, vmax+1))
               plt.subplots_adjust(right=0.73, top=0.8, bottom=0.15)
               #xT2=['', '2', '', '4', '', '6']
               xT2=['', '1', '', '2', '', '3']
               xL2=['00','', '06','', '12','', '18']
               cb.set_label(colorbar_Label, fontsize=30)
               ax.tick_params(axis='both', which='major', labelsize=25)    
               cb.ax.tick_params(labelsize=35) 
               ax.set_yticklabels(xT2, fontsize=30)
               ax.set_xticklabels(xL2, fontsize=30)
               ax.grid(True, lw=4)
               plt.draw()
               subdir_name=directory
               if not os.path.exists(subdir_name):
                   os.umask(0) # unmask if necessary
                   os.makedirs(subdir_name) 
               os.chdir(subdir_name)#
               fig.set_size_inches(13,9)
#               plt.savefig(label+'.png')
               plt.close(fig)
               os.chdir('..')
               
               #
               #
               #
               # make a line plot at L = 2.25 or so
               import numpy.ma as ma
               fig=plt.figure()
               ax=fig.add_subplot(111)
               plt.subplots_adjust(right=0.7, top=0.9, bottom=0.15, left=.2)
               ax.set_yscale('log')
               ax.set_xlim(0,24)
               font = {'family' : 'normal',
                       'weight' : 'bold',
                       'size'   : 22}
               plt.rc('font', **font)
               
               ax.set_ylabel('Resonance Time/'+'\n'+'Bounce Period')
               colors=['seagreen','darkviolet', 'DeepPink', 'Orange', 'Turquoise']
               colorsLight=['lightcoral', 'lightgreen', 'violet', 'navajowhite', 'paleturquoise']
               lshells=[2,3,4,5,6]
               #
               # set up the heating rate
                         # plot the figure 
               #ax2=fig.add_subplot(212)
               #plt.subplots_adjust(right=0.75, top=0.9, bottom=0.15, left=.2)
               ax.set_yscale('log')
               ax.set_xlim(0,24)
               font = {'family' : 'normal',
                       'weight' : 'bold',
                       'size'   : 20}
               plt.rc('font', **font)
               ax.set_xlabel('MLT')
               ax.set_ylabel('eV/hour')          
               mlt=np.linspace(0,24,mltbins)
               parray=[]
               itemCount=0
              
               for item in range(Lbins):
                    temp=np.array(np.swapaxes(heatingTimes[iPA][iEN],1,0)[item])
                    temp2=np.array(np.swapaxes(cover[iPA][iEN],1,0)[item])
                    temp3=np.array(np.swapaxes(amplitude,1,0)[item])
                    
                    HeatingRate=np.zeros(mltbins)
                    Covers=np.zeros(mltbins)
                    CoverValues=np.zeros(mltbins)
                   
                    for imlt in range(mltbins):
                      if imlt+24 < mltbins:
                        HeatingRate[imlt+24]=temp[imlt]
                        CoverValues[imlt+24]=temp2[imlt]*temp[imlt]
                        Covers[imlt+24]=temp2[imlt]
                        
                        if Covers[imlt+24] !=1 :
                            SummedHeating[iPA][iEN][imlt+24][item]+=HeatingRate[imlt+24]
                            WaveAmplitudes[ifreq][imlt+24][item]=temp3[imlt]
                      else:
                        HeatingRate[imlt-24]=temp[imlt]
                        Covers[imlt-24]=temp2[imlt]
                        CoverValues[imlt-24]=temp2[imlt]*temp[imlt]
                        
                        if Covers[imlt-24] != 1:
                            SummedHeating[iPA][iEN][imlt-24][item]+=HeatingRate[imlt-24]
                            WaveAmplitudes[ifreq][imlt-24][item]=temp3[imlt]
                        
                    # plot it
                    if item in lshells:
                        parray += ax.plot(mlt,HeatingRate, lw=3, c=colors[itemCount])
                        ax.plot(mlt, CoverValues, lw=3, c='White')
                        ax.plot(mlt, CoverValues, lw=2,ls='--', c=colors[itemCount])
                        itemCount+=1

               a=ax.get_xticks().tolist()
               for ilist in range(len(a)):
                   a[ilist]+=12
                   if a[ilist] >= 24:
                       a[ilist]=int(a[ilist]-24)
                   else:
                       a[ilist]=int(a[ilist])
               #ax.set_xticklabels([])
               ax.set_xticklabels(a)
               ax.set_ylim(1e-3, 1e-0)
               ax.axvspan(13,16, alpha=0.1, color='blue')
               ax.axvline(13, lw=2, ls='--',c='Silver')
               ax.axvline(16, lw=2, ls='--',c='Silver')
               plt.legend(parray,['L=2.0', "L=2.25","L=2.5", "L=2.75", "L=3.0"],bbox_to_anchor=[1.58, 0.7], fontsize=20,scatterpoints=1)
               #os.chdir('HeatingRateLinePlots')
               subdir_name='Freq='+str(FREQUENCIES[ifreq])
               #if not os.path.exists(subdir_name):
               #    os.umask(0) # unmask if necessary
               #    os.makedirs(subdir_name) 
               #os.chdir(subdir_name)#
               #plt.figure(figsize=(13,9))
#               plt.savefig('linePlot_freq='+str(FREQUENCIES[ifreq])+'_PA='+str(PA)+'_EN='+str(particleEnergy)+'.png')
               os.chdir('..')
               os.chdir('..')
                
                #
                ### FIx WHY HEATING RATES ARE DIFFERENT
     print('Summed Heating by energy' + str(SummedHeating[1][:][0][0]))                                                     
     totalDen=[[[] for i in range(Lbins)] for j in range(mltbins)]
     medDen=[[ 0 for i in range(Lbins)] for j in range(mltbins)]
     PDensity=[[ np.nan for i in range(Lbins)] for j in range(mltbins)]
     os.chdir('/Users/loisks/Documents/ResearchProjects/PlasmaWaves/PlasmaDensity10eV')
     aFiles=glob.glob('*'+'A_species=P')
     bFiles=glob.glob('*'+'B_species=P')
     for iFile in range(len(aFiles)):
          dataA=pickle.load(open(aFiles[iFile], 'rb'))
          dataB=pickle.load(open(bFiles[iFile], 'rb'))
          for iMLT in range(mltbins):
              for iL in range(Lbins):
                  try:
                      totalDen[iMLT][iL]+=list(np.array(dataA[iMLT][iL]))+list(np.array(dataB[iMLT][iL]))
                  except:
                      totalDen[iMLT][iL]=np.nan
      # now get the median
     for iMLT in range(mltbins):
          for iL in range(Lbins):            
                  medDen[iMLT][iL]=np.nanmedian(totalDen[iMLT][iL])
     os.chdir('..')
     for imlt in range(mltbins):
       for iL in range(Lbins):
          if imlt+24 < mltbins:
              PDensity[imlt+24][iL]=medDen[imlt][iL]
          else:
              PDensity[imlt-24][iL]=medDen[imlt][iL]
      
      
     
     import numpy.ma as ma
     fig=plt.figure()
     ax=fig.add_subplot(111)
     plt.subplots_adjust(right=0.57, top=0.9, bottom=0.5, left=.15)
     ax.set_yscale('log')
     ax.set_xlim(0,24)
     font = {'family' : 'normal',
             'weight' : 'bold',
             'size'   : 22}
     plt.rc('font', **font)
     
     colors=['seagreen','darkviolet', 'DeepPink', 'Orange', 'Turquoise']
     colorsLight=['lightcoral', 'lightgreen', 'violet', 'navajowhite', 'paleturquoise']
     lshells=[2,3,4,5,6]
     #
     # set up the heating rate
               # plot the figure 
     #ax2=fig.add_subplot(212)
     #plt.subplots_adjust(right=0.7, top=0.9, bottom=0.15, left=.2)
     ax.set_yscale('log')
     ax.set_xlim(0,24)
     font = {'family' : 'normal',
             'weight' : 'bold',
             'size'   : 20}
     plt.rc('font', **font)
     ax.set_xlabel('MLT', fontsize=18)
     ax.set_ylabel('eV/hour/Hz', fontsize=18)          
     mlt=np.linspace(0,24,mltbins)
     parray=[]
     itemCounter=0
     labels=[]
     labelsList=['L=2.0', "L=2.25","L=2.5", "L=2.75", "L=3.0", "H$^{+}$ Density"]
     for item in lshells:
         tempArr=np.swapaxes(SummedHeating[-2][2],0,1)[item]
         tempArr[tempArr<1e-2]=np.nan
         try:
             parray += ax.plot(mlt,tempArr, lw=3, c=colors[itemCounter])
             labels.append(labelsList[itemCounter])
         except:
             print("no heating yet")
             print tempArr
         itemCounter+=1
     labels.append(labelsList[-1])
     a=ax.get_xticks().tolist()
     for ilist in range(len(a)):
         a[ilist]+=12
         if a[ilist] >= 24:
             a[ilist]=int(a[ilist]-24)
         else:
             a[ilist]=int(a[ilist])
     ax.set_xticklabels(a)
     ax.set_ylim(1e-2, 1e1)
     ax.axvspan(13,16, alpha=0.1, color='blue')
     ax.axvline(13, lw=2, ls='--',c='Silver')
     ax.axvline(16, lw=2, ls='--',c='Silver')
     
     
     ax2=ax.twinx()
     parray+= ax2.plot(mlt,np.swapaxes(PDensity,1,0)[4], c='Black', ls='--', lw=2)
     ax2.set_yscale('log')
     ax2.set_ylabel('Density [cm$^{-3}$]', color='Black', fontsize=18)
     ax2.set_xticklabels(a)
     
     plt.legend(parray,labels,bbox_to_anchor=[2.05, 0.9], fontsize=16,scatterpoints=1)
     
     os.chdir('HeatingRateLinePlots')
     subdir_name='Freq='+str(FREQUENCIES[ifreq])
     if not os.path.exists(subdir_name):
           os.umask(0) # unmask if necessary
           os.makedirs(subdir_name) 
     os.chdir(subdir_name)#
     #plt.figure(figsize=(13,9))
     os.chdir('..')
#     plt.savefig('IntegratedLinePlot_freq='+str(FREQUENCIES[ifreq])+'_PA=89_EN=1eV_time='+str(allowableTime)+'.png')
#     os.chdir('..')
     
     Heating2=[ [ 0 for i in range(Lbins)] for j in range(mltbins)]
     for imlt in range(mltbins):
       for iL in range(Lbins):
          if imlt+24 > mltbins-1:
              Heating2[imlt-24][iL]=SummedHeating[-1][1][imlt][iL]
          else:
              Heating2[imlt+24][iL]=SummedHeating[-1][1][imlt][iL]
     # line plot of Heating rate by frequency and n
     MBINS=np.linspace(0,24,mltbins)
     LBINS=np.linspace(1.5, 4, Lbins)
 
     plots.log_box_plot(np.swapaxes(FreqHeatingRate[-2][2][lowerB:upperB],0,2)[4], MBINS,FREQUENCIES[lowerB:upperB], 0,24, 1, 1000, 'MLT', 'Frequency [Hz]',-2,0, 'eV/Hour/Hz', '', 'HeatingRateMap/PA='+str(89)+'_EN=1eV', 'FreqHeatingRate', cmap='autumn',yscale='log')
     os.chdir('..')      
#
# print out the Heating rate
print('Heating rate at low energy: ' + str(heatingTimes[1][0][5][5]))
print('Heating rate at high energy: ' + str(heatingTimes[1][4][5][5]))

    
#plots.fancy_plot(Heating2, mltbins, Lbins, -1, 1, 'eV/Hour', 'SummedHeatingMaps', 'low', 'PA='+str(89), 'SummedHeatingMaps/PA='+str(89)+'_EN=1', 'log', 'autumn')
#os.chdir('..')
 
 # 
 #

     
#plots.log_box_plot(np.swapaxes(FreqN[lowerB:upperB],0,1)[10],FREQUENCIES[lowerB:upperB],LBINS,  FREQUENCIES[lowerB], FREQUENCIES[upperB-1],1.5, 4,  'Frequency [Hz]','L-Shell', 0, 1, 'n$^{th}$ Harmonic', '', 'HeatingRateMap/PA='+str(PA)+'_EN='+str(particleEnergy), 'NHeatingRate', cmap='autumn', yscale='linear')
#os.chdir('..')

#
#
# now plot in terms of PA and Energy at L = 2.5, bin=3 and 1 eV

fig=plt.figure()
ax=fig.add_subplot(111)
plt.subplots_adjust(right=0.57, top=0.9, bottom=0.5, left=.15)
ax.set_yscale('log')
ax.set_xlim(0,24)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)

colors=['seagreen','darkviolet', 'DeepPink', 'Orange', 'Turquoise']
colorsLight=['lightcoral', 'lightgreen', 'violet', 'navajowhite', 'paleturquoise']
lshells=[2,3,4,5,6]

ax.set_yscale('log')
ax.set_xlim(0,24)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 20}
plt.rc('font', **font)
ax.set_xlabel('MLT', fontsize=18)
ax.set_ylabel('eV/hour', fontsize=18)          
mlt=np.linspace(0,24,mltbins)
parray=[]
legLab=[]
itemCounter=0
PADATA=np.swapaxes(SummedHeating, 0, 3)[4][1]
for item in range(len(paRange)):
    temp=list(np.swapaxes(SummedHeating[item][2],1,0)[4])
    temp2=np.array(temp)
    temp2[temp2<1e-2]=np.nan
    parray += ax.plot(mlt,temp2, lw=3, c=colors[itemCounter])
    legLab.append( 'PA = ' + str(paRange[item]) +'$^{o}$')
    itemCounter+=1
a=ax.get_xticks().tolist()
for ilist in range(len(a)):
    a[ilist]+=12
    if a[ilist] >= 24:
        a[ilist]=int(a[ilist]-24)
    else:
        a[ilist]=int(a[ilist])
ax.set_xticklabels(a)
ax.set_ylim(1e-2, 1e1)
ax.axvspan(13,16, alpha=0.1, color='blue')
ax.axvline(13, lw=2, ls='--',c='Silver')
ax.axvline(16, lw=2, ls='--',c='Silver')

ax2=ax.twinx()
parray+= ax2.plot(mlt,np.swapaxes(PDensity,1,0)[4], c='Black', ls='--', lw=2)
ax2.set_yscale('log')
ax2.set_ylabel('Density [cm$^{-3}$]', color='Black', fontsize=18)
ax2.set_xticklabels(a)
legLab.append("H$^{+}$ Density")

plt.legend(parray,legLab,bbox_to_anchor=[2.05, 0.9], fontsize=16,scatterpoints=1)

os.chdir('HeatingRateLinePlots')
subdir_name='All'
if not os.path.exists(subdir_name):
      os.umask(0) # unmask if necessary
      os.makedirs(subdir_name) 
os.chdir(subdir_name)#
#plt.figure(figsize=(13,9))
#plt.savefig('IntegratedLinePlot_PA_time='+str(allowableTime)+'.png')
#os.chdir('..')
#os.chdir('..')



#
# NOW ENERGY 
#

fig=plt.figure()
ax=fig.add_subplot(111)
plt.subplots_adjust(right=0.57, top=0.9, bottom=0.5, left=.15)
ax.set_yscale('log')
ax.set_xlim(0,24)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)

colors=['seagreen','darkviolet', 'DeepPink', 'Orange', 'Turquoise']
colorsLight=['lightcoral', 'lightgreen', 'violet', 'navajowhite', 'paleturquoise']
ax.set_yscale('log')
ax.set_xlim(0,24)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 20}
plt.rc('font', **font)
ax.set_xlabel('MLT', fontsize=18)
ax.set_ylabel('eV/hour', fontsize=18)          
mlt=np.linspace(0,24,mltbins)
parray=[]
legLab=[]
itemCounter=0
for item in range(len(energyRange)):
    temp=list(np.swapaxes(SummedHeating[-2][item],1,0)[4])
    temp2=np.array(temp)
    temp2[temp2<1e-2]=np.nan
    parray += ax.plot(mlt,temp2, lw=3, c=colors[itemCounter])
    legLab.append( str(energyRange[item]) +' eV')
    itemCounter+=1
a=ax.get_xticks().tolist()
for ilist in range(len(a)):
    a[ilist]+=12
    if a[ilist] >= 24:
        a[ilist]=int(a[ilist]-24)
    else:
        a[ilist]=int(a[ilist])
ax.set_xticklabels(a)
ax.set_ylim(1e-2, 1e1)
ax.axvspan(13,16, alpha=0.1, color='blue')
ax.axvline(13, lw=2, ls='--',c='Silver')
ax.axvline(16, lw=2, ls='--',c='Silver')

ax2=ax.twinx()
parray+= ax2.plot(mlt,np.swapaxes(PDensity,1,0)[4], c='Black', ls='--', lw=2)
ax2.set_yscale('log')
ax2.set_ylabel('Density [cm$^{-3}$]', color='Black', fontsize=18)
ax2.set_xticklabels(a)
legLab.append("H$^{+}$ Density")

plt.legend(parray,legLab,bbox_to_anchor=[2.05, 0.9], fontsize=16,scatterpoints=1)

os.chdir('HeatingRateLinePlots')
subdir_name='All'
if not os.path.exists(subdir_name):
      os.umask(0) # unmask if necessary
      os.makedirs(subdir_name) 
os.chdir(subdir_name)#
#plt.figure(figsize=(13,9))
#plt.savefig('IntegratedLinePlot_ENERGY_time='+str(allowableTime)+'.png')
#os.chdir('..')
#os.chdir('..')


for iL in range(2, 6):
  for iEN in range(5):
      for iPA in range(5):

          ww=np.swapaxes(SummedHeating[iPA][iEN],1,0)[iL]
          pp=np.swapaxes(PDensity,1,0)[iL]
          a=np.corrcoef(np.log10(ww[~np.isnan(pp)]),np.log10(pp[~np.isnan(pp)]))
          print('PA = '+str(paRange[iPA])+' L = ' + str(lshells[iL-2]) + ' Energy = '+str(energyRange[iEN]) + '\n' +
                'R = '+ str(a[0][1]))


# from 50 to 315, r = 0.69


#
#
# plot each L density and PSD separately
itemCounter=0
colors=['seagreen','darkviolet', 'DeepPink', 'Orange', 'Turquoise']
colorsLight=['lightcoral', 'lightgreen', 'violet', 'navajowhite', 'paleturquoise']
for iL in range(2, 7):

    fig=plt.figure()
    ax=fig.add_subplot(111)
    plt.subplots_adjust(right=0.57, top=0.9, bottom=0.5, left=.15)
    ax.set_yscale('log')
    ax.set_xlim(0,24)
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 22}
    plt.rc('font', **font)
    

    ax.set_yscale('log')
    ax.set_xlim(0,24)
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 20}
    plt.rc('font', **font)
    ax.set_xlabel('MLT', fontsize=18)
    ax.set_ylabel('eV/hour', fontsize=18)          
    mlt=np.linspace(0,24,mltbins)
    parray=[]
    legLab=[]

    temp=list(np.swapaxes(SummedHeating[2][-2],1,0)[iL])
    temp2=np.array(temp)
    temp2[temp2<1e-2]=np.nan
    parray = ax.plot(mlt,temp2, lw=3, c=colors[itemCounter])
    itemCounter+=1       

    a=ax.get_xticks().tolist()
    for ilist in range(len(a)):
        a[ilist]+=12
        if a[ilist] >= 24:
            a[ilist]=int(a[ilist]-24)
        else:
            a[ilist]=int(a[ilist])
    ax.set_xticklabels(a)
    ax.set_ylim(1e-2, 1e1)
    ax.axvspan(13,16, alpha=0.1, color='blue')
    ax.axvline(13, lw=2, ls='--',c='Silver')
    ax.axvline(16, lw=2, ls='--',c='Silver')
    
    ax2=ax.twinx()
    parray+= ax2.plot(mlt,np.swapaxes(PDensity,1,0)[iL], c='Black', ls='--', lw=2)
    ax2.set_yscale('log')
    ax2.set_ylabel('Density [cm$^{-3}$]', color='Black', fontsize=18)
    ax2.set_xticklabels(a)
    legLab.append("H$^{+}$ Density")
    
    plt.legend(parray,legLab,bbox_to_anchor=[2.05, 0.9], fontsize=16,scatterpoints=1)
    
    os.chdir('HeatingRateLinePlots')
    subdir_name='DensityLCompare'
    if not os.path.exists(subdir_name):
          os.umask(0) # unmask if necessary
          os.makedirs(subdir_name) 
    os.chdir(subdir_name)#
    #plt.figure(figsize=(13,9))
#    plt.savefig('IntegratedLinePlot_L='+str(iL)+'_time='+str(allowableTime)+'.png')
#    os.chdir('..')
#    os.chdir('..')

# shift to get correlation coefficient
papa=np.swapaxes(SummedHeating[-2][2], 1,0)[4]
newArr=np.zeros(len(papa))
newArr[0:-2]=papa[2:]
newArr[-2:]=papa[0:2]
a=np.corrcoef(np.log10(newArr[~np.isnan(papa)]),np.log10(papa[~np.isnan(papa)]))
print(a)
