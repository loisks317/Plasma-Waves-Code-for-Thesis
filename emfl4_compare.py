# emfl4_compare.py
#
# just look at the fraction of events from 614 days of L4 processed data
#
# LKS, March 2016, Tromso, Norway. 
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

# globals

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

import pickle
os.chdir('PolarizationL4results')
LeftHand=pickle.load(open('leftHandPol.p', 'rb'))
RightHand=pickle.load(open('rightHandPol.p', 'rb'))
Linear=pickle.load(open('linearPol.p', 'rb'))
os.chdir('..')
#
# now go into line plot
f=plt.figure()
ax=f.add_subplot(111)
plt.subplots_adjust(right=0.7, top=0.9, bottom=0.15)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
ax.plot(FREQUENCIES,RightHand, lw=3, c='DarkOrange')
ax.plot(FREQUENCIES,LeftHand, lw=3, c='b')
ax.plot(FREQUENCIES,Linear, lw=3, c='turquoise')
ax.set_xlim(10,10000)
ax.set_xscale('log')
ax.axvline(x=150, c='k', lw=2, ls='--')
ax.axvline(x=600, c='k', lw=2, ls='--')
ax.set_ylabel('Number of Events')
ax.set_xlabel('Frequencies [Hz]')
ax.set_yscale('log')
ax.set_title('Polarization between 2 < L < 3')


plt.legend(['Right Hand', 'Left Hand', 'Linear'], bbox_to_anchor=[1.3, 0.7])
subdir_name='Comparison_Pol'
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name) 
os.chdir(subdir_name)#
f.set_size_inches(13,9)
plt.savefig('polarization_L4.pdf')
plt.close(f)
os.chdir('..')
