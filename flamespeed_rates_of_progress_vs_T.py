# -*- coding: utf-8 -*-
"""
Created on Sat Feb 03 16:33:19 2018

@author: HP USER
"""

import matplotlib.pyplot as plt
import numpy as np


plt.plot(results[0].solution['T'],results[0].net_rates_of_progress[176,:],'b')
plt.xlabel('Temperature (K)')
plt.ylabel(r'Net rates of progress $\omega$')


plt.plot(results[11].solution['T'],results[11].net_rates_of_progress[176,:],'r')
plt.xlabel('Temperature (K)')
plt.ylabel(r'Net rates of progress $\omega$')

plt.annotate('Phi=1.0',xy=(1400,1.25),xytext=(1600,1.25),color='r',arrowprops=dict(arrowstyle='->',color='r'))
plt.annotate('Phi=0.6',xy=(1000,1.1),xytext=(700,0.25),color='b',arrowprops=dict(arrowstyle='->',color='b'))
