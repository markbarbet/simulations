# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 18:17:43 2018

@author: HP USER
"""

"""
Constant-volume, adiabatic kinetics simulation.
This is how most people model an idealized shock tube.
"""

import sys
import numpy as np

import cantera as ct

gas = ct.Solution('gri30.xml')

gas.TPX = 1001.0, ct.one_atm, 'H2:2,O2:1,N2:4'
r = ct.IdealGasReactor(gas)

sim = ct.ReactorNet([r])
time = 0.0
times = np.zeros(100)
data = np.zeros((100,4))
for i in range(gas.n_reactions):
            r.add_sensitivity_reaction(i)
print('%10s %10s %10s %14s' % ('t [s]','T [K]','P [Pa]','u [J/kg]'))
for n in range(100):
    time += 1.e-5
    
    sim.advance(time)
    times[n] = time # time in s
    data[n,0] = r.T
    data[n,1:] = r.thermo['OH','H','H2'].X
    print('%10.3e %10.3f %10.3f %14.6e' % (sim.time, r.T,
                                           r.thermo.P, r.thermo.u))
test=sim.sensitivities()
# Plot the results if matplotlib is installed.
# See http://matplotlib.org/ to get it.
if '--plot' in sys.argv[1:]:
    import matplotlib.pyplot as plt
    plt.clf()
    plt.subplot(2, 2, 1)
    plt.plot(times*1e3, data[:,0])
    plt.xlabel('Time (ms)')
    plt.ylabel('Temperature (K)')
    plt.subplot(2, 2, 2)
    plt.plot(times*1e3, data[:,1])
    plt.xlabel('Time (ms)')
    plt.ylabel('OH Mole Fraction')
    plt.subplot(2, 2, 3)
    plt.plot(times*1e3, data[:,2])
    plt.xlabel('Time (ms)')
    plt.ylabel('H Mole Fraction')
    plt.subplot(2, 2, 4)
    plt.plot(times*1e3,data[:,3])
    plt.xlabel('Time (ms)')
    plt.ylabel('H2 Mole Fraction')
    plt.tight_layout()
    plt.show()
else:
    print("To view a plot of these results, run this script with the option --plot")