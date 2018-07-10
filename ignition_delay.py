# -*- coding: utf-8 -*-
"""
Created on Tue Jan 09 11:16:38 2018

@author: Mark Barbet
"""

"""
A simplistic approach to mechanism reduction which demonstrates Cantera's
features for dynamically manipulating chemical mechanisms.

Here, we use the full GRI 3.0 mechanism to simulate adiabatic, constant pressure
ignition of a lean methane/air mixture. We track the maximum reaction rates for
each reaction to determine which reactions are the most important, according to
a simple metric based on the relative net reaction rate.

We then create a sequence of reduced mechanisms including only the top reactions
and the associated species, and run the simulations again with these mechanisms
to see whether the reduced mechanisms with a certain number of species are able
to adequately simulate the ignition delay problem.
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import os
import efficiency_manipulate as em
import copy
filename=os.getcwd()+'\\TMRP codes\\FFCM-1\\FFCM1.cti'
gas = ct.Solution(filename)
initial_state = 1200, 1 * ct.one_atm, 'CH3OH:1.6, O2:1.0'

# Run a simulation with the full mechanism
gas.TPX = initial_state
r = ct.IdealGasConstPressureReactor(gas)
sim = ct.ReactorNet([r])

tt = []
TT = []
t = 0.0
# Rmax is the maximum relative reaction rate at any timestep
Rmax = np.zeros(gas.n_reactions)
while t < 0.02:
    t = sim.step()
    tt.append(1000 * t)
    TT.append(r.T)
    rnet = abs(gas.net_rates_of_progress)
    rnet /= max(rnet)
    Rmax = np.maximum(Rmax, rnet)
    
gas2=em.efficiency_rate_swap(ct.Solution(filename))

import soln2cti as ctiw
new_file=ctiw.write(gas2)
gas2=ct.Solution(os.getcwd()+'\\pym_gas.cti')
#initial_state2=copy.deepcopy(initial_state)
gas2.TPX=initial_state
r2=ct.IdealGasConstPressureReactor(gas2)
sim2=ct.ReactorNet([r2])

tt2=[]
TT2=[]
t2=0.0
Rmax2=np.zeros(gas2.n_reactions)
while t2<0.02:
    t2=sim2.step()
    tt2.append(1000*t2)
    TT2.append(r2.T)
    rnet2=abs(gas2.net_rates_of_progress)
    rnet2 /= max(rnet2)
    Rmax2=np.maximum(Rmax2,rnet2)
    
#calculate derivatives of temperature with time
dTdt=np.gradient(TT,tt)



plt.plot(tt, TT, label='Baseline', color='k', lw=3, zorder=100)
#plt.figure()
plt.plot(tt2,TT2,label='Manipulated', color='r',lw=2,zorder=100)
plt.savefig('ig_delays.pdf',dpi=1200)

#Now find the ignition delays


# Get the reaction objects, and sort them so the most active reactions are first
R = sorted(zip(Rmax, gas.reactions()), key=lambda x: -x[0])

# Test reduced mechanisms with different numbers of reactions
#C = plt.cm.winter(np.linspace(0,1,5))
#for i,N in enumerate([40,50,60,70,80]):
#    # Get the N most active reactions
#    reactions = [r[1] for r in R[:N]]
#
#    # find the species involved in these reactions. At a minimum, include all
#    # species in the reactant mixture
#    species_names = {'N2', 'CH4', 'O2'}
#    for reaction in reactions:
#        species_names.update(reaction.reactants)
#        species_names.update(reaction.products)
#
#    # Get the species objects
#    species = [gas.species(name) for name in species_names]
#
#    # create the new reduced mechanism
#    gas2 = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
#                       species=species, reactions=reactions)
#
#    # Re-run the ignition problem with the reduced mechanism
#    gas2.TPX = initial_state
#    r = ct.IdealGasConstPressureReactor(gas2)
#    sim = ct.ReactorNet([r])
#
#    t = 0.0
#
#    tt = []
#    TT = []
#    while t < 0.02:
#        t = sim.step()
#        tt.append(1000 * t)
#        TT.append(r.T)
#
#    plt.plot(tt,TT, lw=2, color=C[i],
#             label='K={0}, R={1}'.format(gas2.n_species, N))
#    plt.xlabel('Time (ms)')
#    plt.ylabel('Temperature (K)')
#    plt.legend(loc='upper left')
#    plt.title('Reduced mechanism ignition delay times\n'
#              'K: number of species; R: number of reactions')
#    plt.xlim(0, 20)
#    plt.tight_layout()
#
#plt.show()