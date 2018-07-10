"""
A burner-stabilized lean premixed hydrogen-oxygen flame at low pressure.
"""

import cantera as ct
import numpy as np
import csv
import os
import pandas as pd
phi = 1.5
burner_d = 30e-3
p = 0.1026*ct.one_atm
p=0.0289474*ct.one_atm
#p = 0.05*ct.one_atm
tburner = 500.0
#mdot = 0.02655058
#reactants = 'H2:'+str(phi)+', O2:0.5, N2:1.88'  # premixed gas composition
#reactants = 'O2:0.12883, N2:0.48466, H2:0.38650'  # premixed gas composition, no quenching
# Next line is reactants but enforcing radical quenching
reactants = 'O2:0.12883, N2:0.48466, H2:0.38650, H:0, O:0, OH:0,HO2:0'
reactants = 'O2:0.12883, N2:0.48466, H2:0.38650'
reactants = {'CH3OCHO':25.2,'O2':50.2,'AR':24.6}
sensSpecies = 'CO'

initial_grid = np.linspace(0.0, 0.03, 10) # m
#initial_grid = np.linspace(0.0, 0.01805882, 10) # m
tol_ss = [1.0e-5, 1.0e-13]  # [rtol atol] for steady-state problem
tol_ts = [1.0e-4, 1.0e-10]  # [rtol atol] for time stepping
loglevel = 1  # amount of diagnostic output (0 to 5)
filename = 'C:\\Users\\HP USER\\Google Drive\\Burke Group\\Codes\\TMRP codes\\dooley\\chem'
gas = ct.Solution(filename+'.cti')
gas.TPX = tburner, p, reactants
v=0.909
gas2=ct.Solution(filename+'.cti')
gas2.TPX=300,p,reactants
mdot=v*gas2.density
#mdot=0.021
f = ct.BurnerFlame(gas, initial_grid)
f.burner.mdot = mdot


f.set_initial_guess()
# read temperature vs. position data from a file.
# The file is assumed to have one z, T pair per line, separated by a comma.

#zloc, tvalues = np.genfromtxt('C:\\Users\\HP USER\\Google Drive\\Burke group\\codes\\paperFigs\\mf1p0_T_profile.csv', delimiter=',', comments='#').T
#zloc, tvalues = np.genfromtxt('C:\\Users\\HP USER\\Desktop\\kevin\\t_profile_38torr.txt', delimiter=',', comments='#').T
data=pd.read_csv('C:\\Users\\HP USER\\Google Drive\\Burke group\\codes\\paperFigs\\mf1p0_T_profile.csv',header=3)
zloc=data['z']
tvalues=data['T']
zloc /= max(zloc)

#print(tvalues)

f.flame.set_fixed_temp_profile(zloc, tvalues)  #sets a fixed temperature profile for the flame simulation.  May come from a measurement. Requires no energy conservation

f.flame.set_steady_tolerances(default=tol_ss)  #Set steady tolerances
f.flame.set_transient_tolerances(default=tol_ts) #Set transient tolerances
f.show_solution()

f.energy_enabled = False  #This must be set to false for a burner stabilized flame with known T-profile

import matplotlib.pyplot as plt
#plt.plot(zloc,tvalues)

f.transport_model = 'Multi'   #Sets to multicomponent transport for simulation.  Needs to be set this way to use Soret effect
f.set_max_jac_age(10, 10)       #Age limits on Jacobian-leave as is for best results
f.solve(loglevel, refine_grid=False)  #Solve for initial estimate without grid refinement
f.save('h2_burner_flame1.xml', 'no_energy',
       'solution with the energy equation disabled')
f.soret_enabled = True          #Enable Soret effect.  Remember transport must be set to multi.  Mix causes failure

f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.5)  #Establishes refinement criteria for grid

#print('mixture-averaged flamespeed = ', f.u[0])

f.transport_model = 'Multi'         #This block solves problem again with grid refinement on
f.solve(loglevel, refine_grid=True)
f.show_solution()
print('multicomponent flamespeed = ', f.u[1])
f.save('h2_burner_flame.xml','energy_multi',
       'solution with the energy equation enabled and multicomponent transport')

#f.write_csv(filename+'_'+str(p)+'_burner_flame.csv', quiet=False)
f.write_csv(filename+'_'+str(p)+'quench'+'_burner_flame.csv', quiet=False)

##Begin section to calculate sensitivities
dk = 0.01
NO_init0 = f.X[f.flame.component_index('CO')-5,len(f.grid)-1]

fooname=filename+'_'+str(p)+'_'+sensSpecies+'.txt'
#sensitivities = np.zeros((gas.n_reactions,len(f.grid)))
with open(fooname,'w') as foo:
        for m in range(gas.n_reactions):
            print m
            gas.set_multiplier(1.0)
            gas.set_multiplier(1 + dk, m)
            f.solve(loglevel = 1, refine_grid = False)
            NO_init = f.X[f.flame.component_index('CO')-4,len(f.grid)-1]
            foo.write('{},{: 10.3e},{} \n'.format(m, (NO_init-NO_init0)/(NO_init0*dk), gas.reaction_equation(m)))
            #np.savetxt(filename = 'phi_' + str(i) + '.txt', output, delimiter = ',', header = 'rxn#, sensetivity(du/dk), reaction')
            #sensitivities[m,:] = np.divide(NO_init-NO_init0,dk*NO_init0)
            
        ###

#np.savetxt("sensitivitiesKLIPP.csv", sensitivities, delimiter="   ")
sens=list(csv.reader(open(fooname,'rb'),delimiter=','))
sens.sort(key=lambda x: -abs(float(x[1])))