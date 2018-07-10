# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 12:55:48 2017

@author: Mark Barbet
"""

import cantera as ct
import numpy as np
import pandas as pd

class flame_model_data:
    def __init__(self,kinetic_sens=np.array(()),physical_sens=np.array(()),Solution=[],Index=[],pIndex=[]):
        self.k_sens=kinetic_sens
        self.p_sens=physical_sens
        self.solution=Solution
        self.Index=Index
        self.pIndex=pIndex
        self.sensitivities = self.all_sensitivities()
        self.overall_index = self.sens_index()
        
    def species_slice_ksens(self,species_name):
        position = self.Index[2].index(species_name)
        return self.k_sens[:,:,position]
    def reaction_slice_ksens(self,reaction):
        position = self.Index[1].index(reaction)
        return self.k_sens[:,position,:]
    def position_slice_ksens(self,location):
        for j in np.arange(len(self.Index[0])):
            if self.Index[0][j]==location:
                return self.k_sens[j,:,:]
            elif j<(len(self.Index[0])-1):
                if self.Index[0][j]<location and self.Index[0][j+1]>location:
                    result = (self.k_sens[j+1,:,:]-self.k_sens[j,:,:])*np.divide((location-self.Index[0][j]),(self.Index[0][j+1]-self.Index[0][j]))+self.k_sens[j,:,:]
                    print('Returning an interpolated result')
                    return result
            else:
                print('Invalid grid location: check position and submit a location in range of flame solution')
    def species_slice_psens(self,species_name):
        position = self.pIndex[2].index(species_name)
        return self.p_sens[:,:,position]
    def parameter_slice_psens(self,parameter):
        position = self.pIndex[1].index(parameter)
        return self.p_sens[:,position,:]
    def position_slice_psens(self,location):
        for j in np.arange(len(self.pIndex[0])):
            if self.pIndex[0][j]==location:
                return self.p_sens[j,:,:]
            elif j<(len(self.pIndex[0])-1):
                if self.pIndex[0][j]<location and self.pIndex[0][j+1]>location:
                    result = (self.p_sens[j+1,:,:]-self.p_sens[j,:,:])*np.divide((location-self.pIndex[0][j]),(self.pIndex[0][j+1]-self.pIndex[0][j]))+self.p_sens[j,:,:]
                    print('Returning an interpolated result')
                    return result
            else:
                print('Invalid grid location: check position and submit a location in range of flame solution')

    def speciation_profile(self,species_name):
        profile = self.solution.X[self.solution.flame.component_index(species_name)-4,:]
        if bool(self.Index)==False:
            try:
                independentVar = np.asarray(self.pIndex[0])
            except:
                print('No index of independent variable provided, please check simulation results manually')
        elif bool(self.Index):            
            independentVar = np.asarray(self.Index[0])
        return np.vstack((independentVar,profile))
    
    def sensitivity_profile(self,observable,parameter):
        position = self.overall_index[2].index(observable)
        position_param = self.overall_index[1].index(parameter)
        sens = self.sensitivities[:,position_param,position]
        independentVar = np.asarray(self.overall_index[0])
        sens = np.flatten(sens)
        return np.vstack((independentVar,sens))
    
    def all_sensitivities(self):
        
        if self.p_sens.any()==False and self.k_sens.any()!=False:
            return self.k_sens
        if self.p_sens.any()!=False and self.k_sens.any()==False:
            return self.p_sens
        if self.p_sens.any()!=False and self.k_sens.any()!=False:
            return np.hstack((self.p_sens,self.k_sens))
        else:
            return np.array(())
        
    def sens_index(self):
        if self.p_sens.any()==False and self.k_sens.any()!=False:
            return self.Index
        if self.p_sens.any()!=False and self.k_sens.any()==False:
            return self.pIndex
        if self.p_sens.any()!=False and self.k_sens.any()!=False:
            return [self.Index[0],self.pIndex[1]+self.Index[1],self.Index[2]]
        else:
            return self.Index
        
    def species_slice_sens(self,species_name):
        position = self.overall_index[2].index(species_name)
        return self.p_sens[:,:,position]
    def parameter_slice_sens(self,parameter):
        position = self.overall_index[1].index(parameter)
        return self.p_sens[:,position,:]
    def position_slice_sens(self,location):
        for j in np.arange(len(self.overall_index[0])):
            if self.overall_index[0][j]==location:
                return self.p_sens[j,:,:]
            elif j<(len(self.overall_index[0])-1):
                if self.overall_index[0][j]<location and self.overall_index[0][j+1]>location:
                    result = (self.p_sens[j+1,:,:]-self.p_sens[j,:,:])*np.divide((location-self.overall_index[0][j]),(self.overall_index[0][j+1]-self.overall_index[0][j]))+self.p_sens[j,:,:]
                    print('Returning an interpolated result')
                    return result
            else:
                print('Invalid grid location: check position and submit a location in range of flame solution')
    
#This function is for pre-mixed burner flame
def burner_flame(gas,grid,mdot,data=pd.DataFrame(columns=['z','T']),kinetic_sens=0,physical_sens=0,observables=[],physical_params=[],energycon=False,soret=True):
    #when energycon is off treat flame as burner stabilized, with known T-profile
   
    baseConditions=gas.TPX
    tol_ss = [1.0e-5, 1.0e-13]  # [rtol atol] for steady-state problem
    tol_ts = [1.0e-4, 1.0e-10]  # [rtol atol] for time stepping
    loglevel = 1  # amount of diagnostic output (0 to 5)
   
    f = ct.BurnerFlame(gas, grid)
    f.burner.mdot = mdot
    
    f.set_initial_guess()
    # read temperature vs. position data from a file.
    # The file is assumed to have one z, T pair per line, separated by a comma.
    if data.empty==False and energycon==False:
        zloc=data['z']
        tvalues=data['T']
        #zloc, tvalues = np.genfromtxt('C:\\Users\\HP USER\\Desktop\\kevin\\harrington_78torr_Tprofile_chemkin.txt', delimiter=',', comments='#').T
        #zloc, tvalues = np.genfromtxt('C:\\Users\\HP USER\\Desktop\\kevin\\t_profile_38torr.txt', delimiter=',', comments='#').T    
        zloc /= max(zloc)

        #print(tvalues)

        f.flame.set_fixed_temp_profile(zloc, tvalues)  #sets a fixed temperature profile for the flame simulation.  May come from a measurement. Requires no energy conservation
    elif data.empty==False and energycon!='off':
        raise Exception('User has supplied fixed temperature dataset but energy conservation is not off.  Remove dataset or turn energy conservation off')
        
        
    f.flame.set_steady_tolerances(default=tol_ss)  #Set steady tolerances
    f.flame.set_transient_tolerances(default=tol_ts) #Set transient tolerances
    f.show_solution()

    f.energy_enabled = energycon  #This must be set to false for a burner stabilized flame with known T-profile

    
    
    f.transport_model = 'Multi'   #Sets to multicomponent transport for simulation.  Needs to be set this way to use Soret effect
    f.set_max_jac_age(10, 10)       #Age limits on Jacobian-leave as is for best results
    f.solve(loglevel, refine_grid=False)  #Solve for initial estimate without grid refinement
    
    f.soret_enabled = soret          #Enable Soret effect.  Remember transport must be set to multi.  Mix causes failure

    f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.1)  #Establishes refinement criteria for grid

    #print('mixture-averaged flamespeed = ', f.u[0])

    f.transport_model = 'Multi'         #This block solves problem again with grid refinement on
    f.solve(loglevel, refine_grid=True)
    f.show_solution()
    print('multicomponent flamespeed = ', f.u[0])
    

    
    #solution = f
    ##Begin section to calculate sensitivities
    dk = 0.010
    solution = f.X
    if kinetic_sens==1 and bool(observables):
        #Calculate kinetic sensitivities
        sensIndex = [f.grid.tolist(),gas.reaction_equations(),observables]
        
        S = np.zeros((len(f.grid),gas.n_reactions,len(observables)))
        #print(solution.X[solution.flame.component_index(observables[0])-4,len(f.grid)-1])
        #a=solution.X[solution.flame.component_index(observables[0])-4,len(f.grid)-1]
        for m in range(gas.n_reactions):
            gas.set_multiplier(1.0)
            gas.set_multiplier(1+2*dk,m)
            f.solve(loglevel=1,refine_grid=False)
            for i in np.arange(len(observables)):
                for k in np.arange(len(f.grid)):                    
                    S[k,m,i]=np.log10(solution[f.flame.component_index(observables[i])-4,k])-np.log10(f.X[f.flame.component_index(observables[i])-4,k])
                    #print(solution.X[solution.flame.component_index(observables[i])-4,k])
                    #print(f.X[f.flame.component_index(observables[i])-4,k])
                    S[k,m,i]=np.divide(S[k,m,i],np.log(dk))             
                    
                    
    
    if physical_sens==1 and bool(observables):
        #Calculate physical sensitivities
        gas.set_multiplier(1.0)
        
        psensIndex = [f.grid.tolist(),physical_params,observables]
        pS = np.zeros((len(f.grid),len(physical_params),len(observables)))
        for m in range(len(physical_params)):
            gas.TPX=baseConditions
            if physical_params[m]=='T':
                gas.TPX=baseConditions[0]+dk,baseConditions[1],baseConditions[2]
            elif physical_params[m]=='P':
                gas.TPX=baseConditions[0]+dk,baseConditions[1]+dk,baseConditions[2]
            f.solve(loglevel=1,refine_grid=False)
            for i in np.arange(len(observables)):
                for k in np.arange(len(f.grid)):
                    pS[k,m,i] =np.log10(solution[f.flame.component_index(observables[i])-4,k])-np.log10(f.X[f.flame.component_index(observables[i])-4,k])
                    pS[k,m,i] = np.divide(pS[k,m,i],np.log10(dk))
                    
                    
                    
    
    elif kinetic_sens==1 and bool(observables)==False:
        raise Exception('Please supply a list of observables in order to run kinetic sensitivity analysis')
    elif physical_sens==1 and bool(observables)==False:
        raise Exception('Please supply a list of observables in order to run physical sensitivity analysis')
    gas.set_multiplier(1.0)
    gas.TP = baseConditions[0],baseConditions[1]  
    f.solve(loglevel=1,refine_grid=False)
    
    
    
    
    if kinetic_sens==1 and bool(observables) and physical_sens!=1:
        results = flame_model_data(kinetic_sens=S,Solution=f,Index=sensIndex)        
        return results
    elif physical_sens==1 and bool(observables) and kinetic_sens!=1:
        results = flame_model_data(Solution=f,pIndex=psensIndex,physical_sens=pS)
        return results
    elif kinetic_sens==1 and physical_sens==1 and bool(observables):
        results = flame_model_data(Solution=f,pIndex=psensIndex,Index=sensIndex,physical_sens=pS,kinetic_sens=S)
        return results
    elif kinetic_sens!=1 and physical_sens!=1:
        Index = [f.grid.tolist()] 
        results = flame_model_data(Solution=f,Index=Index)
        return results
    else:
        print('Something went wrong with the parameters given.  Kinetic and physical sens may be set to either 0 or 1, and require a list of observables')
            
            
    
    