# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 11:53:46 2017

@author: Mark Barbet
"""

#Free Flame simulator for use in TMRP project

import numpy as np
import cantera as ct
import pandas as pd
import os
import free_flame as ff
import efficiency_manipulate as em
#import simulations as sim
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
import matplotlib.pyplot as plt
if 'results' not in locals().keys():
    results=[]

phi=[0.6,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.2,3.6,4.0,4.4,4.8,5.0]
phi=np.arange(0.4,3.0,0.1)
phi=np.arange(0.4,5.0,0.2)
#phi=np.arange(0.6,2.5,0.1)
#phi=[0.4,0.6,0.8,1.0,1.4,1.8,2.0]
#phi=[1]
fuels=['NC10H22','C3H8','C4H6','C3H6','pC3H4','aC3H4','cC3H4']
fuels=['CH4']
fuels=['H2','C2H6','C2H4',{'H2':1,'CO':1}]
fuels=['C3H8','C4H6','C3H6','pC3H4','aC3H4','cC3H4']
#fuels=['NC10H22']
fuels=['H2',{'H2':1,'CO':1}]
fuels=['H2']
fuels=['NC10H22']
fuels=['C2H2']
fuels=['CH3OH']
#fuels=[{'H2':1,'CO':1}]
#fuels=['CH4']
#fuels=['C3H8']
#fuels=['C2H6']
oxidizers=[{'O2':0.5,'N2':1.88},{'O2':1, 'N2':3.76}]
oxidizers=[{'O2':16.5,'N2':16.5*3.76}]
#oxidizers=[{'O2':15.5, 'N2':58.28}]
oxidizers=[{'O2':1.5, 'N2':1.5*3.76}]
#oxidizers=[{'O2':5.0, 'N2':5*3.76}]
#oxidizers = [{'O2':2.5,'N2':2.5*3.76}]
temps=[float(600)]
#temps=[float(400)]
pressures=[float(1)] #in atmospheres
efficiency_swap=True
#oxidizers=[{'O2':2,'N2':1.5*3.76}]
fuels=['ch3och3']
oxidizers=[{'o2':3.0,'n2':3.76*3.0}]
fuels=['H2','CH3OH','C2H2']
fuels=['CH3OH','C2H2']
fuels=[{'H2':1,'CO':1}]
fuels=['C3H8']
fuels=['H2']
oxidizers=[{'O2':0.5,'N2':1.88},{'O2':2.0,'N2':2.0*3.76},{'O2':2.5,'N2':2.5*3.76}]
oxidizers=[{'O2':2.0,'N2':2.0*3.76},{'O2':2.5,'N2':2.5*3.76}]
oxidizers=[{'O2':1.0,'N2':3.76}]
oxidizers=[{'O2':5.0,'N2':5.0*3.76}]
oxidizers=[{'O2':0.5,'N2':0.5*3.76}]
mechanisms=[os.getcwd()+'\\TMRP codes\\JetSurf\\jetsurf2.cti',
            os.getcwd()+'\\TMRP codes\\JetSurf\\appendedinput_beta001_forward_chem.cti']
mechanisms=[os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_top27_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_top27_FFCM1.cti']
mechanisms=[os.getcwd()+'\\TMRP codes\\FFCM-1\\FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_cutrxns01_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_cutrxns01_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_top29_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_top29_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_top30_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_top30_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_FFCM1.cti']
mechanisms=[os.getcwd()+'\\TMRP codes\\FFCM-1\\FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_modify948_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_modify948_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_top30_FFCM1.cti']
mechanisms=[os.getcwd()+'\\TMRP codes\\FFCM-1\\FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_top30_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_top30_vinoxy_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_vinoxy_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_top30_vinoxy_FFCM1.cti']
mechanisms=[os.getcwd()+'\\TMRP codes\\Metcalfe\\chem.cti',
            os.getcwd()+'\\TMRP codes\\Metcalfe\\appendedinput_beta001_forward_top30_chem.cti']
mechanisms=[os.getcwd()+'\\TMRP codes\\FFCM-1\\FFCM1.cti']
#mechanisms=[os.getcwd()+'\\TMRP codes\\JetSurf\\jetsurf2.cti']
#mechanisms=[os.getcwd()+'\\TMRP codes\\USCMech\\uscmech.cti',
#            os.getcwd()+'\\TMRP codes\\USCMech\\appendedinput_beta001_uscmech.cti',
#            os.getcwd()+'\\TMRP codes\\USCMech\\appendedinput_beta001_forward_uscmech.cti']
#mechanisms=[os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_FFCM1.cti']
#mechanisms=[os.getcwd()+'\\TMRP codes\\FFCM-1\\FFCM1_H.cti',
#            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_FFCM1_H.cti',
#            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_FFCM1_H.cti',
#            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_specified2500.0_FFCM1_H.cti',
#            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_specified2500.0_FFCM1_H.cti']
#mechanisms=[os.getcwd()+'\\TMRP codes\\FFCM-1\\FFCM1_H.cti']
#mechanisms=[os.getcwd()+'\\TMRP codes\\USCMech\\uscmech.cti']
mechanisms=[os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem.cti',
            os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem_top1_manipulated.cti',
            os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem_top2_manipulated.cti',
            os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem_top3_manipulated.cti',
            os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem_top4_manipulated.cti',
            os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem_top5_manipulated.cti',
            os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem_top6_manipulated.cti',
            os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem_top7_manipulated.cti',
            os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem_top8_manipulated.cti',
            os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem_top9_manipulated.cti',
            os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem_top10_manipulated.cti',
            os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem_top11_manipulated.cti',
            os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem_top12_manipulated.cti',
            os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem_top13_manipulated.cti',
            os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem_top14_manipulated.cti']
mechanisms=[os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem.cti']
mechanisms=[os.getcwd()+'\\TMRP codes\\FFCM-1\\FFCM1.cti']
mechanisms=[os.getcwd()+'\\Mechanisms\\CH4_DME\\chem.cti']
mechanisms=[os.getcwd()+'\\Mechanisms\\FFCM-1\\FFCM1.cti',
            os.getcwd()+'\\Mechanisms\\FFCM-1\\FFCM1_natureTermoleculars.cti',
            os.getcwd()+'\\Mechanisms\\FFCM-1\\FFCM1_natureTermolecularsx10.cti']
mechanisms=[os.getcwd()+'\\Mechanisms\\FFCM-1\\FFCM1.cti',
            os.getcwd()+'\\Mechanisms\\FFCM-1\\FFCM1_plusHCO.cti',
            os.getcwd()+'\\Mechanisms\\FFCM-1\\FFCM1_natureTermoleculars.cti']
mechanisms=[os.getcwd()+'\\Mechanisms\\Aramco2.0\\chem.cti',
            os.getcwd()+'\\Mechanisms\\Aramco2.0\\modified_chem.cti']
mechanisms=[os.getcwd()+'\\Mechanisms\\UdatedH2Model110725\\chem.cti',
            os.getcwd()+'\\Mechanisms\\UdatedH2Model110725\\modified_chem.cti']
#mechanisms=[os.getcwd()+'\\Mechanisms\\dryer_methanol\\chem.cti',
#            os.getcwd()+'\\Mechanisms\\dryer_methanol\\chem_natureTermolecularsfo.cti',
#            os.getcwd()+'\\Mechanisms\\dryer_methanol\\chem_natureTermolecularsx10fo.cti']
#mechanisms=[os.getcwd()+'\\Mechanisms\\FFCM-1\\FFCM1.cti',
#            os.getcwd()+'\\Mechanisms\\FFCM-1\\FFCM1_natureTermolecularsfo.cti',
#            os.getcwd()+'\\Mechanisms\\FFCM-1\\FFCM1_natureTermolecularsx10fo.cti']
efficiency_manipulate=[True,False]
efficiency_manipulate=[False]
global conditionsTups
width=0.03
conditionsTups=[]
for i in np.arange(len(phi)):
    for filename in np.arange(len(mechanisms)):
        for j in np.arange(len(fuels)):
                for t in np.arange(len(temps)):
                    for p in np.arange(len(pressures)):
                        for truth in efficiency_manipulate:
                            conditionsTups.append([mechanisms[filename],phi[i],fuels[j],oxidizers[j],temps[t],pressures[p],truth])
if True in efficiency_manipulate:
    gas=ct.Solution(mechanisms[0])
    gas=em.efficiency_rate_swap(gas)
    gas.name=fuels[0]+'_'+mechanisms[0].rstrip('.cti').split('\\')[-1]
    import soln2cti as ctiw
    new_file=ctiw.write(gas)
#print(conditionsTups)
def solver(conditionsTup):
    try:
        if not conditionsTup[6]:
            fuels=conditionsTup[2]
            #oxidizer={'O2':0.5, 'N2':1.88}
            oxidizer=conditionsTup[3]
            gas=ct.Solution(conditionsTup[0])
            
           
            
            pressures=conditionsTup[5]
            #temps=[298]
            temps=conditionsTup[4]
            width=0.03
            results={}
            gas.TP=temps,pressures*ct.one_atm
            print(conditionsTup)
            results=ff.free_flame(conditionsTup[1],fuels,oxidizer,gas,width,kinetic_sens=0,energycon=True,flamespeed_sens=1,soret=False)
            results.add_mechanism(conditionsTup[0])
            results.tags('efficiency_manipulated:False')

        
        elif conditionsTup[6]:
            fuels=conditionsTup[2]
            #oxidizer={'O2':0.5, 'N2':1.88}
            oxidizer=conditionsTup[3]
            #gas=ct.Solution(conditionsTup[0])
            #If you're doing tests with collision efficiencies use the following function, otherwise comment it out
            #gas=em.efficiency_rate_swap(gas)
            #import soln2cti as ctiw
            #new_file=ctiw.write(gas)
            gas=ct.Solution(os.getcwd()+'\\pym_'+fuels+'_'+conditionsTup[0].rstrip('.cti').split('\\')[-1]+'.cti')
            pressures=conditionsTup[5]
            #temps=[298]
            temps=conditionsTup[4]
            width=0.03
            results={}
            gas.TP=temps,pressures*ct.one_atm
            print(conditionsTup)
            results=ff.free_flame(conditionsTup[1],fuels,oxidizer,gas,width,kinetic_sens=0,energycon=True,flamespeed_sens=1,soret=False)
            results.add_mechanism(conditionsTup[0])
            results.tags('efficiency_manipulated:True')
        return results
    except:
        results='failed'
        return results


pool = ThreadPool(4) 
#results = results+pool.map(solver,conditionsTups)
results=pool.map(solver,conditionsTups)
for j in results:
    if j!='failed':
        if 'dict' in str(type(j.fuel)):
            tempstr=''
            for i in j.fuel.keys():
                tempstr=tempstr+i+'+'
            tempstr=tempstr.rstrip('+')
            j.fuel=tempstr
import shelve
import os
#T='Hiya'
#val=[1,2,3]

filename=os.getcwd()+'\\tmp\\flamespeeds_'+results[0].fuel+'_FFCM1_efficiency_swap.out'
my_shelf = shelve.open(filename,'n') # 'n' for new

for key in dir():
    if key=='results':
        print(key)
        my_shelf[key] = globals()[key]
    #except TypeError:
        #
        # __builtins__, my_shelf, and imported modules can not be shelved.
        #
        #print('ERROR shelving: {0}'.format(key))
my_shelf.close()
#runfile('C:/Users/HP USER/Google Drive/Burke Group/Codes/TMRP_analysis.py', wdir='C:/Users/HP USER/Google Drive/Burke Group/Codes')
#a=[]
#resultsSubset=[]
#for i in results:
#    if i!='failed':
#        resultsSubset.append(i)
#for d in resultsSubset:
#    a.append(d)
#phi_list=[]
#for i in np.arange(len(resultsSubset)):
#    phi_list.append(resultsSubset[i].phi)
#b=pd.DataFrame(data=phi_list)
#b.columns=['Phi']
#u=[]
#for i in np.arange(len(resultsSubset)):
#    u.append(resultsSubset[i].solution['u'][0])
#b['u']=u
##b['solution']=results.values()
#
#mech=[]
#
#for j in np.arange(len(resultsSubset)):
#    mech.append(resultsSubset[j].mechanism)
#mechList=[]
#for i in mech:
#    if i not in mechList:
#        mechList.append(i)
#solutions=[]
#for i in np.arange(len(resultsSubset)):
#    solutions.append(resultsSubset[i])
#b['mechanism']=mech
#b['solution']=solutions
#b=b.sort_values(by='Phi', ascending=1)
#for k in mechList:
#    subset=b[b['mechanism']==k]
#    plt.plot(subset['Phi'],subset['u'])
#for i in fuels:
#    b[i]=[False]*len(resultsSubset)
#for i in np.arange(len(resultsSubset)):
#    for k in fuels:
#        if b['solution'][i].fuel==k:
#            b[k][i]=True
#markerlist=['r:','g-','b--','k+','k']    
#for j in fuels:
#    plt.figure()
#    for k in mechList:
#        subset=b[b['mechanism']==k]
#        subset=subset[subset[j]==True]
#        plt.plot(subset['Phi'],subset['u'],markerlist[mechList.index(k)],linewidth=2,label=k.split('\\')[-1])
#        plt.xlabel('Phi')
#        plt.ylabel('Flame Speed (m/s)')
#        
#        plt.title('Flame Speeds for '+j)
#    plt.legend(loc=4,frameon=True,prop={'size':10})
#    plt.savefig(j+'_flamespeed.pdf',dpi=600)
#    
#y=[]
#for g in results:
#    if g.phi==0.6:
#        y.append(g)
#sortbool=[False]*len(y[1].flamespeed_sens)
#for i in y[3:]:
#    for j in np.arange(len(i.flamespeed_sens)):
#        if i.flamespeed_sens.index[j] in y[0].flamespeed_sens.index:
#            sortbool[j]=False
#        else:
#            sortbool[j]=True
#    termolecularSubset=i.flamespeed_sens[sortbool]
#    threshold = 0.01
#    firstColumn = termolecularSubset.columns[0]
#
#    # For plotting, collect only those steps that are above the threshold
#    # Otherwise, the y-axis gets crowded and illegible
#    sensitivitiesSubset =termolecularSubset[termolecularSubset[firstColumn].abs() > threshold]
#    #print(sensitivitiesSubset)
#    indicesMeetingThreshold =sensitivitiesSubset[firstColumn].abs().sort_values(ascending=False).index
#    sensitivitiesSubset.loc[indicesMeetingThreshold].plot.barh(title="Sensitivities for "+i.mechanism.split('\\')[-1].split('.')[0]+'_'+i.fuel+" at phi="+str(i.phi),
#                                                          legend=None)
#    plt.gca().invert_yaxis()
#
#    plt.rcParams.update({'axes.labelsize': 20})
#    plt.xlabel(r'Sensitivity: $\frac{\partial\:\ln{S_{u}}}{\partial\:\ln{k}}$');
#
#    # Uncomment the following to save the plot. A higher than usual resolution (dpi) helps
#    plt.savefig(i.mechanism.split('\\')[-1].split('.')[0]+'_'+i.fuel+'_phi06'+'.pdf'
#    ,dpi=600,bbox_inches="tight")
