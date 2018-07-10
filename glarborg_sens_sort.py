# -*- coding: utf-8 -*-
"""
Created on Sat Jun 23 16:55:25 2018

@author: HP USER
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

y=results
mechs=['konnov','lam','glarborg_garrison','glarborg_manual']
for i in np.arange(len(y)):
    for o in np.arange(len(y[i].k_sens[0][0])):
        dataframe=pd.DataFrame(pd.DataFrame(columns=['rxn','sens','abs']))
        dataframe['rxn']=y[i].Index[1]
        dataframe['sens']=y[i].k_sens[0,:,o]
        dataframe['abs']=np.abs(y[i].k_sens[0,:,o])
        newdata=dataframe.nlargest(10,'abs')
        newdata=newdata.sort_values(by='abs',ascending=False)
        
        del(newdata['abs'])
        temp=newdata.set_index('rxn')
        print(temp)
        temp.plot.barh(title="Sensitivities for "+mechs[i]+'_'+y[i].Index[2][o],legend=None)
        #threshold = 0.001
        #firstColumn = termolecularSubset.columns[0]
        
        # For plotting, collect only those steps that are above the threshold
        # Otherwise, the y-axis gets crowded and illegible
        #sensitivitiesSubset =termolecularSubset[termolecularSubset[firstColumn].abs() > threshold]
        #print(sensitivitiesSubset)
        #indicesMeetingThreshold =sensitivitiesSubset[firstColumn].abs().sort_values(ascending=False).index
        #sensitivitiesSubset.loc[indicesMeetingThreshold].plot.barh(title="Sensitivities for "+i.mechanism.split('\\')[-1].split('.')[0]+'_'+i.fuel+" at phi="+str(i.phi),
        #                                                      legend=None)
        plt.gca().invert_yaxis()
    
        plt.rcParams.update({'axes.labelsize': 20})
        plt.xlabel(r'Sensitivity: $\frac{\partial\:\ln{S_{u}}}{\partial\:\ln{k}}$');
        plt.savefig('figures\\burner_flame\\'+mechs[i]+'_'+y[i].Index[2][o]+'_38torr.pdf',dpi=1200,bbox_inches='tight')
        