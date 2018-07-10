# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 18:52:52 2017

@author: HP USER
"""

import shelve
import os
#T='Hiya'
#val=[1,2,3]

filename=os.getcwd()+'\\tmp\\burner_glarborg_78.out'
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
#To restore:

