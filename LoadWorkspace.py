# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 18:53:08 2017

@author: HP USER
"""
import shelve
import os
#import simulations as sim
filename=os.getcwd()+'\\tmp\\konnov_38torr'+'.out'
my_shelf = shelve.open(filename)
for key in my_shelf:
    
    globals()[key]=my_shelf[key]
my_shelf.close()

#print(T)
# Hiya
#print(val)
# [1, 2, 3]