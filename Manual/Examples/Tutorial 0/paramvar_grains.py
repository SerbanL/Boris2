"""
This script is part of Boris Computational Spintronics v2.8

@author: Serban Lepadatu, 2020
"""

import os
import sys
from NetSocks import NSClient

#setup communication with server
ns = NSClient('localhost')

########################################

#the working directory : same as this script file
directory = os.path.dirname(sys.argv[0]) + "/"
#restore program to default state
ns.default()
ns.chdir(directory)

########################################

ns.meshrect([800e-9, 800e-9, 5e-9])

#setup a grain structure, marking the grains with values in uniform random number distribution between 0 and 1
ns.setparamvar('permalloy', 'J1', 'vor2D', [0.0, 1.0, 20e-9, 1])
#save it to a file and load it here so we can edit it
ns.saveovf2param(paramname = 'J1', filename = 'grains.ovf')
vec, n, rect = ns.Read_OVF2('grains.ovf')

#snap grain values to 0 or 1 : on average half will be 0, the other half 1
for idx in range(len(vec)):    
    if vec[idx] < 0.5: vec[idx] = 0.0
    else: vec[idx] = 1.0

#write it to a file, then load the new grain values
ns.Write_OVF2('grains_snapped.ovf', vec, n, rect)    
ns.setparamvar('permalloy', 'J1', 'ovf2', 'grains_snapped.ovf')

#setup display so we can see the grain structure
ns.display('ParamVar')
ns.setdisplayedparamsvar('permalloy', 'J1')
ns.displaydetail(5e-9)

    
    
    

    




