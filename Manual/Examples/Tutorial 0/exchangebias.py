"""
This script is part of Boris Computational Spintronics v2.8

@author: Serban Lepadatu, 2020
"""

import os
import sys
from NetSocks import NSClient
import matplotlib.pyplot as plt

#setup communication with server
ns = NSClient('localhost')

########################################

#the working directory : same as this script file
directory = os.path.dirname(sys.argv[0]) + "/"
#restore program to default state
ns.default()
ns.chdir(directory)

########################################

ns.setafmesh('Antiferromagnet', [320e-9, 320e-9, 10e-9])
ns.cellsize([5e-9, 5e-9, 5e-9])
ns.addmodule('Antiferromagnet', 'aniuni')
ns.addmodule('Antiferromagnet', 'surfexchange')
#set sub-lattice A magnetisation to result in biasing towards +ve side
ns.setangle(90, 180)
#Add Fe mesh on top of the antiferromagnet
ns.addmaterial('Fe', [0, 0, 10e-9, 320e-9, 320e-9, 12e-9])
ns.meshfocus('Fe')
#need smaller cellsize for Fe (in Boris meshes can be independently discretised)
ns.cellsize([2.5e-9, 2.5e-9, 2e-9])
ns.pbc('Fe', 'x', 10)
ns.pbc('Fe', 'y', 10)
ns.addmodule('Fe', 'anicubi')
#Now both the Antiferromagnet and Fe meshes have surfexchange module enabled, so exchange bias field will be included in computations
ns.addmodule('Fe', 'surfexchange')
#set bilinear surface exchange coupling value - exchange bias is proportional to this
ns.setparam('Fe', 'J1', 0.2e-3)

ns.setode('LLGStatic', 'RKF45')

ns.setstage('Hpolar_seq', 'supermesh')
ns.editstagevalue(0, [-50e3, 90, 5, 100e3, 90, 5, 50])
ns.editstagestop(0, 'mxh', 1e-4)
ns.editdatasave(0, 'step')
ns.addstage('Hpolar_seq', 'supermesh')
ns.editstagevalue(1, [100e3, 90, 5, -50e3, 90, 5, 50])
ns.editstagestop(1, 'mxh', 1e-4)
ns.editdatasave(1, 'step')

ns.setdata('Ha')
ns.adddata('<M>', 'Fe')
ns.savedatafile('exchangebias.txt')

ns.cuda(1)
ns.Run()

#we should really project along the 5 degree direction, but will keep this simple
data = ns.Get_Data_Columns('exchangebias.txt', [0, 3])

plt.axes(xlabel = 'H (A/m)', ylabel = 'M (A/m)', title = 'Exchange bias')
plt.plot(data[0], data[1])
plt.show()


