"""
This script is part of Boris Computational Spintronics v2.8

@author: Serban Lepadatu, 2020
"""

import os
import sys
from NetSocks import NSClient
import numpy as np
import matplotlib.pyplot as plt

#setup communication with server. By default sent messages are not displayed in console. 
#To enable verbose mode use e.g.: NSClient('localhost', True)
ns = NSClient('localhost')

########################################

#the working directory : same as this script file, typically expecting simulation file to be in same directory as this script file
directory = os.path.dirname(sys.argv[0]) + "/"
#make sure Boris is reset to default in case a ready set simulation file is not loaded below (if so this can be commented out)
ns.default()
#set working directory same as this script file
ns.chdir(directory)

########################################

ns.setode('sLLB', 'TEuler')

#Co layer
ns.setmaterial('Co/Pt', [512e-9, 512e-9, 2e-9])
ns.cellsize([1e-9, 1e-9, 2e-9])
ns.scellsize([4e-9, 4e-9, 2e-9])
ns.setdtstoch(20e-15)
ns.addmodule('Co/Pt', 'iDMexchange')
ns.addmodule('Co/Pt', 'aniuni')
ns.addmodule('Co/Pt', 'heat')
ns.tcellsize([2e-9, 2e-9, 2e-9])
ns.tmodel(2, 'Co/Pt')
ns.curietemperature(500)
ns.pbc('Co/Pt', 'x', 10)
ns.pbc('Co/Pt', 'y', 10)

#Pt layer
ns.addmaterial('Pt', [0.0, 0.0, -8e-9, 512e-9, 512e-9, 0.0])
ns.meshfocus('Pt')
ns.addmodule('Pt', 'heat')
ns.tcellsize([4e-9, 4e-9, 4e-9])
ns.tmodel(2, 'Pt')

#SiO2 layer
ns.addmaterial('SiO2', [0.0, 0.0, -48e-9, 512e-9, 512e-9, -8e-9])
ns.meshfocus('SiO2')
ns.addmodule('SiO2', 'heat')
ns.tcellsize([8e-9, 8e-9, 8e-9])

#general settings
ns.setangle(0, 0)
ns.setfield(100e3, 0, 0)
ns.temperature(300)

#simulation stage
ns.setstage('Qequation', 'Co/Pt')
ns.editstagevalue(0, 'Q0 * exp(-sqrt((x/Lx - 0.5)^2 + (y/Ly - 0.5)^2) / ((d0/Lx)^2/(4*ln(2)))) * exp(-(t-2*tau)^2/(tau^2/(4*ln(2))))')
ns.equationconstants('d0', 400e-9)
ns.equationconstants('tau', 100e-15)
ns.equationconstants('Q0', 4e21)

#output data
ns.setdata('time')
ns.adddata('Q_topo', 'Co/Pt')
ns.adddata('<T>', 'Co/Pt', [255e-9, 255e-9, 0.0, 256e-9, 256e-9, 2e-9])
ns.editdatasave(0, 'time', 10e-15)
ns.savedatafile('ufsky.txt')

#set-up display
ns.meshfocus('Co/Pt')
ns.display('M', 'Co/Pt')
ns.vecrep('Co/Pt', 3)

#definitely need cuda for this one!
ns.cuda(1)

#next time you can just load this and run it
ns.savesim('ufsky_fm')

#0 to 10ps
ns.editstagestop(0, 'time', 10e-12)
ns.setdt(1e-15)
ns.setheatdt(1e-15)
ns.Run()

#10ps to 20ps
ns.editstagestop(0, 'time', 20e-12)
ns.setdt(2e-15)
ns.setheatdt(2e-15)
ns.Run()

#20ps to 60ps
ns.editstagestop(0, 'time', 60e-12)
#mid-simulation re-meshing! 1nm cellsize only needed around the Curie temperature.
#finely tuned simulations of this type means you can run thousands of these events in a reasonable time-scale and still keep accuracy
ns.cellsize([2e-9, 2e-9, 2e-9])
ns.setdt(10e-15)
ns.setheatdt(2e-15)
ns.Run()

#60ps to 800ps
ns.editstagestop(0, 'time', 800e-12)
ns.tcellsize([4e-9, 4e-9, 1e-9])
ns.setdt(50e-15)
#could do with Crank-Nicolson method in next version
ns.setheatdt(5e-15)
ns.editdatasave(0, 'time', 250e-15)
ns.Run()

ns.savemeshimage('ufsky_final')

#now plot |Q| as a function of time
data = ns.Get_Data_Columns('ufsky.txt', [0, 1])
time_ps = [t/1e-12 for t in data[0]]
Qmod = [np.abs(Qval) for Qval in data[1]]

plt.axes(xlabel = 'Time (ps)', ylabel = '|Q|')
plt.xscale('log')
plt.yscale('log')
plt.plot(time_ps, Qmod)
plt.xlim(0.1)
plt.ylim(1e-3)
plt.savefig('ufsky_plot.png', dpi = 600)
plt.show()
