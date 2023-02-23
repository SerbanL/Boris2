"""
This script is part of BORIS
@author: Dr. Serban Lepadatu, 2023

Skyrmion creation using ultrafast laser pulse in a Co/Pt bilayer on SiO2 substrate.
"""

from NetSocks import NSClient
import matplotlib.pyplot as plt
import numpy as np

ns = NSClient(); ns.configure(True)

########################################
#Setup trilayer

#Co layer
FM = ns.Material('Co/Pt', [512e-9, 512e-9, 2e-9], [1e-9, 1e-9, 2e-9])
FM.modules(['Zeeman', 'demag', 'iDMexchange', 'aniuni', 'heat'])
FM.scellsize([4e-9, 4e-9, 2e-9])
FM.tcellsize([2e-9, 2e-9, 1e-9])
FM.tmodel(2)
FM.curietemperature(500)
FM.pbc('x', 10)
FM.pbc('y', 10)
FM.setangle(0, 0)
FM.setfield(100e3, 0, 0)

#Pt layer
HM = ns.Material('Pt', [0.0, 0.0, -8e-9, 512e-9, 512e-9, 0.0], [4e-9, 4e-9, 4e-9])
HM.modules(['heat'])
HM.tmodel(2)

#SiO2 layer
Substrate = ns.Material('SiO2', [0.0, 0.0, -48e-9, 512e-9, 512e-9, -8e-9], [8e-9, 8e-9, 8e-9])
Substrate.modules(['heat'])

ns.temperature(300)

########################################
#Solver

#RK4 with quartic polynomial extrapolation of demag field
ns.setode('sLLB', 'RK4')
ns.evalspeedup(5)
ns.setdtstoch(20e-15)

########################################
#Set-up display
FM.display('M')
FM.vecrep(3)
ns.displaydetail(1e-9)

########################################
#Ouput data and heat source

ns.setsavedata('ufsky.txt', ['time'], ['Q_topo', FM], ['<T>', FM, [255e-9, 255e-9, 0.0, 256e-9, 256e-9, 2e-9]])

#heat source due to laser spot in Co layer
FM.param.Q = 2e21
FM.param.Q.setparamvar('equation', 'exp(-sqrt((x/Lx - 0.5)^2 + (y/Ly - 0.5)^2) / ((d0/Lx)^2/(4*ln(2)))) * exp(-(t-2*tau)^2/(tau^2/(4*ln(2))))')
ns.equationconstants('d0', 400e-9)
ns.equationconstants('tau', 100e-15)

########################################
#Simulate

ns.cuda(1)

#0 to 10ps
ns.setdt(1e-15)
ns.setheatdt(1e-15)
ns.Relax(['time', 10e-12, 'time', 10e-15])

#10ps to 20ps
ns.setdt(4e-15)
ns.setheatdt(2e-15)
ns.Relax(['time', 20e-12, 'time', 10e-15])

#20ps to 60ps
#mid-simulation re-meshing! 1nm cellsize only needed around the Curie temperature.
FM.cellsize([2e-9, 2e-9, 2e-9])
ns.setdt(20e-15)
ns.setheatdt(2e-15)
ns.Relax(['time', 60e-12, 'time', 10e-15])

#60ps to 800ps
FM.tcellsize([4e-9, 4e-9, 1e-9])
ns.setdt(100e-15)
ns.setheatdt(5e-15)
ns.Relax(['time', 800e-12, 'time', 250e-15])

########################################
#Plotting

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
