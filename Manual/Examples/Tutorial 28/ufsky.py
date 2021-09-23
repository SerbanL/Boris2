"""
This script is part of Boris Computational Spintronics v3.0

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient, customize_plots
import matplotlib.pyplot as plt
import numpy as np

#setup communication with server
ns = NSClient(); ns.configure(True); customize_plots()

########################################

#for a script running this program entirely from scratch (no loadsim) see Tutorial 0
ns.loadsim('ufsky_fm')
ns.displaydetail(1e-9)

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
