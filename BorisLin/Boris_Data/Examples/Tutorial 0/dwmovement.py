"""
This script is part of BORIS
@author: Dr. Serban Lepadatu, 2023

Computing field-driven DW velocity and Walker breakdown.
"""

from NetSocks import NSClient
import matplotlib.pyplot as plt
import numpy as np

ns = NSClient(); ns.configure(True)

########################################

#This is based on Exercise 5.1, done entirely using a Python script

Py = ns.Ferromagnet([320e-9, 80e-9, 20e-9], [5e-9, 5e-9, 5e-9])
Py.modules(['demag', 'exchange', 'Zeeman'])

#setup the moving mesh algorithm for a transverse domain wall along the x axis: 
#1. remove end magnetic charges using two dipole meshes (enables strayfield module)
#2. freeze x-axis ends spins
#3. set a transverse domain wall
#4. enable moving mesh algorithm which keeps the dw centered
ns.preparemovingmesh()

#relax dw in zero field to |mxh| < 10^-5
ns.setode('LLGStatic', 'SDesc')
ns.Relax(['mxh', 1e-5])
ns.reset()

#set fixed time-step RK4 method with 500fs time step
ns.setode('LLG', 'RK4')
ns.setdt(500e-15)

#save time (s) and dw shift (m) data
ns.setsavedata('dwmovement_temp.txt', ['stime'], ['dwshift'])

Hrange = np.arange(100, 2400, 200)
dwvelocity = np.array([])

for H in Hrange:
    
    #first stage achieves steady state movement
    ns.Hxyz([Py, [H, 0, 0], 'time', 5e-9])
    #second stage captures data
    ns.Hxyz([Py, [H, 0, 0], 'time', 10e-9, 'time', 1e-12])
    ns.reset()
    
    #process data to extract domain wall velocity
    ns.dp_load('dwmovement_temp.txt', [0, 1, 0, 1])
    ns.dp_replacerepeats(1)
    dwdata = ns.dp_linreg(0, 1)
    dwvelocity = np.append(dwvelocity, dwdata[0])
    
    print('H (A/m) = %f, DW velocity (m/s) = %0.4f' % (H, dwdata[0]))

plt.axes(xlabel = 'H (A/m)', ylabel = 'DW Velocity (m/s)', title = 'DW Velocity and Walker Breakdown')
plt.plot(Hrange, dwvelocity, 'o-')
plt.show()
