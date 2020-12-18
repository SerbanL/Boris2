"""
This script is part of Boris Computational Spintronics v3.0

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient
import matplotlib.pyplot as plt
import numpy as np

ns = NSClient()
ns.configure(True)

########################################

#This is based on Exercise 5.1, done entirely using a Python script

ns.meshrect([320e-9, 80e-9, 20e-9])
ns.cellsize([5e-9, 5e-9, 5e-9])

#setup the moving mesh algorithm for a transverse domain wall along the x axis: 
#1. remove end magnetic charges using two dipole meshes (enables strayfield module)
#2. freeze x-axis ends spins
#3. set a transverse domain wall
#4. enable moving mesh algorithm which keeps the dw centered
ns.preparemovingmesh()

#relax dw in zero field to |mxh| < 10^-5
ns.editstagestop(0, 'mxh', 1e-5)
ns.Run()

#setup 2 field stages, each 5 ns long, but only second one saves data at 1 ps time intervals
ns.setstage('Hxyz')
ns.editstagestop(0, 'time', 5e-9)
ns.addstage('Hxyz')
ns.editstagestop(1, 'time', 5e-9)
ns.editdatasave(1, 'time', 1e-12)

#save time (s) and dw shift (m) data
ns.setdata('stime')
ns.adddata('dwshift')
ns.savedatafile('dwmovement_temp.txt')

#set fixed time-step RK4 method with 500fs time step
ns.setode('LLG', 'RK4')
ns.setdt(500e-15)

#save setup simulation file (next time you can just load it using ns.loadsim('dwmovement'))
ns.savesim('dwmovement')

Hrange = np.arange(100, 2400, 200)
dwvelocity = np.array([])

for H in Hrange:

    ns.reset()
    #first stage achieves steady state movement
    ns.editstagevalue(0, H)
    #second stage captures data
    ns.editstagevalue(1, H)
    ns.Run()
    
    #process data to extract domain wall velocity
    ns.dp_load('dwmovement_temp.txt', [0, 1, 0, 1])
    ns.dp_replacerepeats(1)
    dwdata = ns.dp_linreg(0, 1)
    dwvelocity = np.append(dwvelocity, dwdata[0])
    
    print('H (A/m) = %f, DW velocity (m/s) = %0.4f' % (H, dwdata[0]))
    
plt.axes(xlabel = 'H (A/m)', ylabel = 'DW Velocity (m/s)', title = 'DW Velocity and Walker Breakdown')
plt.plot(Hrange, dwvelocity, 'o-')
plt.show()