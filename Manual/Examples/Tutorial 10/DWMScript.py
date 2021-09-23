"""
This script is part of Boris Computational Spintronics

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient
import matplotlib.pyplot as plt
import numpy as np

#setup communication with server
ns = NSClient(); ns.configure(True)

########################################

#This is based on Exercise 5.1, done entirely using a Python script

ns.meshrect([320e-9, 80e-9, 10e-9])
ns.cellsize([5e-9, 5e-9, 5e-9])

#set spin polarization and STT non-adiabaticity
ns.setparam('permalloy', 'P', 0.4)
ns.setparam('permalloy', 'beta', 0.04)

#setup the moving mesh algorithm for a transverse domain wall along the x axis: 
#1. remove end magnetic charges using two dipole meshes (enables strayfield module)
#2. freeze x-axis ends spins
#3. set a transverse domain wall
#4. enable moving mesh algorithm which keeps the dw centered
ns.preparemovingmesh()

#relax dw in zero field to |mxh| < 10^-5
ns.setode('LLGStatic', 'SDesc')
ns.editstagestop(0, 'mxh', 1e-5)
ns.Run()

ns.setstage('V')
ns.editstagestop(0, 'time', 5e-9)
ns.addstage('V')
ns.editstagestop(1, 'time', 5e-9)
ns.editdatasave(1, 'time', 10e-12)

#save time (s) and dw shift (m) data
ns.setdata('stime')
ns.adddata('dwshift')
ns.savedatafile('cidwmovement_temp.txt')

#set fixed time-step RK4 method with 500fs time step, enabling STT (LLG-STT equation)
ns.setode('LLG-STT', 'RK4')
ns.setdt(500e-15)

ns.addmodule('permalloy', 'transport')
ns.setdefaultelectrodes()

#save setup simulation file (next time you can just load it using ns.loadsim('dwmovement'))
ns.savesim('cidwm')

Vrange = np.arange(-4.57e-3, -46e-3, -4.113e-3)
Jcrange = np.array([])
dwvelocity = np.array([])

for V in Vrange:

    ns.reset()
    #first stage achieves steady state movement
    ns.editstagevalue(0, V)
    #second stage captures data
    ns.editstagevalue(1, V)
    ns.Run()
    
    #process data to extract domain wall velocity
    ns.dp_load('cidwmovement_temp.txt', [0, 1, 0, 1])
    ns.dp_replacerepeats(1)
    dwdata = ns.dp_linreg(0, 1)
    dwvelocity = np.append(dwvelocity, dwdata[0])
    Jc = ns.showdata('<Jc>')
    Jcrange = np.append(Jcrange, Jc[0])
    
    print('Jc (A/m3) = %f, DW velocity (m/s) = %0.4f' % (Jc[0], dwdata[0]))
    
plt.axes(xlabel = 'Jc (A/m^2)', ylabel = 'DW Velocity (m/s)', title = 'DW Velocity')
plt.plot(Jcrange, dwvelocity, 'o-')
plt.show()