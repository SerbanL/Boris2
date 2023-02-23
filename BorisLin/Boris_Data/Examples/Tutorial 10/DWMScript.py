"""
This script is part of BORIS
@author: Dr. Serban Lepadatu, 2023

Computing current-driven DW velocity.
"""

from NetSocks import NSClient
import matplotlib.pyplot as plt
import numpy as np

ns = NSClient(); ns.configure(True)

########################################

#This is based on Exercise 5.1, done entirely using a Python script

Py = ns.Ferromagnet([320e-9, 80e-9, 20e-9], [5e-9, 5e-9, 5e-9])
Py.modules(['demag', 'exchange', 'Zeeman', 'transport'])
Py.ecellsize([5e-9, 5e-9, 5e-9])
ns.setdefaultelectrodes()

#set spin polarization and STT non-adiabaticity
Py.param.P = 0.4
Py.param.beta = 0.04

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
ns.setode('LLG-STT', 'RK4')
ns.setdt(500e-15)

#save time (s) and dw shift (m) data
ns.setsavedata('cidwmovement_temp.txt', ['stime'], ['dwshift'])

#see the DW shift in data box
Py.addpinneddata('dwshift')

Vrange = np.arange(-4.57e-3, -46e-3, -4.113e-3)
Jcrange = np.array([])
dwvelocity = np.array([])

for V in Vrange:
    
    #first stage achieves steady state movement
    ns.V([V, 'time', 5e-9])
    #second stage captures data
    ns.V([V, 'time', 10e-9, 'time', 1e-12])
    ns.reset()
    
    #process data to extract domain wall velocity
    ns.dp_load('cidwmovement_temp.txt', [0, 1, 0, 1])
    ns.dp_replacerepeats(1)
    dwdata = ns.dp_linreg(0, 1)
    dwvelocity = np.append(dwvelocity, dwdata[0])
    Jc = ns.showdata('<Jc>')
    Jcrange = np.append(Jcrange, Jc[0])
    
    print('Jc (A/m3) = %f, DW velocity (m/s) = %0.4f' % (Jc[0], dwdata[0]))
    
plt.axes(xlabel = 'J$_c$ (TA/m^2)', ylabel = 'DW Velocity (m/s)', title = 'DW Velocity')
plt.plot(np.array(Jcrange)/1e12, dwvelocity, 'o-')
plt.show()