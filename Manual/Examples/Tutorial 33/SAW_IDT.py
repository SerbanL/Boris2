"""
This script is part of BORIS

@author: Serban Lepadatu, 2023

Generate a 200 nm wavelength SAW using interdigitated electrodes.
"""

from NetSocks import NSClient, customize_plots
import matplotlib.pyplot as plt
import numpy as np

#setup communication with server
ns = NSClient(); ns.configure(True); customize_plots()

########################################
#Parameters

F0 = 0.5e9
v = 4500
lam = 200e-9

length = 2000e-9
width = 50e-9
thickness = 40e-9
cellsize = 2.5e-9

########################################
#Definition

FM = ns.Material('Fe', [length, width, thickness], [cellsize])
FM.modules(['Zeeman', 'demag', 'exchange', 'anicubi', 'melastic'])
FM.mcellsize([cellsize])

FM.param.mdamping = 1e10
FM.param.mMEc = [0, 0]
FM.param.mdamping.setparamvar('abl_tanh', [0.0, 200e-9/length, 0, 0, 0, 0, 1, 1e6, 200])

FM.surfacefix('-z')
FM.surfacestress([0, 0, thickness, lam/4, width, thickness], 'F0*sin(2*PI*v*t/lam), 0, 0')
FM.surfacestress([lam/2, 0, thickness, 3*lam/4, width, thickness], '-F0*sin(2*PI*v*t/lam), 0, 0')

ns.equationconstants('F0', F0)
ns.equationconstants('lam', lam)
ns.equationconstants('v', v)

ns.setode('LLG', 'RK4', 100e-15)
ns.seteldt(100e-15)

########################################
#Display

FM.display('u')
FM.vecrep(1)
ns.displaydetail(2.5e-9)

########################################
#Simulation

ns.cuda(1)
ns.Relax(['time', 5e-9])

########################################
#Plotting

FM.dp_getexactprofile([cellsize/2, width/2, thickness/2], [length - cellsize/2, width/2, thickness/2], cellsize, 0)
ns.dp_save('SAW_IDT.txt', [0, 1, 2, 3])

data = ns.Get_Data_Columns('SAW_IDT.txt', [0, 1, 2, 3])

plt.axes(xlabel = 'Position (nm)', ylabel = 'Displacement (pm)')
plt.plot(np.array(data[0])/1e-9, np.array(data[1])/1e-12)
plt.savefig('SAW_IDT.png', dpi = 600)
plt.show()



