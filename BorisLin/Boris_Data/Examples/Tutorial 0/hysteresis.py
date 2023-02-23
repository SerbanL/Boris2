"""
This script is part of BORIS
@author: Dr. Serban Lepadatu, 2023

Example of a simple hysteresis loop.
"""

import matplotlib.pyplot as plt
import numpy as np
from NetSocks import NSClient
ns = NSClient(); ns.configure(True)

########################################

#Make a ferromagnetic mesh with given rectangle, cubic cellsize, and set modules
Py = ns.Ferromagnet([160e-9, 80e-9, 10e-9], [5e-9])
Py.modules(['demag', 'exchange', 'Zeeman'])

#configure output data : applied field and average magnetisation
ns.setsavedata('hysteresis.txt', ['Ha', Py], ['<M>', Py])

########################################

#run two stages to sweep field up and down between -100 kA/m and +100kA/m in 100 steps, slightly off-axis
#each field step relaxed to mxh < 1e-5, and data saved every field step
ns.setode('LLGStatic', 'SDesc')
ns.Hpolar_seq([Py, [-100e3, 90, 1, +100e3, 90, 1, 100], 'mxh', 1e-5, 'step'])
ns.Hpolar_seq([Py, [100e3, 90, 1, -100e3, 90, 1, 100], 'mxh', 1e-5, 'step'])

########################################

#output file has field (x, y, z components) in columns 0, 1, 2, and average magnetisation (x, y, z components) in columns 3, 4, 5
hysteresis_data = ns.Get_Data_Columns('hysteresis.txt', [0, 3])
#plot Mx vs Hx
plt.axes(xlabel = 'H (kA/m)', ylabel = '<M> (kA/m)')
plt.plot(np.array(hysteresis_data[0])/1e3, np.array(hysteresis_data[1])/1e3)
plt.show()

