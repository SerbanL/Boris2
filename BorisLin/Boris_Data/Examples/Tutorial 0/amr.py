"""
This script is part of BORIS
@author: Dr. Serban Lepadatu, 2023

Example of a simple longitudinal AMR loop simulation
"""

from NetSocks import NSClient
import matplotlib.pyplot as plt
import numpy as np

ns = NSClient(); ns.configure(True)

########################################
#This is based on Exercise 8.3, done entirely using a Python script

Py = ns.Ferromagnet([160e-9, 80e-9, 10e-9], [5e-9, 5e-9, 5e-9])
Py.modules(['Zeeman', 'demag', 'exchange', 'transport'])

#set amr percentage of 2%
Py.param.amr = 2.0

#amr loop angle (deg.)
direction_deg = 5
ns.setangle(90, 180.0 + direction_deg)

#set electrodes at x-axis ends with a 1 mV potential drop
ns.setdefaultelectrodes()
ns.setpotential(1e-3)

#save applied field (A/m) and resistance (Ohms)
ns.setsavedata('amr_rawdata.txt', ['Ha', Py], ['R'])

#Run hysteresis loop
ns.setode('LLGStatic', 'SDesc')
ns.Hpolar_seq(
    [Py, [-100e3, 90, direction_deg, 100e3, 90, direction_deg, 200],
     'mxh', 1e-6, 'step'])
ns.Hpolar_seq(
    [Py, [100e3, 90, direction_deg, -100e3, 90, direction_deg, 200],
     'mxh', 1e-6, 'step'])

########################################
#Plot loop

#load all columns from file (0, 1, 2, 3) into internal arrays (0, 1, 2, 4)
ns.dp_load('amr_rawdata.txt', [0, 1, 2, 3, 0, 1, 2, 4])
#get field strength along loop direction and save it in internal array 3
ns.dp_dotprod(0, np.cos(np.radians(direction_deg)), np.sin(np.radians(direction_deg)), 0, 3)
#save field strength and resistance in processed file, then plot it here
ns.dp_save('amr_loop.txt', [3, 4])
amr_data = ns.Get_Data_Columns('amr_loop.txt', [0, 1])
plt.axes(xlabel = 'H (kA/m)', ylabel = 'R (Ohms)', title = 'AMR Loop')
plt.plot(np.array(amr_data[0])/1e3, amr_data[1])
plt.show()
