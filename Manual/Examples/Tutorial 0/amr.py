"""
This script is part of Boris Computational Spintronics

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient
import matplotlib.pyplot as plt
import numpy as np

ns = NSClient(); ns.configure(True)

########################################

#This is based on Exercise 8.3, done entirely using a Python script

ns.meshrect([160e-9, 80e-9, 10e-9])
ns.cellsize([5e-9, 5e-9, 5e-9])

#amr loop angle (deg.)
direction_deg = 5

ns.addmodule('permalloy', 'transport')
#set electrodes at x-axis ends with a 1 mV potential drop
ns.setdefaultelectrodes()
ns.setpotential(1e-3)

ns.setstage('Hpolar_seq')
ns.editstagevalue(0, [-100e3, 90, direction_deg, 100e3, 90, direction_deg, 200])
ns.editstagestop(0, 'mxh', 1e-6)
ns.editdatasave(0, 'step')
ns.addstage('Hpolar_seq')
ns.editstagevalue(1, [100e3, 90, direction_deg, -100e3, 90, direction_deg, 200])
ns.editstagestop(1, 'mxh', 1e-6)
ns.editdatasave(1, 'step')

ns.setangle(90, 180.0 + direction_deg)

#set amr percentage of 2%
ns.setparam('permalloy', 'amr', 2.0)

#save applied field (A/m) and resistance (Ohms)
ns.setdata('Ha')
ns.adddata('R')
ns.savedatafile('amr_rawdata.txt')

ns.setode('LLGStatic', 'SDesc')

ns.Run()

#load all columns from file (0, 1, 2, 3) into internal arrays (0, 1, 2, 4)
ns.dp_load('amr_rawdata.txt', [0, 1, 2, 3, 0, 1, 2, 4])
#get field strength along loop direction and save it in internal array 3
ns.dp_dotprod(0, np.cos(np.radians(direction_deg)), np.sin(np.radians(direction_deg)), 0, 3)
#save field strength and resistance in processed file, then plot it here
ns.dp_save('amr_loop.txt', [3, 4])
amr_data = ns.Get_Data_Columns('amr_loop.txt', [0, 1])
plt.axes(xlabel = 'H (A/m)', ylabel = 'R (Ohms)', title = 'AMR Loop')
plt.plot(amr_data[0], amr_data[1])
plt.show()
