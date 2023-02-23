"""
This script is part of Boris Computational Spintronics

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient, customize_plots
import matplotlib.pyplot as plt
import numpy as np

ns = NSClient(); customize_plots()

########################################

data = ns.Get_Data_Columns('fmr_she_fieldsweepFMRshe_data.txt', [0, 1])

H_st = np.array(data[0])/1e3
FMR_st = np.array(data[1])**2

data = ns.Get_Data_Columns('fmr_sot_fieldsweepFMR_data.txt', [0, 1])

H_sot = np.array(data[0])/1e3
FMR_sot = np.array(data[1])**2

plt.axes(xlabel = 'H (kA/m)', ylabel = 'FMR (a.u.)')
plt.plot(H_st, FMR_st, 'o', label = 'J$_c$ = 1e12 A/m$^2$ ST Solver')
plt.plot(H_sot, FMR_sot, '--', label = 'J$_c$ = 1e12 A/m$^2$ SOT')
plt.legend()
plt.savefig('Tutorial21.png')
plt.show()


