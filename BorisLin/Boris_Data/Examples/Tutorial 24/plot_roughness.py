"""
This script is part of BORIS

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient, customize_plots
import matplotlib.pyplot as plt
import numpy as np

ns = NSClient(); customize_plots()

########################################

data = ns.Get_Data_Columns('roughness_hysteresis_easy.txt', [0, 3])
Hx_e = np.array(data[0])
Mx_e = np.array(data[1])

data = ns.Get_Data_Columns('roughness_hysteresis_hard.txt', [1, 4])
Hx_h = np.array(data[0])
Mx_h = np.array(data[1])

plt.axes(xlabel = 'Field (kA/m)', ylabel = 'Magnetisation (kA/m)')
plt.plot(Hx_e / 1e3, Mx_e / 1e3, '-', label = 'Easy (x)')
plt.plot(Hx_h / 1e3, Mx_h / 1e3, '--', label = 'Hard (y)')
plt.legend()
plt.savefig('hyster_rough.png')
plt.show()


