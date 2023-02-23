"""
This script is part of BORIS

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient, customize_plots
import matplotlib.pyplot as plt
import numpy as np

ns = NSClient(); customize_plots()

########################################

data = ns.Get_Data_Columns('cppgmr_switching.txt', [0, 1, 2])

time = np.array(data[0])
R = np.array(data[1])
Mx = np.array(data[2])

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

ax1.set(xlabel = 'Time (ns)', ylabel = 'R ($\Omega$)')
ax1.plot(time / 1e-9, R, 'o--', label = 'Resistance')
ax2.set(ylabel = 'Magnetisation (kA/m)')
ax2.plot(time / 1e-9, Mx / 1e3, '-', color = 'black', label = 'Magnetisation')
fig.legend(loc='lower right', bbox_to_anchor=(0.93, 0.18))
plt.savefig('gmr.png')
plt.show()


