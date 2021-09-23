"""
This script is part of Boris Computational Spintronics v3.0

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient, customize_plots
import matplotlib.pyplot as plt
import numpy as np

#setup communication with server
ns = NSClient(); ns.configure(True); customize_plots()

########################################

simfile = 'circle_staircase'
outfile = 'circle_staircase_data.txt'

ns.loadsim(simfile)

for angle in range(0, 360):

    ns.setangle(90, angle)
    ns.editstagevalue(0, [1e6 * np.cos(np.radians(angle)), 1e6 * np.sin(np.radians(angle)), 0])
    ns.reset()

    ns.Run()

    e_demag = ns.showdata('e_demag')
    e_rough = ns.showdata('e_rough')

    ns.SaveDataToFile(outfile, [angle, e_demag, e_rough])
    
data = ns.Get_Data_Columns(outfile, [0, 1, 2])
demag_nocorrrections = [data[1][i] for i in range(len(data[1]))]
demag_withcorrrections = [data[1][i] + data[2][i] for i in range(len(data[1]))]

plt.axes(xlabel = 'Angle (deg.)', ylabel = 'Demag energy demsity (J/m^3)')
plt.plot(data[0], demag_nocorrrections, label = 'without corrections')

plt.plot(data[0], demag_withcorrrections, label = 'with corrections')
plt.legend()
plt.savefig('demag_correction.png')
plt.show()

    

    
    


