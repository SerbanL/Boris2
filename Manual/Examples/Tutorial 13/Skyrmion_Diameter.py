"""
This script is part of Boris Computational Spintronics

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient, customize_plots
import matplotlib.pyplot as plt
import numpy as np

#setup communication with server
ns = NSClient(); ns.configure(True)
customize_plots()

########################################

#the simulation file
simfile = 'skyrmion_diameter'

#results output file
outputfile = 'skyrmion_diameter.txt'

#field sweep range
start_field = 2e3
end_field = 20e3
field_step = 1e3

ns.loadsim(simfile)

meshrect = ns.meshrect()
cellsize = ns.cellsize()
length = meshrect[3] - meshrect[0]
width = meshrect[4] - meshrect[1]
thickness = meshrect[5] - meshrect[2]

steps = int((end_field - start_field) / field_step)

for i in range(0, steps + 1):

    field = start_field + i * field_step

    ns.setfield([field, 0, 0])

    ns.Run()

    ns.dp_getexactprofile([0, width/2, thickness/2], [length, width/2, thickness/2], cellsize[0], 0)

    skyrmion = ns.dp_fitskyrmion(0, 3)
    diameter = skyrmion[0] * 2
    diameter_std = skyrmion[4] * 2

    ns.SaveDataToFile(outputfile, [field, diameter, diameter_std])

data = ns.Get_Data_Columns(outputfile, [0, 1])
plt.axes(xlabel = 'H (kA/m)', ylabel = 'Skyrmion Diameter (nm)')
plt.plot(np.array(data[0]) / 1e3, np.array(data[1]) / 1e-9, 'o-')
plt.savefig('Figure13.2.png')
plt.show()

    

    

    

    

    

    


    



























