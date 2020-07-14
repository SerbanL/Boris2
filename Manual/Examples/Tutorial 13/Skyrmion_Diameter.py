"""
This script is part of Boris Computational Spintronics v2.8

@author: Serban Lepadatu, 2020
"""

import os
import sys
from NetSocks import NSClient
import matplotlib.pyplot as plt

#setup communication with server
ns = NSClient('localhost')

########################################

#the working directory : same as this script file
directory = os.path.dirname(sys.argv[0]) + "/"
ns.default()
ns.chdir(directory)

########################################

#the simulation file
simfile = 'skyrmion_diameter'

#results output file
outputfile = 'skyrmion_diameter.txt'

#field sweep range
start_field = 2e3
end_field = 20e3
field_step = 1e3

ns.loadsim(directory + simfile)

meshrect = ns.meshrect()
length = meshrect[3] - meshrect[0]
width = meshrect[4] - meshrect[1]
thickness = meshrect[5] - meshrect[2]

steps = int((end_field - start_field) / field_step)

for i in range(0, steps + 1):

    field = start_field + i * field_step

    ns.setfield([field, 0, 0])

    ns.Run()

    ns.dp_getprofile([0, width/2, thickness/2], [length, width/2, thickness/2], 0)

    skyrmion = ns.dp_fitskyrmion(0, 3)
    diameter = skyrmion[0] * 2
    diameter_std = skyrmion[4] * 2

    ns.SaveDataToFile(outputfile, [field, diameter, diameter_std])

data = ns.Get_Data_Columns(outputfile, [0, 1])
plt.axes(xlabel = 'H (A/m)', ylabel = 'diameter (m)', title = 'skyrmion diameter vs H')
plt.plot(data[0], data[1], 'o-')
plt.show()

    

    

    

    

    

    


    



























