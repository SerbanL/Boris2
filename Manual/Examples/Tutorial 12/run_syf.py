"""
This script is part of Boris Computational Spintronics v2.8

@author: Serban Lepadatu, 2020
"""

import os
import sys
from NetSocks import NSClient
import matplotlib.pyplot as plt
import numpy as np

#setup communication with server
ns = NSClient('localhost')

########################################

#the working directory : same as this script file
directory = os.path.dirname(sys.argv[0]) + "/"
ns.chdir(directory)

########################################

ns.loadsim('syf_bilayer')

ns.Run()

u = [np.cos(np.radians(1)), np.sin(np.radians(1)), 0]
ns.dp_load('syf_hysteresis', [0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 5, 6, 7, 8])

ns.dp_dotprod(0, u[0], u[1], u[2], 10)
ns.dp_dotprod(3, u[0], u[1], u[2], 11)
ns.dp_dotprod(6, u[0], u[1], u[2], 12)
ns.dp_mul(11, 1.0/3)
ns.dp_mul(12, 2.0/3)
ns.dp_adddp(11, 12, 13)
ns.dp_save('syf_hysteresis_loop.txt', [10, 13])

loop = ns.Get_Data_Columns('syf_hysteresis_loop.txt', [0, 1])
plt.axes(xlabel = 'H (A/m)', ylabel = 'M (A/m)', title = 'SyF Hysteresis Loop')
plt.plot(loop[0], loop[1])
plt.show()


    
    

    


    



























