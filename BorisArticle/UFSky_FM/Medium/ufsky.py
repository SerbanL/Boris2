"""
This script is part of Boris Computational Spintronics v2.8

@author: Serban Lepadatu, 2020
"""

import os
import sys
from NetSocks import NSClient
import numpy as np
import matplotlib.pyplot as plt

#setup communication with server. By default sent messages are not displayed in console. 
#To enable verbose mode use e.g.: NSClient('localhost', True)
ns = NSClient('localhost')

########################################

#the working directory : same as this script file, typically expecting simulation file to be in same directory as this script file
directory = os.path.dirname(sys.argv[0]) + "/"
#make sure Boris is reset to default in case a ready set simulation file is not loaded below (if so this can be commented out)
ns.default()
#set working directory same as this script file
ns.chdir(directory)

########################################

def run_laser_pulse(pulse_num):

    #for a script running this program entirely from scratch (no loadsim) see Tutorial 0
    ns.loadsim('ufsky_fm')
    
    data_file = 'ufsky_%d.txt' % pulse_num
    ns.savedatafile(data_file)
    
    #0 to 10ps
    ns.editstagestop(0, 'time', 10e-12)
    ns.setdt(1e-15)
    ns.setheatdt(1e-15)
    ns.Run()
    
    #10ps to 20ps
    ns.editstagestop(0, 'time', 20e-12)
    ns.setdt(2e-15)
    ns.setheatdt(2e-15)
    ns.Run()
    
    #20ps to 60ps
    ns.editstagestop(0, 'time', 60e-12)
    #mid-simulation re-meshing! 1nm cellsize only needed around the Curie temperature.
    #finely tuned simulations of this type means you can run thousands of these events in a reasonable time-scale and still keep accuracy
    ns.cellsize([2e-9, 2e-9, 2e-9])
    ns.setdt(10e-15)
    ns.setheatdt(2e-15)
    ns.Run()
    
    #60ps to 800ps
    ns.editstagestop(0, 'time', 800e-12)
    ns.tcellsize([4e-9, 4e-9, 1e-9])
    ns.setdt(50e-15)
    #could do with Crank-Nicolson method in next version
    ns.setheatdt(5e-15)
    ns.editdatasave(0, 'time', 250e-15)
    ns.Run()
    
    ns.savemeshimage('ufsky_final_%d' % pulse_num)
    
    #now plot |Q| as a function of time
    data = ns.Get_Data_Columns(data_file, [0, 1])
    time_ps = [t/1e-12 for t in data[0]]
    Qmod = [np.abs(Qval) for Qval in data[1]]
    
    plt.axes(xlabel = 'Time (ps)', ylabel = '|Q|')
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(time_ps, Qmod)
    plt.xlim(0.1)
    plt.ylim(1e-3)
    plt.savefig('ufsky_plot_%d.png' % pulse_num, dpi = 600)
    plt.show()
    
    Qval = ns.dp_topocharge()
    
    return Qval, np.round(Qval)

########################################

Qvalues_file = 'Qcount.txt'

for pulse in range(0, 1):
    
    Q, Qrounded = run_laser_pulse(pulse)
    
    ns.SaveDataToFile(Qvalues_file, [pulse, Q, Qrounded])
