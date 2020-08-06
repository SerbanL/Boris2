"""
This script is part of Boris Computational Spintronics Article

Generates data in Figure 2.

@author: Serban Lepadatu, 2020
"""

import os
import sys
from NetSocks import NSClient
import matplotlib.pyplot as plt
import matplotlib as mpl

#setup communication with server
ns = NSClient('localhost')

########################################

#plots customizations; see https://matplotlib.org/3.2.1/tutorials/introductory/customizing.html

#axes
mpl.rcParams['axes.linewidth'] = 2.0
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['font.family'] = 'Arial'

#ticks
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15

mpl.rcParams['xtick.major.size'] = 6.0
mpl.rcParams['xtick.major.width'] = 2.0
mpl.rcParams['xtick.minor.size'] = 6.0
mpl.rcParams['xtick.minor.width'] = 2.0

mpl.rcParams['ytick.major.size'] = 6.0
mpl.rcParams['ytick.major.width'] = 2.0
mpl.rcParams['ytick.minor.size'] = 6.0
mpl.rcParams['ytick.minor.width'] = 2.0

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

mpl.rcParams['xtick.major.pad'] = 5.0
mpl.rcParams['xtick.minor.pad'] = 5.0

#legend
mpl.rcParams['patch.linewidth'] = 2.0
mpl.rcParams['legend.fontsize'] = 12.0
mpl.rcParams['legend.edgecolor'] = 'black'
mpl.rcParams['legend.title_fontsize'] = 15
mpl.rcParams['legend.labelspacing'] = 0.1
mpl.rcParams['legend.borderpad'] = 0.3

#math type
mpl.rcParams['mathtext.default'] = 'regular'

mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['savefig.bbox'] = 'tight'

########################################

#the working directory : same as this script file
directory = os.path.dirname(sys.argv[0]) + "/"
#restore program to default state
ns.default()
ns.chdir(directory)

########################################

mumag4sim_field1 = 'muMAG4 field1'
mumag4sim_field2 = 'muMAG4 field2'

ns.loadsim(mumag4sim_field1)
datafile1 = ns.savedatafile()
Ms = ns.setparam('permalloy', 'Ms')
ns.Run()

ns.loadsim(mumag4sim_field2)
datafile2 = ns.savedatafile()
ns.Run()

########################################

data = ns.Get_Data_Columns(datafile1, [0, 1, 2, 3])
time_ns = [time / 1e-9 for time in data[0]]
mx = [M / Ms for M in data[1]]
my = [M / Ms for M in data[2]]
mz = [M / Ms for M in data[3]]

plt.axes(xlabel = 'time (ns)', ylabel = 'M / Ms', title = 'muMAG4 field 1')
plt.plot(time_ns, mx, '.', label = 'mx', color = 'r')
plt.plot(time_ns, my, '.', label = 'my', color = 'b')
plt.plot(time_ns, mz, '.', label = 'mz', color = 'g')
plt.legend()
plt.show()

########################################

data = ns.Get_Data_Columns(datafile2, [0, 1, 2, 3])
time_ns = [time / 1e-9 for time in data[0]]
mx = [M / Ms for M in data[1]]
my = [M / Ms for M in data[2]]
mz = [M / Ms for M in data[3]]

plt.axes(xlabel = 'time (ns)', ylabel = 'M / Ms', title = 'muMAG4 field 2')
plt.plot(time_ns, mx, '.', label = 'mx', color = 'r')
plt.plot(time_ns, my, '.', label = 'my', color = 'b')
plt.plot(time_ns, mz, '.', label = 'mz', color = 'g')
plt.legend()
plt.show()

########################################
    
