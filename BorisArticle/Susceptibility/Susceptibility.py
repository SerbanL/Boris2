"""
This script is part of Boris Computational Spintronics Article

@author: Serban Lepadatu, 2020
"""

import os
from NetSocks import NSClient
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
import matplotlib.ticker

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

mpl.rcParams['ytick.major.pad'] = 5.0
mpl.rcParams['ytick.minor.pad'] = 5.0

#legend
mpl.rcParams['patch.linewidth'] = 2.0
mpl.rcParams['legend.fontsize'] = 14
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
directory = os.path.dirname(os.path.realpath(__file__)) + "/"
#restore program to default state
ns.default()
ns.chdir(directory)

########################################
# Settings 

#Curie temeprature (K)
TCurie = 800
#Maximum temperature and step (K)
Tmax, Tstep = 1000, 10

#Atomic moment (muB multiple)
mu = 1

#Field to set (A/m)
Hfield = 1e4

me_input = 'me_input.txt'
susrel_input = 'susrel_input.txt'
me_output = 'me_output.txt'
figure_name = 'susceptibility.png'

########################################
# Compute me output function

ns.meshrect([10e-9, 10e-9, 10e-9])
ns.cellsize([5e-9, 5e-9, 5e-9])

ns.curietemperature(TCurie)
ns.atomicmoment(mu)

ns.setode('LLB', 'RKF45')
ns.setdt(100e-15)

ns.setfield(Hfield, 90, 0)
ns.setangle(90, 0)

chunkiters = 250
iters_timeout = 25000
ns.editstagestop(0, 'iter', chunkiters)

meshname = ns.meshfocus()
Ms = ns.setparam(meshname, 'Ms')
ns.delmodule(meshname, 'demag')

#wipe output file clean
ns.dp_newfile(me_output)

def Run_Chunk():
    
    ns.reset()
    ns.Run()
    Mav = ns.showdata('<M>')
    m_value = np.sqrt(Mav[0]**2 + Mav[1]**2 + Mav[2]**2) / Ms
    return m_value
    
m_value_previous, m_value = 1.0, 0.0
for T in np.arange(0, Tmax + Tstep, Tstep):
    
    ns.temperature(T)
    ns.setdt(10e-15)
    
    iters = 0
    while iters < iters_timeout:
        
        m_value = Run_Chunk()
        if m_value < m_value_previous: m_value_previous = m_value
        else: break
    
        iters += chunkiters
    
    ns.SaveDataToFile(me_output, [T, m_value])

########################################
# Plot everything

#input me function
ns.dp_dumptdep(meshname, 'Ms', Tmax, 0)
ns.dp_save(me_input, 0)
me = ns.Get_Data_Columns(me_input, 0)

#input longitudinal relative susceptibility
ns.dp_dumptdep(meshname, 'susrel', Tmax, 0)
ns.dp_save(susrel_input, 0)
susrel = ns.Get_Data_Columns(susrel_input, 0)

Trange = np.arange(0, Tmax + 1, 1)

#output me function
data_computed = ns.Get_Data_Columns(me_output, [0, 1])

ax = plt.axes(xlabel = 'Temperature (K)', ylabel = 'm')
ax2 = ax.twinx()
ax2.set_ylabel(r'$\tilde{\chi}$ (1/T)')

lns1 = ax.plot(data_computed[0], data_computed[1], 'ob', label = 'Output m$_{e}$')
lns2 = ax.plot(Trange, me, '--g', label = 'Input m$_{e}$')
lns3 = ax2.plot(Trange, susrel, '-r', label = r'$\tilde{\chi}$')
ax2.set_yscale('log')

locmaj = matplotlib.ticker.LogLocator(base=10,numticks=10)
ax2.yaxis.set_major_locator(locmaj)
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=10)
ax2.yaxis.set_minor_locator(locmin)
ax2.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

plt.axvline(TCurie, ls = '--', color = 'black')

lns = lns1+lns2+lns3
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc = 'lower left')
plt.savefig(figure_name)
plt.show()
    



    
    
    
    
    
    




