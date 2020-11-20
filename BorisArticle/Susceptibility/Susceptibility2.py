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

#Curie temperature (K)
TCurie = 565
#Maximum temperature and step (K)
Tmax, Tstep = 1000, 10

#Field to set (A/m)
Hfield = 1e4

#Atomic moment (muB multiple)
#Gd0.25FeCo0.75 material as found in P. Nieves et al., Low Temp. Phys. 41, 739 (2015)
mu = [1.92, 7.63]

#Coupling constants to phase transition temperature as tau11, tau22, tau12, tau22
#tau values calculated from constants in Table 1 of P. Nieves et al., Low Temp. Phys. 41, 739 (2015)

#J0_A = x z J_A, where x = 0.75 is the concentration of species A (TM: FeCo), z is the number of nearest neighbors
#P. Nieves et al., Low Temp. Phys. 41, 739 (2015) gives z*J values, so J0_A is readily calculated
#Then: |J0_A| = 3*kB * tau_A * TCurie gives tau_A = 0.957

#Similarly J0_B = (1-x) z J_B, |J0_B| = 3*kB * tau_B * TCurie, gives tau_B = 0.127 (B is RE: Gd)

#J0_AB = (1-x) z J_AB, where again x is A concentration, so 1-x = 0.25.
#|J0_AB| = 3*kB * tau_AB * TCurie, so tau_AB = 0.111

#J0_BA = x z J_AB
#|J0_BA| = 3*kB * tau_BA * TCurie, so tau_BA = 0.333
tau = [0.958, 0.127, 0.111, 0.333]

me_input = 'me2_input.txt'
susrel_input = 'susrel2_input.txt'
damping_input = 'damping2_input.txt'
me_output = 'me2_output.txt'
susceptibility_figure_name = 'susceptibility2.png'
damping_figure_name = 'damping2.png'

########################################
# Compute me output function

meshname = '2SL'
ns.setafmesh(meshname, [10e-9, 10e-9, 10e-9])
ns.cellsize([5e-9, 5e-9, 5e-9])

ns.curietemperature(TCurie)

meshname = ns.meshfocus()
Ms_AFM = ns.setparam(meshname, 'Ms_AFM')
ns.tau(tau)
ns.atomicmoment(mu, meshname)

ns.setode('LLB', 'RKF45')
ns.setdt(100e-15)

ns.setfield(Hfield, 90, 0)
ns.setangle(90, 0)

chunkiters = 250
iters_timeout = 25000
ns.editstagestop(0, 'iter', chunkiters)

#wipe output file clean
ns.dp_newfile(me_output)

def Run_Chunk():
    
    ns.reset()
    ns.Run()
    Mav = ns.showdata('<M>')
    Mav2 = ns.showdata('<M2>')
    m1_value = np.sqrt(Mav[0]**2 + Mav[1]**2 + Mav[2]**2) / Ms_AFM[0]
    m2_value = np.sqrt(Mav2[0]**2 + Mav2[1]**2 + Mav2[2]**2) / Ms_AFM[1]
    return m1_value, m2_value
    
m1_value_previous, m2_value_previous, m1_value, m2_value = 1.0, 1.0, 0.0, 0.0
for T in np.arange(0, Tmax + Tstep, Tstep):
    
    ns.temperature(T)
    ns.setdt(10e-15)
    
    iters = 0
    while iters < iters_timeout:
        
        m1_value, m2_value = Run_Chunk()
        if m1_value < m1_value_previous or m2_value < m2_value_previous: 
            m1_value_previous = m1_value
            m2_value_previous = m2_value
        else: break
    
        iters += chunkiters
        
    ns.SaveDataToFile(me_output, [T, m1_value, m2_value])

########################################
# Plot me and susceptibilities

#input me function
ns.dp_dumptdep(meshname, 'Ms_AFM', Tmax, 0)
ns.dp_save(me_input, [0, 1])
me2 = ns.Get_Data_Columns(me_input, [0, 1])

#input longitudinal relative susceptibility
ns.dp_dumptdep(meshname, 'susrel_AFM', Tmax, 0)
ns.dp_save(susrel_input, [0, 1])
susrel2 = ns.Get_Data_Columns(susrel_input, [0, 1])

Trange = np.arange(0, Tmax + 1, 1)
    
#output me function
data_computed = ns.Get_Data_Columns(me_output, [0, 1, 2])

ax = plt.axes(xlabel = 'Temperature (K)', ylabel = 'm')
ax2 = ax.twinx()
ax2.set_ylabel(r'$\tilde{\chi}$ (1/T)')

lns1 = ax.plot(data_computed[0], data_computed[1], 'ob', label = 'Output m$_{eA}$')
lns2 = ax.plot(data_computed[0], data_computed[2], 's', color = 'darkviolet', label = 'Output m$_{eB}$')
lns3 = ax.plot(Trange, me2[0], '-g', label = 'Input m$_{eA}$')
lns4 = ax.plot(Trange, me2[1], '--', color = 'orange', label = 'Input m$_{eB}$')
lns5 = ax2.plot(Trange, susrel2[0], '-r', label = r'$\tilde{\chi}_{A}$')
lns6 = ax2.plot(Trange, susrel2[1], '-.', color = 'chocolate', label = r'$\tilde{\chi}_{B}$')
ax2.set_yscale('log')

locmaj = matplotlib.ticker.LogLocator(base=10,numticks=10)
ax2.yaxis.set_major_locator(locmaj)
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=10)
ax2.yaxis.set_minor_locator(locmin)
ax2.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

plt.axvline(TCurie, ls = '--', color = 'black')

lns = lns1+lns2+lns3+lns4+lns5+lns6
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc = 'lower right')
plt.savefig(susceptibility_figure_name)
plt.show()

########################################
# Plot damping

muB = 9.274009994e-24
kB = 1.38064852E-23
TCurie_renorm = 2*TCurie / (tau[0] + tau[1] + np.sqrt((tau[0]-tau[1])**2 + 4*tau[2]*tau[3]))
print("Renormalized Curie temperature: %f K. If setup correctly should be same as Curie temperature (%f K)!" % (TCurie_renorm, TCurie))

ns.dp_dumptdep(meshname, 'damping_AFM', Tmax, 0)
ns.dp_save(damping_input, [0, 1])
damping2 = ns.Get_Data_Columns(damping_input, [0, 1])

ax = plt.axes(xlabel = 'Temperature (K)', ylabel = 'Damping Scaling')
lns1 = ax.plot(Trange, damping2[0], '-g', label = 'A')
lns2 = ax.plot(Trange, damping2[1], '--', color = 'orange', label = 'B')
plt.axvline(TCurie, ls = '--', color = 'black')

lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc = 'lower right')
plt.savefig(damping_figure_name)
plt.show()
    
    
    
    
    




