"""
This script is part of Boris Computational Spintronics v3.0

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit

#setup communication with server
ns = NSClient(); ns.configure(True)

########################################

def Boltzmann_Distribution(mA, prop):
    
    a = -Ms0_A*V / (4 * mu_A * me_A * kB * T)
    
    b = (mA**2 - me_A**2)**2 / me_A
    c = (mu_A + 3*tau_AB*kB*Tn*susrel_B) / (2*susrel_A)
    
    d = (mB**2 - me_B**2) * 3*tau_AB * kB * Tn * mA**2 / me_B
    
    return prop * mA**2 * np.exp(a * (b * c + d))

########################################

#useful constants
kB = 1.3806488e-23
muB = 9.27400968e-24

#Neel and set temperature
Tn = 500
T = 495

tau = [0.5, 0.5, 0.5, 0.5]
mu = [1, 1]
Ms0_AFM = [800e3, 800e3]

########################################

#setup antiferromagnetic mesh
ns.setafmesh('AF', [320e-9, 320e-9, 320e-9])

h = [5e-9, 5e-9, 5e-9]
hs = [5e-9, 5e-9, 5e-9]
V = hs[0]*hs[1]*hs[2]

dT = 2.5e-15
dTstochastic = 10e-15

relax_time = 20e-12

ns.cellsize(h)
ns.scellsize(hs)

ns.setode('sLLB', 'TEuler')
ns.setdt(dT)
ns.setdtstoch(dTstochastic)

ns.editstagestop(0, 'time', relax_time)
ns.curietemperature(Tn)
ns.temperature(T)

ns.tau(tau)
ns.atomicmoment(mu, 'AF')
ns.setparam('AF', 'Ms_AFM', Ms0_AFM)

########################################

#get parameter values

#1. Ms0
#setparam always gets the zero value temperature
#Ms0_AFM = ns.setparam('AF', 'Ms_AFM')
Ms0_A = Ms0_AFM[0]
Ms0_B = Ms0_AFM[1]

#2. get equilibrium magnetisation at set temperature
ns.dp_dumptdep('AF', 'Ms_AFM', Tn, 0)
me_A = ns.dp_get(0, int(T))
me_B = ns.dp_get(1, int(T))

#3. longitudinal relative susceptibility
ns.dp_dumptdep('AF', 'susrel_AFM', Tn, 0)
susrel_A = ns.dp_get(0, int(T))
susrel_B = ns.dp_get(1, int(T))

#4. atomic moments
#mu = ns.atomicmoment()
mu_A = mu[0] * muB
mu_B = mu[1] * muB

#5. tau values (returned as tau11 tau22, tau12 tau21)
#tau = ns.tau()
tau_AB = tau[2]

########################################

#simulate to relax magnetisation for set temperature

ns.cuda(1)
ns.Run()

########################################

#get histogram for sub-lattice A

mA_low = 0.6 * me_A
mA_high = 1.3 * me_A
steps_A = 120

mB_low = 0.6 * me_B
mA_high = 1.3 * me_B
steps_B = 120

mB_range = np.arange(mB_low, mA_high, (mA_high-mB_low)/steps_B)
mA_range = np.array([])

Z = np.array([[]])
Zeq = np.array([[]])

for mB in mB_range:

    ns.dp_histogram2(0, 1, ( mA_high - mA_low) * Ms0_A / steps_A, Ms0_A * mA_low, Ms0_A * mA_high, Ms0_B * mB, Ms0_B * mB * 0.01)
    ns.dp_save('histo.txt', [0, 1])
    histo = ns.Get_Data_Columns('histo.txt', [0, 1])
    
    #computed histogram data for sub-lattice A
    if mA_range.size == 0:
        mA_range = np.array([M / Ms0_A for M in histo[0]])
    
    P_Histogram = [P for P in histo[1]]
    
    #fit proportionality constant for Boltzmann distribution to computed distribution
    prop_guess = np.max(P_Histogram) / Boltzmann_Distribution(mB, 1)
    popt, pcov = curve_fit(Boltzmann_Distribution, mA_range, P_Histogram, p0 = prop_guess)
    P_Boltzmann = [Boltzmann_Distribution(m, popt[0]) for m in mA_range]
    
    if Z.size == 0:
        Z = np.column_stack(P_Histogram)
    else:
        Z = np.vstack((Z, P_Histogram))
        
    if Zeq.size == 0:
        Zeq = np.column_stack(P_Boltzmann)
    else:
        Zeq = np.vstack((Zeq, P_Boltzmann))
    
    #plt.plot(mA_range, P_Histogram, '.')
    #plt.plot(mA_range, P_Boltzmann, '--')
    #plt.show()

X, Y = np.meshgrid(mA_range, mB_range)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlabel('$m_{A}$')
ax.set_ylabel('$m_{B}$')
ax.set_zlabel('Probability')
ax.view_init(30, 110)

surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

fig.colorbar(surf, shrink=0.5, aspect=8)

ax.plot_wireframe(X, Y, Zeq, rstride=5, cstride=5, linewidth = 0.5, linestyle = (0, (5, 1)), color = 'black')

ax.set_zticklabels([])
plt.title('T/Tn = %0.2f' % (T / Tn))

plt.savefig('boltzmann.png', dpi = 600)
plt.show()



    

    


