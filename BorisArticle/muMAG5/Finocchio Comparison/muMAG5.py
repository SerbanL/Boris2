"""
This script is part of Boris Computational Spintronics Article

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient
import numpy as np
import matplotlib.pyplot as plt
from itertools import product

#setup communication with server
ns = NSClient()
ns.configure(True)

########################################

#constant
gmub_2e = 5.795094E-05

########################################

#material parameters
Ms = 8e5
A = 1.3e-11
damping = 0.1
gamma_boris = 2.212761569e5
gamma = 2.211e5
grel = gamma / gamma_boris

########################################

skip_initialization = False

if not skip_initialization:

    #First setup vortex
    
    #dimensions
    meshrect = [100e-9, 100e-9, 10e-9]
    
    #cellsize
    h = 2.5e-9
    
    #set material parameters
    ns.setparam('permalloy', 'Ms', Ms)
    ns.setparam('permalloy', 'A', A)
    ns.setparam('permalloy', 'damping', damping)
    ns.setparam('permalloy', 'grel', grel)
    P = ns.setparam('permalloy', 'P')
    
    #set dimensions
    ns.meshrect(meshrect)
    ns.cellsize([h, h, h])
    
    Nx, Ny, Nz = int(meshrect[0] / h), int(meshrect[1] / h), int(meshrect[2] / h)
    M = [[0,0,0]] * Nx*Ny*Nz
    
    x0 = meshrect[0] / 2
    y0 = meshrect[1] / 2
    z0 = meshrect[2] / 2
    R = 10e-9
    
    #generate m
    for i, j, k in product(range(Nx), range(Ny), range(Nz)):
                
        x, y, z = x0 - (i + 0.5) * h, y0 - (j + 0.5) * h, z0 - (k + 0.5) * h
        f = [y, -x, R]
        fmag = np.sqrt(f[0]**2 + f[1]**2 + f[2]**2)
        M[i + j*Nx + k*Nx*Ny] = [f[0] / fmag, f[1] / fmag, f[2] / fmag]
    
    #load m into Boris, normalizing to Ms
    fileName = 'vortex.ovf'
    ns.Write_OVF2(fileName, M, [Nx, Ny, Nz], [0.0, 0.0, 0.0, Nx*h, Ny*h, Nz*h])
    ns.loadovf2mag(Ms, fileName)
    
    #relax loaded vortex
    ns.editstagestop(0, 'mxh', 1e-5)
    ns.cuda(1)
    ns.Run()
    ns.reset()
    
    #configure for STT
    ns.setode('LLG-STT', 'RK4')
    ns.setdt(100e-15)
    ns.addmodule('permalloy', 'transport')
    ns.ecellsize([h, h, h])
    ns.setdefaultelectrodes()
    
    #simulate for 10ns, saving every 100ps
    ns.editstagestop(0, 'time', 10e-9)
    ns.editdatasave(0, 'time', 100e-12)
    
    #save time and average M
    ns.setdata('time')
    ns.adddata('<M>')
    
    #don't need to update display often
    ns.iterupdate(1000)
    
    #save this state so we don't have to do it again
    ns.savesim('vortex_relaxed')

########################################

#from u and beta return current value we need to set for transport module
def get_I(u, beta):
    
    #current density
    J = u * Ms * (1 + beta**2) / (gmub_2e * 0.4)
    #current to set
    return J * meshrect[1] * meshrect[2]

########################################

def plot_solution(beta):

    savefile = 'muMAG5_beta_%.2f.txt' % beta
    
    data = ns.Get_Data_Columns(savefile, [0, 1, 2])
    time_ns = [t / 1e-9 for t in data[0]]
    Mx = [M for M in data[1]]
    My = [M for M in data[2]]

    solution = ns.Get_Data_Columns('solution_beta_%.2f.txt' % beta, [0, 1, 2])
    
    plt.axes(xlabel = 'time (ns)', ylabel = 'M (A/m)', title = 'beta = %0.2f' % beta)
    plt.plot(time_ns, Mx, 'o', label = 'Mx')
    plt.plot(solution[0], solution[1], '--')
    
    plt.plot(time_ns, My, 'o', label = 'My')
    plt.plot(solution[0], solution[2], '--', label = 'Finocchio et al.')
    plt.legend()
    plt.savefig('solution_beta_%.2f.png' % beta, dpi = 600)
    plt.show()
    
########################################

#u value set in LLG-STT in Boris
#u = J*P*g*muB / (2*e*Ms*(1 + beta^2))

#u and beta values to set
    
#This is what Finocchio et al used as stated in their submitted solution:
u_beta = [[-72.45, 0.0], [-72.45, 0.05], [-72.45, 0.1], [-72.45, 0.5]]

for vals in u_beta:

    ns.loadsim('vortex_relaxed')
    
    u = vals[0]
    beta = vals[1]
    
    I = get_I(u, beta)
    
    ns.setcurrent(I)
    ns.setparam('permalloy', 'beta', beta)
    
    #output data file
    savefile = 'muMAG5_beta_%.2f.txt' % beta
    ns.savedatafile(savefile)
    
    #make sure current density is set accurately
    ns.tsolverconfig(1e-9, 10000)
    ns.computefields()
    #now we have J set we can disable transport solver iteration
    ns.statictransportsolver(1)
    
    #compute vortex orbit
    ns.cuda(1)
    ns.Run()

    plot_solution(beta)



