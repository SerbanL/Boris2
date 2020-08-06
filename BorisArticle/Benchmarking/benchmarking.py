"""
This script is part of Boris Computational Spintronics Article

Generates data in Figure 11.

@author: Serban Lepadatu, 2020
"""

import os
import sys
from NetSocks import NSClient
import numpy as np

#setup communication with server
ns = NSClient('localhost')

########################################

#the working directory : same as this script file
directory = os.path.dirname(sys.argv[0]) + "/"
#restore program to default state
ns.default()
ns.chdir(directory)

ns.loadsim('benchmarking')
h = ns.cellsize()

#Use this for Euler
#ns.setode('LLG', 'Euler')
#ns.setdt(1e-15)
#bench_file = 'Boris_benchmark_Euler_GTX980Ti.txt'
#evals_per_iter = 1

#Use this for RK4
ns.setode('LLG', 'RK4')
ns.setdt(100e-15)
bench_file = 'Boris_benchmark_RK4_GTX980Ti.txt'
evals_per_iter = 4

Nx0, Ny0 = 64, 512
iters0 = 25600

for z in range(7):

    Nz = 2 ** z
    Ny = Ny0 / Nz
    Nx = Nx0
    
    for x in range(0,9):
        
        Nx = Nx0 * (2**x)
        iters = iters0 / (2**x)
        
        ns.meshrect([Nx*h[0], Ny*h[1], Nz*h[2]])
        ns.cellsize(h)
        ns.setangle(90, 0)
        ns.setfield(0.01*1e7/(4*np.pi), 90, 90)
        
        ns.reset()
        ns.editstagestop(0, 'iter', iters)
        
        ns.cuda(1)
        ns.iterupdate(0)
        ns.Run()
        
        benchtime = ns.benchtime()
        
        ns.SaveDataToFile(bench_file, [Nx, Ny, Nz, Nx*Ny*Nz, iters * evals_per_iter, benchtime, benchtime / (iters * evals_per_iter)])


    
    
    

    




