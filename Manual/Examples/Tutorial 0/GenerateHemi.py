"""
This script is part of Boris Computational Spintronics v3.0

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient
import numpy as np
from itertools import product

ns = NSClient()
ns.configure(True)

########################################

#Mesh with Nxy cells along x and y, and Nz cells along z.
#Easier to generate mesh this way but the OVF2 file can be loaded in an arbitrarily shaped mesh in Boris (mapped to dimensions).
Nxy, Nz = 32, 16

#M list to write in OVF2 file : has Nxy*Nxy*Nz cells, initialised with empty cells
M = np.zeros((Nxy**2*Nz,3))

#setup hollow hemisphere values
origin = [Nxy/2, Nxy/2, Nz]
inner_ratio, outer_ratio = 0.75, 1.0
rad_inner, rad_outer = Nxy * inner_ratio / 2, Nxy * outer_ratio / 2

for i, j, k in product(range(Nxy), range(Nxy), range(Nz)):
    
    rad = np.sqrt((i - origin[0])**2 + (j - origin[1])**2 + (k - origin[2])**2)
    
    if rad >= rad_inner and rad <= rad_outer:
        #Mark these cells as non-empty
        M[i + j*Nxy + k*Nxy*Nxy] = [1.0,0.0,0.0]

#Write M to OVF2 file ready to load into Boris
#to load into Boris use loadovf2mag command, e.g. as 
#loadovf2mag 8e5 HHemi 
#(specified 8e5 to renormalize to Ms = 8e5 A/m)
fileName = 'HHemi.ovf'
#cellsize used to generate mesh rectangle
h = 5e-9
ns.Write_OVF2(fileName, M, [Nxy, Nxy, Nz], [0.0, 0.0, 0.0, Nxy*h, Nxy*h, Nz*h])

#can also load HHemi through this Python script if Boris is on (set mesh dimensions in Boris first as required)
ns.meshrect([0.0, 0.0, 0.0, Nxy*h, Nxy*h, Nz*h])
ns.loadovf2mag(8e5, fileName)