"""
This script is part of Boris Computational Spintronics Article

@author: Serban Lepadatu, 2020
"""

import os
from NetSocks import NSClient
import numpy as np
from numpy import random
from itertools import product

#setup communication with server
ns = NSClient('localhost')

########################################

#the working directory : same as this script file
directory = os.path.dirname(os.path.realpath(__file__)) + "/"
#restore program to default state
ns.default()
ns.chdir(directory)

########################################
# Simulation function

def Simulate_SpinIce(d, w, l, h, t, max_dim, results_fileName):

    u = int(np.floor(max_dim / d))
    if u % 2 == 1: u -= 1
    
    u_x, u_y = u, u
    hx, hy, hz = h[0], h[1], h[2]
    
    #spacing (should be > w)
    s = d - l

    fileName = 'SpinIce.ovf'
    output_fileName = 'Relaxed_' + fileName
    
    ########################################
    # Make spin ice structure
     
    #Number of cells in mesh and M
    Nx, Ny, Nz = int(np.round(u_x * d / hx)), int(np.round(u_y * d / hy)), int(np.round(t / hz))
    M = np.zeros((Nx*Ny*Nz,3))
    
    #number of cells in unit, and unit list
    nx, ny, nz = int(np.round(d/hx)), int(np.round(d/hy)), int(np.round(t/hz))
    U = np.zeros((nx*ny*nz,3))
    
    Dirs = np.zeros((u_x*u_y,2))
    for i, j in product(range(u_x), range(u_y)):
        
        horiz_dir = random.randint(2) * 2 - 1
        vert_dir = random.randint(2) * 2 - 1
        Dirs[i + j * u_x] = [horiz_dir, vert_dir]
    
    def Make_Unit(x, y):
        
        horiz_dir = Dirs[x + y * u_x][0]
        vert_dir = Dirs[x + y * u_x][1]
        horizcap_dir = Dirs[(x + 1)%u_x + y * u_x][0]
        vertcap_dir = Dirs[x + ((y + 1)%u_y) * u_x][1]
        
        for i, j, k in product(range(nx), range(ny), range(nz)):
            
            px, py = (i + 0.5) * hx, (j + 0.5) * hy
            idx = i + j * nx + k * nx*ny
            
            #horizontal rectangle
            if px < l and py >= d - s/2 - w/2 and py < d - s/2 + w/2: U[idx] = [horiz_dir,0,0]
            
            #vertical rectangle
            if py < l and px >= d - s/2 - w/2 and px < d - s/2 + w/2: U[idx] = [0,vert_dir,0]
            
            #horizontal rectangle right cap
            if px >= l and np.sqrt((px - l)**2 + (py - (d - s/2))**2) < w/2: U[idx] = [horiz_dir,0,0]
            
            #horizontal rectangle left cap
            if np.sqrt((px - d)**2 + (py - (d - s/2))**2) < w/2: U[idx] = [horizcap_dir,0,0]
            
            #vertical rectangle top cap
            if py >= l and np.sqrt((py - l)**2 + (px - (d - s/2))**2) < w/2: U[idx] = [0,vert_dir,0]
            
            #vertical rectangle bottom cap
            if np.sqrt((py - d)**2 + (px - (d - s/2))**2) < w/2: U[idx] = [0,vertcap_dir,0]
    
    #Now copy over unit to M the indicated number of times in each direction
    def Copy_Unit(nsx, nsy, M, U):
        
        for i, j, k in product(range(nx), range(ny), range(nz)):
            
            I, J, K = nsx + i, nsy + j, k        
            M[I + J * Nx + K * Nx*Ny] = U[i + j * nx + k * nx*ny]
            
        return M
    
    Make_Unit(0, 0)
    
    for i, j in product(range(u_x), range(u_y)):
        
        nsx, nsy = i * nx, j * ny
        
        #Make_Unit(i, j)
        M = Copy_Unit(nsx, nsy, M, U)
        
    ########################################
    # Set spin ice structure and relax it
        
    ns.Write_OVF2(fileName, M, [Nx, Ny, Nz], [0.0, 0.0, 0.0, Nx*hx, Ny*hy, Nz*hz])
    
    ns.cuda(0)
    ns.setstage('Relax')
    meshName = ns.meshfocus()
    Ms = ns.setparam(meshName, 'Ms')
    
    ns.cellsize([hx, hy, hz])
    ns.meshrect([0.0, 0.0, 0.0, Nx*hx, Ny*hy, Nz*hz])
    ns.cellsize([hx, hy, hz])
    ns.loadovf2mag(Ms, fileName)
    
    ns.pbc(meshName, 'x', 10)
    ns.pbc(meshName, 'y', 10)
    
    ns.displaydetail(500e-9)
    ns.vecrep(meshName, 4)
    
    #Relax from random state : first stage without SDesc
    ns.random()
    ns.setode('LLGStatic', 'AHeun')
    ns.editstagestop(0, 'iter', 100)
    ns.cuda(1)
    ns.Run()
    
    #Second relaxation stage : SDesc to a not so low mxh value (relaxing to very low mxh values doesn't make any difference to the final configuration as it tends to "freeze" at relatively high mxh values)
    ns.setode('LLGStatic', 'SDesc')
    ns.editstagestop(0, 'mxh', 1e-2)
    ns.reset()
    ns.Run()
    
    #At this point there will be some islands not in single domain, but there will be a bias towards one direction or another
    #In practice these are removed using a rotating magnetic field, but this is not feasible in a micromagnetic simulation due to prohibitively long simulation times.
    #Instead we force them into single domain by computing the current preferred direction.
    ns.saveovf2mag(output_fileName)
    
    #Third stage: read M in and check for any islands not in single-domain state: if any found force them into single domain state depending on current bias
    M, n, meshrect = ns.Read_OVF2(output_fileName)
    
    def Make_SingleDomain(M, Ms, sx, sy, ex, ey, ez):
        
        M_av = 0.0
        total = (ex - sx) * (ey - sy) * ez
        
        def Set_Island_Direction(M, Value, sx, sy, ex, ey, ez):
            
            for i, j, k in product(range(ex - sx), range(ey - sy), range(ez)):
        
                idx = sx + i + (sy + j) * Nx + k * Nx*Ny
                M[idx] = Value
        
        for i, j, k in product(range(ex - sx), range(ey - sy), range(ez)):
        
            idx = sx + i + (sy + j) * Nx + k * Nx*Ny
            
            #horizontal : x component
            if (ex-sx) > (ey-sy): M_av += M[idx][0] / total
            #vertical: y component
            else: M_av += M[idx][1] / total
            
        if M_av / Ms < 1/2 and M_av / Ms > -1/2:
            
            #not purely deterministic
            m = (M_av / Ms)*12
            p = random.random()*2 - 1
            
            if (ex-sx) > (ey-sy):
                if m > p: Set_Island_Direction(M, [Ms, 0, 0], sx, sy, ex, ey, ez)
                else: Set_Island_Direction(M, [-Ms, 0, 0], sx, sy, ex, ey, ez)
            else:
                if m > p: Set_Island_Direction(M, [0, Ms, 0], sx, sy, ex, ey, ez)
                else: Set_Island_Direction(M, [0, -Ms, 0], sx, sy, ex, ey, ez)
    
    for i, j in product(range(u_x), range(u_y)):
    
        #positition of vertex center
        px, py = i * d + l + s/2 + hx/2, j * d + l + s/2 + hy/2
        
        Make_SingleDomain(M, Ms, 
                          int((px - s/2 - l) / hx), int((py - w/2) / hy),
                          int((px - s/2) / hx), int((py + w/2) / hy), int(np.round(t / hz)))
        
        Make_SingleDomain(M, Ms, 
                          int((px - w/2) / hx), int((py - s/2 - l) / hy),
                          int((px + w/2) / hx), int((py - s/2) / hy), int(np.round(t / hz)))
        
    #Write M, then relax again
    ns.Write_OVF2(fileName, M, [Nx, Ny, Nz], [0.0, 0.0, 0.0, Nx*hx, Ny*hy, Nz*hz])
    ns.loadovf2mag(Ms, fileName)
    
    ns.setode('LLGStatic', 'RK4')
    ns.setdt(1e-12)
    ns.editstagestop(0, 'mxh', 5e-2)
    ns.reset()
    ns.Run()

    #Relaxed spin ice state, all in single domain (there might still be exceptions, but negligible)
    ns.saveovf2mag(output_fileName)
            
    ########################################
    # Analyze vertices
    
    M, n, meshrect = ns.Read_OVF2(output_fileName)
    
    #L R T B
    
    #Type 1:
    #+ - + -
    #- + - +
    type_1 = 0
    
    #Type 2:
    #+ + + +
    #+ + - -
    #- - + +
    #- - - -
    type_2 = 0
    
    #Type 3:
    #- + + +
    #- + - -
    #+ - + +
    #+ - - -
    #+ + + -
    #- - + -
    #+ + - +
    #- - - +
    type_3 = 0
    
    #Type 4:
    #+ - - +
    #- + + -
    type_4 = 0
    
    for i, j in product(range(u_x), range(u_y)):
        
        #positition of vertex center
        px, py = i * d + l + s/2, j * d + l + s/2
        
        #left, right, top, bottom indexes
        idx_l = int(np.floor((px - s/2) / hx + (py / hy) * Nx))
        idx_r = int(np.floor((px + s/2) / hx + (py / hy) * Nx))
        idx_t = int(np.floor(px / hx + ((py + s/2 - hy) / hy) * Nx))
        idx_b = int(np.floor(px / hx + ((py - s/2) / hy) * Nx))
        
        lft_pve = M[idx_l][0] > 0
        rgt_pve = M[idx_r][0] > 0
        top_pve = M[idx_t][1] > 0
        bot_pve = M[idx_b][1] > 0
        
        if ((lft_pve and not rgt_pve and top_pve and not bot_pve) or 
            (not lft_pve and rgt_pve and not top_pve and bot_pve)):
            type_1 += 1
        
        elif ((lft_pve and rgt_pve and top_pve and bot_pve) or 
            (lft_pve and rgt_pve and not top_pve and not bot_pve) or
            (not lft_pve and not rgt_pve and top_pve and bot_pve) or 
            (not lft_pve and not rgt_pve and not top_pve and not bot_pve)):
            type_2 += 1
            
        elif ((lft_pve and not rgt_pve and not top_pve and bot_pve) or
              (not lft_pve and rgt_pve and top_pve and not bot_pve)):
            type_4 += 1
            
        else: type_3 += 1
        
    total = type_1 + type_2 + type_3 + type_4

    print("Type 1 : %0.1f, Type 2 : %0.1f, Type 3 : %0.1f, Type 4 : %0.1f" % 
          (type_1 * 100 / total, type_2 * 100 / total, type_3 * 100 / total, type_4 * 100 / total))
    print("Excess Type 1 : %0.1f, Excess Type 2 : %0.1f, Excess Type 3 : %0.1f, Excess Type 4 : %0.1f" % 
          (type_1 * 100 / total - 12.5, type_2 * 100 / total - 25, type_3 * 100 / total - 50, type_4 * 100 / total - 12.5))
    
    #save lattice spacing, width, length, number of units cells per planar dimension, excess type 1, excess type 2, excess type 3, excess type 4
    ns.SaveDataToFile(results_fileName, [d, w, l, u, type_1 * 100 / total - 12.5, type_2 * 100 / total - 25, type_3 * 100 / total - 50, type_4 * 100 / total - 12.5])

########################################
# Simulation

#width and length
w = 80e-9
l = 220e-9 - w

#cellsize
hx, hy, hz = 5e-9, 5e-9, 25e-9

#thickness
t = 25e-9

#maximum mesh dimension
max_dim = 24e-6

for d_nm in range(300, 950, 50):
    
    #lattice constant
    d = d_nm * 1e-9
    
    Simulate_SpinIce(d, w, l, [hx, hy, hz], t, max_dim, 'SpinIce_Results.txt')
    
    



    
    
    
    
    
    




