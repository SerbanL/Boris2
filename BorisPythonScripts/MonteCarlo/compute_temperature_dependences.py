import os, sys
sys.path.append(os.path.expanduser('~/Documents/Boris Data/'))
sys.path.append(os.path.expanduser('~/Documents/Boris_Data/'))
from NetSocks import NSClient, customize_plots

ns = NSClient(); ns.configure(True); customize_plots()

from MonteCarlo import simulate_m_scaling
from MonteCarlo import simulate_susceptibilities
from MonteCarlo import simulate_k_uniaxial_scaling
from MonteCarlo import simulate_k_cubic_scaling
from MonteCarlo import simulate_m_k_uniaxial_scaling
from MonteCarlo import simulate_m_k_cubic_scaling
from MonteCarlo import simulate_exchange_stiffness
from MonteCarlo import simulate_exchange_stiffness_spring
from MonteCarlo import simulate_RKKY

import numpy as np

######################################################
# Setup mesh, modules, and parameters - USER EDIT

#atomistic magnetic moment as multiple of muB (Bohr magneton)
mu_s = 1.3
#atomistic exchange (units J)
J = 6e-21
#uniaxial anisotropy (units J)
K = 5e-24
#surface exchange coupling (units J)
Js = 6e-21

#mesh side (nm)
mesh_x_nm = 40
mesh_y_nm = 40
mesh_z_nm = 40

#lattice constant (nm)
unit_cell_nm = 0.25

######################################################
# Apply setup

def set_mesh(meshName, rectangle, direction = [90, 0], pbc = [1, 1, 1], method_set_mesh = True):

    if method_set_mesh: ns.setameshcubic(meshName, rectangle)
    else: ns.addameshcubic(meshName, rectangle)
    ns.meshfocus(meshName)
    ns.cellsize([unit_cell_nm*1e-9])
    
    #no demag and set pbc in all directions : change if needed
    ns.delmodule(meshName, 'demag')
    ns.pbc(meshName, 'x', pbc[0])
    ns.pbc(meshName, 'y', pbc[1])
    ns.pbc(meshName, 'z', pbc[2])
    
    #set magnetic moment (units of muB)
    ns.setparam(meshName, 'mu_s', mu_s)
    #set atomistic exchange (units J)
    ns.setparam(meshName, 'J', J)
    #anisotropy constant (units J) - set uniaxial anisotropy by default, but routines below will change this if needed
    ns.addmodule(meshName, 'aniuni')
    ns.setparam(meshName, 'K1', K)
    #damping set to 1 in case sLLG is used
    ns.setparam(meshName, 'damping', 1)
    
    ns.setangle(direction[0], direction[1], meshName)
    
    #don't update display often to save CPU power
    ns.iterupdate(1000)
    
def set_mesh_rkky(meshName, rectangle, direction = [90, 0], method_set_mesh = True):

    set_mesh(meshName, rectangle, direction, [1, 1, 0], method_set_mesh)
    ns.setparam(meshName, 'ea1', [0, 1, 0])
    ns.addmodule(meshName, 'surfexchange')
    theta, phi = np.radians(direction[0]), np.radians(direction[1])
    ns.mcconstrain([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)], meshName)
    ns.setparam(meshName, 'Js', Js)

######################################################
# Simulate. Uncomment as needed.

#Use this mesh setting function for routines below
#set_mesh('asc', [mesh_x_nm*1e-9, mesh_y_nm*1e-9, mesh_z_nm*1e-9])    

#Use this to compute m scaling only.
#simulate_m_scaling(ns, 'mscaling.txt', Tmax = 700, Hext = 0.0)

#Use this to compute susceptibilities - also m scaling obtained.
#simulate_susceptibilities(ns, 'susceptibilities.txt', Tmax = 1000, Hext = 0.0)

#Use these to simulate uniaxial, resp. cubic, anisotropies using CMC (slow, but no assumptions required) - also m scaling obtained.
#simulate_k_uniaxial_scaling(ns, 'CMC_uni.txt', Tmax = 700, hard_axis = [90, 90], easy_axis = [90, 0])
#simulate_k_cubic_scaling(ns, 'CMC_cubi.txt', Tmax = 700, hard_axis = [np.degrees(np.arctan(np.sqrt(2))), 45], easy_axis1 = [90, 0], easy_axis2 = [90, 90])

#Use these to simulate uniaxial, resp. cubic, anisotropies using MC (much faster than CMC, but assumptions made - see function descriptions) - also m scaling obtained.
#simulate_m_k_uniaxial_scaling(ns, 'MC_uniaxial_m_k.txt', Tmax = 1000, hard_axis = [90, 90], easy_axis = [90, 0])
#simulate_m_k_cubic_scaling(ns, 'MC_cubic_m_k.txt', Tmax = 700, hard_axis = [np.degrees(np.arctan(np.sqrt(2))), 45], easy_axis1 = [90, 0], easy_axis2 = [90, 90])

#Use this to compute a scaling - also m and k (uniaxial) scalings obtained. Should set large K value for this : K = 5e-23 J is good.
#simulate_exchange_stiffness(ns, 'ascaling.txt', Tmax = 700, Hext = 0.0)
    
"""
#use this to compute interlayer exchange coupling temperature scaling due to spin-waves
set_mesh_rkky('bottom', [mesh_x_nm*1e-9, mesh_y_nm*1e-9, mesh_z_nm*1e-9], [90, 90])
set_mesh_rkky('top', [0, 0, mesh_z_nm*1e-9, mesh_x_nm*1e-9, mesh_y_nm*1e-9, 2*mesh_z_nm*1e-9], [90, 270], False)
outputFile = 'rkkyscaling_%0.2fnm_Js%0.2e.txt' % (mesh_z_nm, Js)
simulate_RKKY(ns, 'bottom', 'top', outputFile, Tmax = 700, averaging_level = 2, use_H_method = False, skip_normal_method = False)
"""





