from NetSocks import NSClient
from NetSocks import customize_plots

from MonteCarlo import simulate_m_scaling
from MonteCarlo import simulate_susceptibilities
from MonteCarlo import simulate_k_uniaxial_scaling
from MonteCarlo import simulate_k_cubic_scaling
from MonteCarlo import simulate_m_k_uniaxial_scaling
from MonteCarlo import simulate_m_k_cubic_scaling
from MonteCarlo import simulate_exchange_stiffness

ns = NSClient()
ns.configure(True)
customize_plots()

######################################################
# Setup mesh, modules, and parameters - USER EDIT

#atomistic magnetic moment as multiple of muB (Bohr magneton)
mu_s = 1.3
#atomistic exchange (units J)
J = 6e-21
#uniaxial anisotropy (units J)

#anisotropy energy
K = 5e-23

#mesh side (nm)
mesh_x_nm = 40
mesh_y_nm = 40
mesh_z_nm = 40

#lattice constant (nm)
unit_cell_nm = 0.25

######################################################
# Apply setup

meshname = 'asc'
ns.setameshcubic(meshname, [mesh_x_nm*1e-9, mesh_y_nm*1e-9, mesh_z_nm*1e-9])
ns.cellsize([unit_cell_nm*1e-9])

#no demag and set pbc in all directions : change if needed
ns.delmodule(meshname, 'demag')
ns.pbc(meshname, 'x', 1)
ns.pbc(meshname, 'y', 1)
ns.pbc(meshname, 'z', 1)

#set magnetic moment (units of muB)
ns.setparam(meshname, 'mu_s', mu_s)
#set atomistic exchange (units J)
ns.setparam(meshname, 'J', J)
#anisotropy constant (units J) - set uniaxial anisotropy by default, but routines below will change this if needed
ns.addmodule(meshname, 'aniuni')
ns.setparam(meshname, 'K1', K)
#damping set to 1 in case sLLG is used
ns.setparam(meshname, 'damping', 1)

#don't update display often to save CPU power
ns.iterupdate(1000)

######################################################
# Simulate. Uncomment as needed.

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





