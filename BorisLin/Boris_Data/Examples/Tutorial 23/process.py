"""
This script is part of BORIS

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient, customize_plots
import matplotlib.pyplot as plt
import numpy as np

#setup communication with server
ns = NSClient(); ns.configure(True); customize_plots()

########################################

files = ['PtCo_skyrmion', 'PtCo_skyrmion_noSHE', 'PtCo_skyrmion_SOT']

#obtain polar coordinates movement path for all files above and plot
for file in files:

    ns.dp_load(file, [1, 2, 1, 2])
    ns.dp_replacerepeats(1, 1)
    ns.dp_replacerepeats(2, 2)
    ns.dp_removeoffset(1, 1)
    ns.dp_removeoffset(2, 2)
    ns.dp_cartesiantopolar(1, 2, 3, 4)
    output_file = file + '_rtheta.txt'
    ns.dp_save(output_file, [3, 4])
    
    data = ns.Get_Data_Columns(output_file, [0, 1])
    skypath_r = [r for r in data[0]]
    skypath_theta = [np.radians(theta) for theta in data[1]]
    plt.polar(skypath_theta, skypath_r, label = file)
    
#from first two files obtain the difference path in polar coordinates
#this is spin transport minus diffusion effect path, comparable to that obtained with SOT
ns.dp_load(files[0], [1, 2, 1, 2])
ns.dp_replacerepeats(1, 1)
ns.dp_replacerepeats(2, 2)
ns.dp_removeoffset(1, 1)
ns.dp_removeoffset(2, 2)

ns.dp_load(files[1], [1, 2, 3, 4])
ns.dp_replacerepeats(3, 3)
ns.dp_replacerepeats(4, 4)
ns.dp_removeoffset(3, 3)
ns.dp_removeoffset(4, 4)

ns.dp_subdp(1, 3, 5)
ns.dp_subdp(2, 4, 6)
ns.dp_cartesiantopolar(5, 6, 7, 8)
output_file = files[0] + '_minus_diffusion_rtheta.txt'
ns.dp_save(output_file, [7, 8])

#Show agreement between analytical SOT path and path obtained as the difference between full spin torque and iSTT only (noSHE)
#NOTE : In PtCo_skyrmion_SOT simulation the SHA entered, used by the SOTfield module, is the effective spin Hall angle at the interface, not the intrinsic spin-Hall angle.
#The effective spin Hall angle can be computed using the dp_calcsot command : dp_calcsot Pt Co (see manual for SHAeff formula)
#This will return a value of ~0.03 (intrinsic value is 0.19). When you run the simulation however you must scale it by sigma Pt / sigma Co factor,
#i.e. set a value of 0.03 * 7 / 5 ~= 0.042.
#Reason for this the SOT is generated using the current density in Pt in the self-consistent spin transport solver (as it should be),
#however when you run the SOTfield module the SOT is generated using the current density in Co (SOTfield module is local to the Co mesh so not expected to know about neighboring heavy metal meshes).
#Thus to get the same SOT strength as you do with the spin transport solver, need to account for the different conductivities of Pt and Co set in the simulation.
data = ns.Get_Data_Columns(output_file, [0, 1])
skypath_r = [r for r in data[0]]
skypath_theta = [np.radians(theta) for theta in data[1]]
plt.polar(skypath_theta, skypath_r, '.', label = output_file, zorder = 0)

plt.legend()
plt.savefig('skyrmion.png')
plt.show()



    
    


