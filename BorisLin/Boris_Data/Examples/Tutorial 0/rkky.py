"""
This script is part of BORIS
@author: Dr. Serban Lepadatu, 2023

RKKY hysteresis loop example using a synthetic ferrimagnet.
"""

from NetSocks import NSClient
import matplotlib.pyplot as plt
import numpy as np

ns = NSClient(); ns.configure(True)

########################################

direction_deg = 1.0

FM1 = ns.Ferromagnet([320e-9, 160e-9, 10e-9], [5e-9, 5e-9, 5e-9])
FM1.modules(['Zeeman', 'demag', 'exchange', 'surfexchange'])
#shape mesh as an ellipse (mask file is stretched to mesh aspect ratios)
FM1.loadmaskfile('Circle')

#add a new ferromagnetic mesh (permalloy by default) above the first one with 1 nm separation
FM2 = ns.Ferromagnet([0.0, 0.0, 11e-9, 320e-9, 160e-9, 31e-9], [5e-9, 5e-9, 5e-9])
FM2.modules(['Zeeman', 'demag', 'exchange', 'surfexchange'])
FM2.loadmaskfile('Circle')

#enable multilayered demag field calculation allowing exact and efficient calculation of demag fields, even though the separation between meshes is 1 nm
ns.addmodule('supermesh', 'sdemag')

ns.setsavedata('rkky_hysteresis.txt', ['Ha', FM1], ['<M>', FM1], ['<M>', FM2])
ns.setode('LLGStatic', 'SDesc')
ns.cuda(1)
ns.Hpolar_seq(['supermesh', [-300e3, 90, direction_deg, 300e3, 90, direction_deg, 300], 'mxh', 1e-5, 'step'])
ns.Hpolar_seq(['supermesh', [300e3, 90, direction_deg, -300e3, 90, direction_deg, 300], 'mxh', 1e-5, 'step'])

########################################

u = [np.cos(np.radians(direction_deg)), np.sin(np.radians(direction_deg)), 0]
ns.dp_load('rkky_hysteresis', [0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 5, 6, 7, 8])

ns.dp_dotprod(0, u[0], u[1], u[2], 10)
ns.dp_dotprod(3, u[0], u[1], u[2], 11)
ns.dp_dotprod(6, u[0], u[1], u[2], 12)
ns.dp_mul(11, 1.0/3)
ns.dp_mul(12, 2.0/3)
ns.dp_adddp(11, 12, 13)
ns.dp_save('rkky_hysteresis_loop.txt', [10, 13])

loop = ns.Get_Data_Columns('rkky_hysteresis_loop.txt', [0, 1])
plt.axes(xlabel = 'H (kA/m)', ylabel = 'M (kA/m)', title = 'RKKY Hysteresis Loop')
plt.plot(np.array(loop[0])/1e3, np.array(loop[1])/1e3)
plt.show()

