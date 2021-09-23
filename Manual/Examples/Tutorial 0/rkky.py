"""
This script is part of Boris Computational Spintronics

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient
import matplotlib.pyplot as plt
import numpy as np

ns = NSClient(); ns.configure(True)

########################################

direction_deg = 1.0

ns.meshrect([320e-9, 160e-9, 10e-9])
#shape mesh as an ellipse (mask file is stretched to mesh aspect ratios)
ns.loadmaskfile('Circle')
#add a new ferromagnetic mesh (permalloy by default) above the first one with 1 nm separation
#the 2 meshes still retain a cubic 5 nm cellsize
ns.addmesh('permalloy2', [0.0, 0.0, 11e-9, 320e-9, 160e-9, 31e-9])
ns.meshfocus('permalloy2')
ns.loadmaskfile('Circle')
#enable multilayered demag field calculation allowing exact and efficient calculation of demag fields, even though the separation between meshes is 1 nm
ns.addmodule('supermesh', 'sdemag')
#enable RKKY coupling (surface exchange coupling) keeping default J1 (bilinear) and J2 (biquadratic) values
ns.addmodule('permalloy', 'surfexchange')
ns.addmodule('permalloy2', 'surfexchange')
#set field sequence to apply to both meshes (so set it to the supermesh)
ns.setstage('Hpolar_seq', 'supermesh')
ns.editstagevalue(0, [-300e3, 90, direction_deg, 300e3, 90, direction_deg, 300])
ns.editstagestop(0, 'mxh', 1e-5)
ns.editdatasave(0, 'step')
ns.addstage('Hpolar_seq', 'supermesh')
ns.editstagevalue(1, [300e3, 90, direction_deg, -300e3, 90, direction_deg, 300])
ns.editstagestop(1, 'mxh', 1e-5)
ns.editdatasave(1, 'step')

ns.setdata('Ha')
ns.adddata('<M>', 'permalloy')
ns.adddata('<M>', 'permalloy2')
ns.savedatafile('rkky_hysteresis.txt')

ns.setode('LLGStatic', 'SDesc')
ns.cuda(1)

ns.Run()

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
plt.axes(xlabel = 'H (A/m)', ylabel = 'M (A/m)', title = 'RKKY Hysteresis Loop')
plt.plot(loop[0], loop[1])
plt.show()

