"""
This script is part of BORIS
@author: Dr. Serban Lepadatu, 2023

Exchange bias loop example using an AFM/FM bilayer.
"""

from NetSocks import NSClient
import matplotlib.pyplot as plt
import numpy as np

ns = NSClient(); ns.configure(True)

########################################

AFM = ns.AntiFerromagnet([320e-9, 320e-9, 10e-9], [5e-9, 5e-9, 5e-9])
AFM.modules(['Zeeman', 'exchange', 'aniuni', 'surfexchange'])

#set sub-lattice A magnetisation to result in biasing towards +ve side
AFM.setangle(90, 180)

#Add Fe mesh on top of the antiferromagnet
FM = ns.Material('Fe', [0, 0, 10e-9, 320e-9, 320e-9, 12e-9], [2.5e-9, 2.5e-9, 2e-9])
FM.modules(['Zeeman', 'demag', 'exchange', 'anicubi', 'surfexchange'])
FM.pbc('x', 10)
FM.pbc('y', 10)
#set bilinear surface exchange coupling value - exchange bias is proportional to this
FM.param.J1 = 0.2e-3

########################################

ns.setsavedata('exchangebias.txt', ['Ha', FM], ['<M>', FM])
ns.setode('LLGStatic', 'SDesc')
ns.cuda(1)
ns.Hpolar_seq(['supermesh', [-50e3, 90, 5, 100e3, 90, 5, 50], 'mxh', 1e-5, 'step'])
ns.Hpolar_seq(['supermesh', [100e3, 90, 5, -50e3, 90, 5, 50], 'mxh', 1e-5, 'step'])

########################################

#we should really project along the 5 degree direction, but will keep this simple
data = ns.Get_Data_Columns('exchangebias.txt', [0, 3])
plt.axes(xlabel = 'H (kA/m)', ylabel = 'M (kA/m)', title = 'Exchange bias')
plt.plot(np.array(data[0])/1e3, np.array(data[1])/1e3)
plt.show()


