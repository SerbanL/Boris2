"""
This script is part of Boris Computational Spintronics Article

Generates data in Figure 5.

@author: Serban Lepadatu, 2020
"""

import os
import sys
from NetSocks import NSClient
import matplotlib.pyplot as plt
import matplotlib as mpl

#setup communication with server
ns = NSClient('localhost')

########################################

#plots customizations; see https://matplotlib.org/3.2.1/tutorials/introductory/customizing.html

#axes
mpl.rcParams['axes.linewidth'] = 2.0
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['font.family'] = 'Arial'

#ticks
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15

mpl.rcParams['xtick.major.size'] = 6.0
mpl.rcParams['xtick.major.width'] = 2.0
mpl.rcParams['xtick.minor.size'] = 6.0
mpl.rcParams['xtick.minor.width'] = 2.0

mpl.rcParams['ytick.major.size'] = 6.0
mpl.rcParams['ytick.major.width'] = 2.0
mpl.rcParams['ytick.minor.size'] = 6.0
mpl.rcParams['ytick.minor.width'] = 2.0

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

mpl.rcParams['xtick.major.pad'] = 5.0
mpl.rcParams['xtick.minor.pad'] = 5.0

#legend
mpl.rcParams['patch.linewidth'] = 2.0
mpl.rcParams['legend.fontsize'] = 12.0
mpl.rcParams['legend.edgecolor'] = 'black'
mpl.rcParams['legend.title_fontsize'] = 15
mpl.rcParams['legend.labelspacing'] = 0.1
mpl.rcParams['legend.borderpad'] = 0.3

#math type
mpl.rcParams['mathtext.default'] = 'regular'

mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['savefig.bbox'] = 'tight'

########################################

#the working directory : same as this script file
directory = os.path.dirname(sys.argv[0]) + "/"
ns.default()
ns.chdir(directory)

########################################

output_file = 'rkky.txt'
savesim_file = 'rkky'

thickness1, spacing, thickness2 = 4.9e-9, 0.6e-9, 2.9e-9
repetitions = 10
planesize, planecell = 300e-9, 3e-9
max_field, steps = 500e3, 500

########################################

#make repetitions as [layer1 / spacing / layer2 / spacing]n
def make_repetition(z_height, rep):
    
    #First layer (first repetition set the material to replace existing mesh, for the others add material)
    if z_height == 0.0: layer1 = ns.setmaterial('Co90Fe10/Ru', [0.0, 0.0, z_height, planesize, planesize, z_height + thickness1])
    else: layer1 = ns.addmaterial('Co90Fe10/Ru', [0.0, 0.0, z_height, planesize, planesize, z_height + thickness1])
        
    ns.meshfocus(layer1)
    ns.cellsize([planecell, planecell, thickness1])
    ns.addmodule(layer1, 'aniuni')
    ns.addmodule(layer1, 'surfexchange')
    #setup polycrystalline structure with easy axis direction varying randomly in crystallites
    ns.setparamvar(layer1, 'ea1', 'vorrot2D', [80, 100, -20, 20, 20e-9, 2 * rep])
    
    #Second layer
    layer2 = ns.addmaterial('Co90Fe10/Ru', [0, 0, z_height + thickness1 + spacing, planesize, planesize, z_height + thickness1 + spacing + thickness2])
    ns.meshfocus(layer2)
    ns.cellsize([planecell, planecell, thickness2])
    ns.addmodule(layer2, 'aniuni')
    ns.addmodule(layer2, 'surfexchange')
    #setup polycrystalline structure with easy axis direction varying randomly in crystallites
    ns.setparamvar(layer2, 'ea1', 'vorrot2D', [80, 100, -20, 20, 20e-9, 2 * rep + 1])
    
    #set output data
    
    #save field in first layer only - applied field is the same in all layers
    if z_height == 0.0: ns.setdata('Ha', layer1)
    
    #save average M in all layers
    ns.adddata('<M>', layer1)
    ns.adddata('<M>', layer2)
    
    return z_height + thickness1 + 2 * spacing + thickness2

########################################

z_height = 0.0

for rep in range(repetitions):
    z_height = make_repetition(z_height, rep)

#multilayered demag with pbc
ns.addmodule('supermesh', 'sdemag')
ns.pbc('supermesh', 'x', 10)
ns.pbc('supermesh', 'y', 10)

#setup field sweep
ns.setstage('Hpolar_seq', 'supermesh')
ns.editstagevalue(0, [-max_field, 90, 0, max_field, 90, 0, steps])
ns.editstagestop(0, 'mxh', 1e-5)
ns.editdatasave(0, 'step')

ns.addstage('Hpolar_seq', 'supermesh')
ns.editstagevalue(1, [max_field, 90, 0, -max_field, 90, 0, steps])
ns.editstagestop(1, 'mxh', 1e-5)
ns.editdatasave(1, 'step')

ns.savedatafile(output_file)

ns.setode('LLGStatic', 'SDesc')

#save simulation file so we don't have to set everything up again next time, can just run this instead
ns.savesim(savesim_file)

ns.cuda(1)
#this will take up to 15 to 30 mins to simulate depending on the GPU
ns.Run()

########################################

Mx_columns = [3 + rep * 3 for rep in range(2*repetitions)]
loop = ns.Get_Data_Columns(output_file, [0] + Mx_columns)
Hrange_k = [H / 1e3 for H in loop[0]]

M_k = []
for cell in range(len(Hrange_k)):

    value = 0.0    
    for i in range(repetitions):        
        M1_val = loop[2*i + 1][cell] * thickness1 / (thickness1 + thickness2)
        M2_val = loop[2*i + 2][cell] * thickness2 / (thickness1 + thickness2)        
        value += M1_val + M2_val
    
    M_k.append(value / (1e3 * repetitions))

plt.axes(xlabel = 'H (kA/m)', ylabel = 'M (kA/m)', title = 'Multilayered SyF Hysteresis Loop')
plt.plot(Hrange_k, M_k)
plt.show()



    
    

    


    



























