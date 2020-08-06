"""
This script is part of Boris Computational Spintronics Article

Generates data in Figure 6(b).

@author: Serban Lepadatu, 2020
"""

import os
import sys
from NetSocks import NSClient
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

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
mpl.rcParams['legend.fontsize'] = 10.5
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
#restore program to default state
ns.default()
ns.chdir(directory)

########################################

dimensions = [320e-9, 160e-9, 10e-9]
contact_extension = 100e-9
contact_foothold = 30e-9
contact_thickness = 5e-9

max_field = 600e3
field_steps = 600

use_prepared_sims = True

def setup_and_simulate_amr_loop(direction_deg):

    ns.default()
    ns.chdir(directory)
    
    #Setup permalloy ellipse
    ns.meshrect(dimensions)
    ns.cellsize([5e-9, 5e-9, 5e-9])
    ns.ecellsize([2.5e-9, 2.5e-9, 2.5e-9])
    
    ns.loadmaskfile(filename = 'Circle')
    
    ns.addmodule('permalloy', 'transport')
    
    #add electrical contacts
    contact1 = ns.addmaterial(
            'Ru', 
            [-contact_extension, 0, dimensions[2], contact_foothold, dimensions[1], dimensions[2] + contact_thickness])
    
    ns.meshfocus2(contact1)
    ns.ecellsize([2.5e-9, 2.5e-9, 2.5e-9])
    
    contact2 = ns.addmaterial(
            'Ru', 
            [dimensions[0] - contact_foothold, 0, dimensions[2], dimensions[0] + contact_extension, dimensions[1], dimensions[2] + contact_thickness])
    
    ns.meshfocus2(contact2)
    ns.ecellsize([2.5e-9, 2.5e-9, 2.5e-9])
    
    ns.addmodule(contact1, 'transport')
    ns.addmodule(contact2, 'transport')
    
    #set electrodes at x-axis ends with a 1 mV potential drop
    ns.setdefaultelectrodes()
    ns.setpotential(1e-3)
    
    #set amr percentage of 2%
    ns.setparam('permalloy', 'amr', 2.0)
    
    ns.setode('LLGStatic', 'RKF45')
    #only need current density each step after magnetisation has relaxed
    #without static transport solver switching events are very slow to simulate as the transport solver will be iterating a lot
    ns.statictransportsolver(1)
    
    #1. Relax transport solver to set convergence error
    ns.setstage('Relax')
    ns.editstagestop(0, 'mxh', 1e-5)
    ns.setfield(max_field, 90, 180.0 + direction_deg)
    ns.setangle(90, 180.0 + direction_deg)
    #need a low connvergence error, which precludes using cuda in single precision
    #you can still use cuda, but need double precision to relax transport solver to low enough convergence error
    #the biggest problem with this type of simulations where current density is not uniform is loop opening
    #(different resistance value for start and end)
    #there are several reasons for this:
    #1. current density not computed accurately enough, so low convergence error helps (this is the main reason)
    #2. simulating a minor loop due to maximum field too low
    #3. first simulation loop can show opening even for large maximum field due to very slightly different starting and ending magnetisation states
    #   thus run simulation once, then save ending state and run the simulation again: should be better now; repeat if needed.
    #That being said, good AMR loops can be obtained just by solving first 2 problems.
    ns.tsolverconfig(2e-9, 10000)
    ns.meshfocus2('permalloy')
    
    #relax transport solver without cuda : better in double floating point precision
    #if you run this in cuda mode single precision is not really good enough as solver won't be able to relax to low enough convergence error 
    #this results in significant loop opening; can still run in cuda but will need the double precision compilation.
    ns.cuda(0)
    ns.Run()
    
    #setup AMR loop stages
    ns.setstage('Hpolar_seq', 'permalloy')
    ns.editstagevalue(0, [-max_field, 90, direction_deg, max_field, 90, direction_deg, field_steps])
    ns.editstagestop(0, 'mxh', 1e-5)
    ns.editdatasave(0, 'step')
    ns.addstage('Hpolar_seq', 'permalloy')
    ns.editstagevalue(1, [max_field, 90, direction_deg, -max_field, 90, direction_deg, field_steps])
    ns.editstagestop(1, 'mxh', 1e-5)
    ns.editdatasave(1, 'step')
    
    #save applied field (A/m) and resistance (Ohms)
    ns.setdata('Ha', 'permalloy')
    ns.adddata('R')
    output_file = 'amr_rawdata_deg_%d.txt' % direction_deg
    ns.savedatafile(output_file)
    
    #get loop
    ns.Run()
    
    #load all columns from file (0, 1, 2, 3) into internal arrays (0, 1, 2, 4)
    ns.dp_load(output_file, [0, 1, 2, 3, 0, 1, 2, 4])
    #get field strength along loop direction and save it in internal array 3
    ns.dp_dotprod(0, np.cos(np.radians(direction_deg)), np.sin(np.radians(direction_deg)), 0, 3)
    #save field strength and resistance in processed file, then plot it here
    loop_file = 'amr_loop_deg_%d.txt' % direction_deg
    ns.dp_save(loop_file, [3, 4])
    amr_data = ns.Get_Data_Columns(loop_file, [0, 1])
    plt.axes(xlabel = 'H (A/m)', ylabel = 'R (Ohms)', title = 'AMR Loop')
    plt.plot(amr_data[0], amr_data[1])
    plt.show()
    
def simulate_amr_loop(simfile, direction_deg):
    
    ns.loadsim(simfile)
    
    output_file = ns.savedatafile()
    
    #get loop
    ns.Run()
    
    #load all columns from file (0, 1, 2, 3) into internal arrays (0, 1, 2, 4)
    ns.dp_load(output_file, [0, 1, 2, 3, 0, 1, 2, 4])
    #get field strength along loop direction and save it in internal array 3
    ns.dp_dotprod(0, np.cos(np.radians(direction_deg)), np.sin(np.radians(direction_deg)), 0, 3)
    #save field strength and resistance in processed file, then plot it here
    loop_file = 'amr_loop_deg_%d.txt' % direction_deg
    ns.dp_save(loop_file, [3, 4])
    amr_data = ns.Get_Data_Columns(loop_file, [0, 1])
    plt.axes(xlabel = 'H (A/m)', ylabel = 'R (Ohms)', title = 'AMR Loop')
    plt.plot(amr_data[0], amr_data[1])
    plt.show()
    
    
########################################

for direction_deg in [1, 31, 61, 91]:
    
    simfile = 'amr_%d' % direction_deg
    
    if not use_prepared_sims:
        setup_and_simulate_amr_loop(direction_deg)
        
        #save simulation : next time you can just run this instead of setting everything up and relaxing transport solver - faster.
        ns.reset()    
        ns.savesim(simfile)
        
    else:
        simulate_amr_loop(simfile, direction_deg)

    
