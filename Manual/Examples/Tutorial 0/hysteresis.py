"""
This script is part of Boris Computational Spintronics v3.0

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient

ns = NSClient()
ns.configure(True)

########################################

#This is similar to Exercise 2.1, done entirely using a Python script

ns.meshrect([160e-9, 80e-9, 10e-9])
ns.cellsize([5e-9, 5e-9, 10e-9])

#setup two stages to sweep field up and down between -100 kA/m and +100kA/m in 100 steps, slightly off-axis.
#setstage sets a single stage, replacing the default stage
ns.setstage('Hxyz_seq')
ns.editstagevalue(0, [-100e3, 1e3, 0, +100e3, 1e3, 0, 100])
#add new stage
ns.addstage('Hxyz_seq')
ns.editstagevalue(1, [100e3, 1e3, 0, -100e3, 1e3, 0, 100])

#stop each field step using |mxh| < 10^-5 condition
ns.editstagestop(-1, 'mxh', 1e-5)
ns.editdatasave(-1, 'step')

#output data : applied field and average magnetisation
ns.setdata('Ha')
ns.adddata('<M>')
ns.savedatafile('hysteresis.txt')

#solve using LLGStatic equation (damping set to 1 and no precession)
ns.setode('LLGStatic', 'RKF45')

#run program
ns.Run()

#output file has field (x, y, z components) in columns 0, 1, 2, and average magnetisation (x, y, z components) in columns 3, 4, 5
hysteresis_data = ns.Get_Data_Columns('hysteresis.txt', [0, 3])
#plot Mx vs Hx
ns.Plot_Data(hysteresis_data[0], hysteresis_data[1], xlabel = 'H (A/m)', ylabel = 'M (A/m)', title = 'Hysteresis Loop')

