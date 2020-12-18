"""
This script is part of Boris Computational Spintronics v3.0

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient
import matplotlib.pyplot as plt
import numpy as np

#setup communication with server
ns = NSClient()
ns.configure(True)

########################################

def load_and_prepare_sim(filename, rawdata_file):
    
    ns.loadsim(filename)

    #set 2 stages : the first will achieve steady state motion but not save data, the second will save data
    #Note : addstage simply adds a generic stage; we will need to edit the stop and save conditions, as well as the set stage values
    ns.setstage('V')
    ns.addstage('V')

    #configure the stop conditions for the 2 stages
    ns.editstagestop(0, 'time', 5e-9)
    ns.editstagestop(1, 'time', 5e-9)

    #set a save condition of every 10ps for the second V stage only
    ns.editdatasave(1, 'time', 10e-12)

    ns.savedatafile(rawdata_file)

######################################

#the simulation file
filename = 'cidwm_withOe'

#the final output data will be saved here
outputdata_file = 'dwvelocity_vs_Jc_withOe.txt'

#using this requires the DWMScript_noOe was already run
noOe_outputdata_file = 'dwvelocity_vs_Jc_noOe.txt'

#we'll save temporary data to this file so we can perform linear regression on it
rawdata_file = 'dwvelocity_rawdata'

vstart = -4.57e-3
vend = -45.7e-3
steps = 10

for polarity in range(-1, 2, 2):

    #simulate from low to high current strength for each polarity
    load_and_prepare_sim(filename, rawdata_file)

    for step in range(0, steps + 1):

        #the voltage value to set for this step
        voltage = vstart + (vend-vstart) * step / steps

        voltage = voltage * polarity

        #set the voltage values for the 2 stages
        ns.editstagevalue(0, voltage)
        ns.editstagevalue(1, voltage)

        #make sure to reset before simulating with this voltage value
        ns.reset()

        #wait for the 2 stages to finish
        ns.Run()

        #load time vs dwshift raw data : the simulation file is configured so these are in the first 2 columns
        ns.dp_load(rawdata_file, [0, 1, 0, 1])

        #linear regression on shift vs time data : linregdata will contain in this order : g, g_err, c, c_err
        linregdata = ns.dp_linreg(0,1)
        #we need the gradient (g)
        dwvelocity = linregdata[0]
        #uncertainty
        dwvel_err = linregdata[1]

        #get the current density
        Jc = ns.showdata('<Jc>')
        
        print('Jc (A/m2) = %f, dw velocity (m/s) = %.2f +/- %0.2f' % (Jc[0], dwvelocity, dwvel_err))

        #append new entry to output data : current density (along x) and domain wall velocity
        ns.SaveDataToFile(outputdata_file, [Jc[0], dwvelocity, dwvel_err])
   
data_Oe = ns.Get_Data_Columns(outputdata_file, [0, 1, 2])

#u = |J|*P*muB/(e*Ms*(1+beta^2))
P = ns.setparam('permalloy', 'P')
Ms = ns.setparam('permalloy', 'Ms')
beta = ns.setparam('permalloy', 'beta')
muB = 9.274009994e-24
e = 1.60217662e-19
u_data = [-J*P*muB/(e*Ms*(1 + beta**2)) for J in data_Oe[0]]

plt.axes(xlabel = 'u (m/s)', ylabel = 'DW Velocity (m/s)', title = 'with Oe')
plt.plot(u_data, data_Oe[1], 'o')
plt.show()

umag = []
vdiff = []
err = []

for i in range(int(len(u_data)/2)):
    
    vdiff.append(np.abs(data_Oe[1][i+int(len(u_data)/2)]) - np.abs(data_Oe[1][i]))
    umag.append(np.abs(u_data[i]))
    err.append(np.sqrt(data_Oe[2][i]**2 + data_Oe[2][int(len(u_data)/2)]**2))

plt.axes(xlabel = 'u (m/s)', ylabel = 'DW Velocity asymmetry (m/s)', title = 'Oe field effect')
plt.errorbar(umag, vdiff, yerr = err, fmt = 'o')
plt.show()



