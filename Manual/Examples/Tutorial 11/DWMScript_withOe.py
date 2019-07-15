import os
from WinSocks import *

#setup communication with server. By default sent messages are not displayed in console. 
#To enable verbose mode use e.g.: WSClient('localhost', True)
ws = WSClient('localhost')

def load_and_prepare_sim(filename_with_path, rawdata_file):
    
    ws.SendCommand('loadsim', [filename_with_path])

    #set 2 stages : the first will achieve steady state motion but not save data, the second will save data
    #Note : addstage simply adds a generic stage; we will need to edit the stop and save conditions, as well as the set stage values
    ws.SendCommand('addstage V')
    ws.SendCommand('addstage V')

    #delete the first stage : the loaded simulation file has 1 stage set which we don't need
    ws.SendCommand('delstage 0')

    #configure the stop conditions for the 2 stages
    ws.SendCommand('editstagestop', [0, 'time', 4e-9])
    ws.SendCommand('editstagestop', [1, 'time', 4e-9])

    #set a save condition of every 50ps for the second V stage only
    ws.SendCommand('editdatasave', [1, 'time', 50e-12])

    ws.SendCommand('savedatafile', [rawdata_file])

######################################

#the working directory : same as the simulation file
directory = os.path.dirname(sys.argv[0]) + "\\"

#the simulation file
filename = 'cidwm_withOe'

#the final output data will be saved here
outputdata_file = 'dwvelocity_vs_Jc.txt'

#we'll save temporary data to this file so we can perform linear regression on it
rawdata_file = 'dwvelocity_rawdata'

vstart = 2.28e-3
vend = 22.8e-3
steps = 10

for polarity in range(-1, 2, 2):

    #simulate from low to high current strength for each polarity
    load_and_prepare_sim(directory + filename, rawdata_file)

    for step in range(0, steps + 1):

        #the voltage value to set for this step
        voltage = vstart + (vend-vstart) * step / steps

        voltage = voltage * polarity

        #set the voltage values for the 2 stages
        ws.SendCommand('editstagevalue', [0, voltage])
        ws.SendCommand('editstagevalue', [1, voltage])

        #make sure to reset before simulating with this voltage value
        ws.SendCommand('reset')

        #wait for the 2 stages to finish
        ws.Run()

        #load time vs dwshift raw data : the simulation file is configured so these are in the first 2 columns
        ws.SendCommand('dp_load', [rawdata_file, 0, 1, 0, 1])

        #linear regression on shift vs time data : linregdata will contain in this order : g, g_err, c, c_err
        linregdata = ws.SendCommand('dp_linreg', [0,1])
        #we need the gradient (g)
        dwvelocity = Get(linregdata, 0)
        #uncertainty
        dwvel_err = Get(linregdata, 1)

        #get the current density
        Jc = ws.SendCommand('showdata <Jc>')

        #append new entry to output data : current density (along x) and domain wall velocity
        ws.SaveDataToFile(outputdata_file, [Get(Jc, 0), dwvelocity, dwvel_err])

    
    

    


    



























