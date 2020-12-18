"""
This script is part of Boris Computational Spintronics v3.0

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient
import matplotlib.pyplot as plt

#setup communication with server
ns = NSClient()
ns.configure(True)

########################################

def func_lorentz(x, y0, S, w, x0):
    return y0 + S * w / (4*(x-x0)**2 + w**2)

########################################

#this script simulates a field-swept FMR peak with bias field along -y direction and r.f. field along x direction

#################################################################################################

#Prepare FMR simulation for given simulation file and ferromagnetic mesh name
#Geometry, parameters and correct modules must be prepared already in the simulation file
#This only changes required stages, output data and data file
def PrepareFMRSimulation(simulationFile, ferromagnetic_meshName):

    ns.loadsim(simulationFile)

    #temporary data file
    ns.savedatafile('fmrcycle.txt')
    
    #set an fmr stage
    ns.setstage('Hfmr', ferromagnetic_meshName)
    
    #add a second fmr stage (this will be the one that save to the temporary output data file)
    ns.addstage('Hfmr', ferromagnetic_meshName)
    ns.editdatasave(1, 'step')

    #save data - applied field : columns 1, 2, 3
    ns.setdata('Ha', ferromagnetic_meshName)

    #save data - average magnetisation : columns 4, 5, 6
    ns.adddata('<M>', ferromagnetic_meshName)

    #FMR bias field is along -y direction so set starting angle for magnetisation
    ns.setangle(90, 270)

    #save file
    ns.savesim(simulationFile)

#################################################################################################

#run a single FMR step for given bias field strength (along y axis in A/m), at rf frequency in Hz
#apply the fields for the named ferromagnetic mesh
#Save fmr data in the given file
#Before calling this the simulation file should be loaded and correctly prepared (PrepareFMRSimulation helps)
def RunFMRStep(biasH, rfFreq, ferromagnetic_meshName, fileName):

    #maximum number of cycles per bias field step
    cyclesMax = 200

    #number of cycles per chunk
    cyclesChunk = 20

    #rf field amplitude (A/m)
    rfH = 100

    #steps per cycle to use
    steps_per_cycle = 20

    #threshold for accepting amplitude
    threshold = 0.001

    #the stopping time per step required to result in the required fmr frequency (s)
    stopping_time = (1/rfFreq) / 20

    #setup stage values - assume simulation already prepared correctly (e.g. use PrepareFMRSimulation)
    ns.editstagevalue(0, [0, biasH, 0, rfH, 0, 0, steps_per_cycle, cyclesChunk])
    ns.editstagestop(0, 'time', stopping_time)

    ns.editstagevalue(1, [0, biasH, 0, rfH, 0, 0, steps_per_cycle, 1])
    ns.editstagestop(1, 'time', stopping_time)

    #determine stable oscillation amplitude for the given bias field
    previous_amplitude = 0.0
    cyclesSimulated = 0

    while cyclesSimulated < cyclesMax:

        ns.reset()
        ns.dp_clearall()

        #run a chunk of cycles to get starting oscillation value
        ns.Run()

        ns.dp_load('fmrcycle.txt', [3, 0])
        new_amplitude = ns.dp_getampli(0, steps_per_cycle)

        cyclesSimulated += cyclesChunk
   
        if new_amplitude:
            
            change = abs( (new_amplitude - previous_amplitude) / new_amplitude )

            if change < threshold:
                break

        previous_amplitude = new_amplitude
    
    #save it to file - append bias field and amplitude pair of points
    ns.SaveDataToFile(fileName, [abs(biasH), new_amplitude])
                         
##################################################################################################

#the simulation file correctly prepared for fmr simulations
simulation_file = 'fmr_she'

output_file = simulation_file + '_fieldsweepFMRshe_data.txt'

#the meshname in which you want to apply the fmr fields
ferromagnetic_meshName = 'permalloy'

#field-swept fmr : start, stop, step fields (A/m)
Hstart = 320e3
Hend = 420e3
Hstep = 1000

#the FMR rf frequency (Hz)
rfFreq = 20e9

#load simulation and prepare
PrepareFMRSimulation(simulation_file, ferromagnetic_meshName)

#sweep field
for step in range(int((Hend-Hstart)/Hstep) + 1):

    Hbias = Hstart + step * Hstep
    print("Bias field = ", Hbias)
    RunFMRStep(-Hbias, rfFreq, ferromagnetic_meshName, output_file)
    
##################################################################################################

data = ns.Get_Data_Columns(output_file, [0, 1])
Hrange = [H for H in data[0]]
fmr_signal = [amplitude**2 for amplitude in data[1]]

ns.dp_load(output_file, [0, 1, 0, 1])
ns.dp_muldp(1, 1, 1)
lorentz_params = ns.dp_fitlorentz(0, 1)
lorentz = [func_lorentz(H, lorentz_params[3], lorentz_params[0], lorentz_params[2], lorentz_params[1]) for H in Hrange]

data = ns.Get_Data_Columns(output_file, [0, 1])
fmr_signal = [amplitude**2 for amplitude in data[1]]

plt.axes(xlabel = 'H (A/m)', ylabel = 'FMR Signal (a.u.)', title = 'Field-swept FMR - SOT from SHE')
plt.plot(Hrange, fmr_signal, 'o', label = 'simulated FMR')
plt.plot(Hrange, lorentz, '--', label = 'Lorentz fit')
plt.legend()
plt.show()