import os
from WinSocks import *
ws = WSClient('localhost')

#this script simulates a field-swept FMR peak with bias field along -y direction and r.f. field along x direction


#################################################################################################

#Prepare FMR simulation for given simulation file and ferromagnetic mesh name
#Geometry, parameters and correct modules must be prepared already in the simulation file
#This only changes required stages, output data and data file
def PrepareFMRSimulation(simulationFile_withPath, ferromagnetic_meshName):

    ws.SendCommand('loadsim', [simulationFile_withPath])

    #temporary data file
    ws.SendCommand('savedatafile', ['fmrcycle.txt'])
    
    #delete all stages : only one will be left
    ws.SendCommand('delstage -1')
    #set it to be an fmr stage
    ws.SendCommand('editstage', [0, 'Hfmr', ferromagnetic_meshName])
    
    #add a second fmr stage (this will be the one that save to the temporary output data file)
    ws.SendCommand('addstage', ['Hfmr', ferromagnetic_meshName])
    ws.SendCommand('editdatasave', [1, 'step'])

    #prepare required output data
    #first get rid of everything, leaving just a default time stage : column 0
    ws.SendCommand('deldata -1')

    #applied field : columns 1, 2, 3
    ws.SendCommand('adddata', ['Ha', ferromagnetic_meshName])

    #average magnetisation : columns 4, 5, 6
    ws.SendCommand('adddata', ['<M>', ferromagnetic_meshName])

    #FMR bias field is along -y direction so set starting angle for magnetisation
    ws.SendCommand('setangle', [90, 270])

    #save file
    ws.SendCommand('savesim', [simulationFile_withPath])

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
    threshold = 0.005

    #the stopping time per step required to result in the required fmr frequency (s)
    stopping_time = (1/rfFreq) / 20

    #setup stage values - assume simulation already prepared correctly (e.g. use PrepareFMRSimulation)
    ws.SendCommand('editstagevalue', [0, 0, biasH, 0, rfH, 0, 0, steps_per_cycle, cyclesChunk])
    ws.SendCommand('editstagestop', [0, 'time', stopping_time])

    ws.SendCommand('editstagevalue', [1, 0, biasH, 0, rfH, 0, 0, steps_per_cycle, 1])
    ws.SendCommand('editstagestop', [1, 'time', stopping_time])

    #determine stable oscillation amplitude for the given bias field
    previous_amplitude = 0.0
    cyclesSimulated = 0

    while True:

        ws.SendCommand('reset')
        ws.SendCommand('dp_clearall')

        #run a chunk of cycles to get starting oscillation value
        ws.Run()

        ws.SendCommand('dp_load', ['fmrcycle.txt', 4, 0])
        new_amplitude = ws.SendCommand('dp_getampli', [0, steps_per_cycle])

        cyclesSimulated += cyclesChunk
   
        if new_amplitude:
            
            change = abs( (new_amplitude - previous_amplitude) / new_amplitude )

            if change < threshold or cyclesSimulated >= cyclesMax:
                break

        previous_amplitude = new_amplitude
    
    #save it to file - append bias field and amplitude pair of points
    ws.SaveDataToFile(fileName, [abs(biasH), previous_amplitude])
                         
##################################################################################################

#the directory for the simulation file
import os 
directory = os.path.dirname(sys.argv[0]) + "\\"

#the simulation file correctly prepared for fmr simulations
simulation_file = 'fmr'

output_file = simulation_file + '_fieldsweepFMR_data'

#the meshname in which you want to apply the fmr fields
ferromagnetic_meshName = 'permalloy'

#field-swept fmr : start, stop, step fields (A/m)
Hstart = 320e3
Hend = 420e3
Hstep = 1000

#the FMR rf frequency (Hz)
rfFreq = 20e9

#load simulation and prepare
PrepareFMRSimulation(directory + simulation_file, ferromagnetic_meshName)

#sweep field
for step in range(int((Hend-Hstart)/Hstep) + 1):

    Hbias = Hstart + step * Hstep

    RunFMRStep(-Hbias, rfFreq, ferromagnetic_meshName, output_file)
    



