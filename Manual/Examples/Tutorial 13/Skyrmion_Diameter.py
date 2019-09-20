import os
from WinSocks import *

#setup communication with server. By default sent messages are not displayed in console. 
#To enable verbose mode use e.g.: WSClient('localhost', True)
ws = WSClient('localhost', True)

#the working directory : same as the simulation file
directory = os.path.dirname(sys.argv[0]) + "\\"

#the simulation file
simfile = 'skyrmion_diameter'

#results output file
outputfile = 'skyrmion_diameter.txt'

#field sweep range
start_field = 2e3
end_field = 20e3
field_step = 1e3

ws.SendCommand('loadsim', [directory + simfile])

meshrect = ws.SendCommand('meshrect')
length = Get(meshrect, 3) - Get(meshrect, 0)
width = Get(meshrect, 4) - Get(meshrect, 1)
thickness = Get(meshrect, 5) - Get(meshrect, 2)

steps = int((end_field - start_field) / field_step)

for i in range(0, steps + 1):

    field = start_field + i * field_step

    ws.SendCommand('setfield', [field, 0, 0])

    ws.Run()

    ws.SendCommand('dp_getprofile', [0, width/2, thickness/2, length, width/2, thickness/2, 0])

    skyrmion = ws.SendCommand('dp_fitskyrmion', [0, 3])
    diameter = Get(skyrmion, 0) * 2
    diameter_std = Get(skyrmion, 4) * 2

    if diameter == 0.0:
        break

    ws.SaveDataToFile(outputfile, [field, diameter, diameter_std])

    

    

    

    

    

    


    



























