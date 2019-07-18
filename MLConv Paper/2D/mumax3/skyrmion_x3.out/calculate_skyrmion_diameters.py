import os
from WinSocks import *

#setup communication with server.
ws = WSClient('localhost', True)

#########################################################################
# Load a set of OVF files computed by mumax3 into Boris,
# then extract skyrmion diameters using in-built fitting function
#########################################################################

def get_skyrmion_diameter():

    #get mesh dimensions
    meshRect = ws.SendCommand('meshrect')
    x_start = Get(meshRect, 0)
    y_start = Get(meshRect, 1)
    z_start = Get(meshRect, 2)
    x_end = Get(meshRect, 3)
    y_end = Get(meshRect, 4)
    z_end = Get(meshRect, 5)
    length = x_end - x_start
    width = y_end - y_start
    thickness = z_end - z_start

    #skyrmion profile through the center
    ws.SendCommand('dp_getprofile', [x_start, y_start + width/2, z_start + thickness/2, x_end, y_start + width/2, z_start + thickness/2, 0])
    #center the profile about zero
    ws.SendCommand('dp_sub', [0, length/2])
    #LMA algorithm fit of skyrmion profile to extract radius (the fitting uncertainty is typically negligible)
    skyrmion_data = ws.SendCommand('dp_fitskyrmion', [0, 3])

    sky_diameter = Get(skyrmion_data, 0) * 2

    return sky_diameter


#########################################################################

#script must be in the same directory as the simulation files
directory = os.path.dirname(sys.argv[0]) + "\\"

#normalize OVF data to this value
Ms = 6e5

#output results here
output_file = 'skyrmion_diameters.txt'

num_files = 15
layers = [0, 4, 8]

#reset Boris and set directory
ws.SendCommand('default')
ws.SendCommand('chdir', [directory])

for i in range(0, num_files):

    diameter = 0.0

    for l in layers:

        file_name = 'm_zrange'

        file_name += str(l) + '_' + str(i).zfill(6)

        ws.SendCommand('loadovf2mag', [Ms, directory + file_name])

        diameter += get_skyrmion_diameter()

    #average diameter
    diameter /= len(layers)

    ws.SaveDataToFile(output_file, [i, diameter])

        
    

#########################################################################




    


    
    


