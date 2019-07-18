import os
from WinSocks import *

#setup communication with server.
ws = WSClient('localhost', True)

#########################################################################

def get_skyrmion_diameter(layers_names):

    sky_diameter = 0.0

    for layer in layers_names:

        #get mesh dimensions
        ws.SendCommand('meshfocus', [layer])
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

        sky_diameter += Get(skyrmion_data, 0) * 2

    #average diameter
    sky_diameter /= len(layers_names)

    return sky_diameter


#########################################################################

#script must be in the same directory as the simulation files
directory = os.path.dirname(sys.argv[0]) + "\\"
ws.SendCommand('chdir', [directory])

#the preconfigured files to simulate
file_names = [
    'Pt_disk_x2',
    'Pt_disk_x3',
    'Pt_disk_x4',
    'Pt_disk_x5',
    'Pt_disk_x6']

#output results here
output_file = 'skyrmion_diameters.txt'

#the layer names in the simulation files
layers_names = [
    ['Co', 'Co2'],
    ['Co', 'Co2', 'Co3'],
    ['Co', 'Co2', 'Co3', 'Co4'],
    ['Co', 'Co2', 'Co3', 'Co4', 'Co5'],
    ['Co', 'Co2', 'Co3', 'Co4', 'Co5', 'Co6']]

#required start fields in A/m for each of the simulation files
start_fields = [4000, 6000, 8000, 10000, 12000]

#field step and last field value field in A/m - applies to all
field_step = 1000
end_field = 20000

################

#do both multiconvolution (1) and supermesh convolution (0)
multiconvolution = [1, 0]

for convtype in multiconvolution:

    sim_idx = 0

    for start_field in start_fields:

        #load simulation
        ws.SendCommand('loadsim', [file_names[sim_idx]])

        #set convolution type
        ws.SendCommand('multiconvolution', [multiconvolution])

        field_steps = (end_field - start_field) / field_step

        #sweep field and calculate skyrmion diameter
        for i in range(0, field_steps + 1):

            field = start_field + i * field_step

            #set field in polar coordinates
            ws.SendCommand('setfield', [field, 0, 0])

            #relax skyrmion to new equilibrium diameter
            ws.Run()

            #now get average skyrmion diameter from all the layers
            diameter = get_skyrmion_diameter(layers_names[sim_idx])

            #output : convolution type - sim file name - field set (A/m) - diameter calculated (m)
            ws.SaveDataToFile(output_file, [convtype, file_names[sim_idx], field, diameter])

        sim_idx+=1

#########################################################################




    


    
    


