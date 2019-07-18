import os
import math
from WinSocks import *
ws = WSClient('localhost')

#################################################################################################

#the directory for the simulation file
directory = os.path.dirname(sys.argv[0]) + "\\"

#the output file
data_file = directory + 'trilayer_switch_mconv'

overall_weight = 50.0 / 52.0

weight_layer_1 = 0.4 * overall_weight
weight_layer_2 = 0.2 * overall_weight
weight_layer_3 = 0.4 * overall_weight
Ms = 8e5

ws.SendCommand('dp_load', [data_file, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

ws.SendCommand('dp_mul', [1, weight_layer_1 / Ms])
ws.SendCommand('dp_mul', [2, weight_layer_1 / Ms])
ws.SendCommand('dp_mul', [3, weight_layer_1 / Ms])

ws.SendCommand('dp_mul', [4, weight_layer_2 / Ms])
ws.SendCommand('dp_mul', [5, weight_layer_2 / Ms])
ws.SendCommand('dp_mul', [6, weight_layer_2 / Ms])

ws.SendCommand('dp_mul', [7, weight_layer_3 / Ms])
ws.SendCommand('dp_mul', [8, weight_layer_3 / Ms])
ws.SendCommand('dp_mul', [9, weight_layer_3 / Ms])

ws.SendCommand('dp_adddp', [1, 4, 10])
ws.SendCommand('dp_adddp', [2, 5, 11])
ws.SendCommand('dp_adddp', [3, 6, 12])
ws.SendCommand('dp_adddp', [7, 10, 13])
ws.SendCommand('dp_adddp', [8, 11, 14])
ws.SendCommand('dp_adddp', [9, 12, 15])

ws.SendCommand('dp_save', [data_file + '_processed', 0, 13, 14, 15])

ws.SendCommand('dp_load', [data_file, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

ws.SendCommand('dp_div', [1, Ms])
ws.SendCommand('dp_div', [2, Ms])
ws.SendCommand('dp_div', [3, Ms])

ws.SendCommand('dp_div', [4, Ms])
ws.SendCommand('dp_div', [5, Ms])
ws.SendCommand('dp_div', [6, Ms])

ws.SendCommand('dp_div', [7, Ms])
ws.SendCommand('dp_div', [8, Ms])
ws.SendCommand('dp_div', [9, Ms])

ws.SendCommand('dp_save', [data_file + '_bottom', 0, 1, 2, 3])
ws.SendCommand('dp_save', [data_file + '_middle', 0, 4, 5, 6])
ws.SendCommand('dp_save', [data_file + '_top', 0, 7, 8, 9])



    



