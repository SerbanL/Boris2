from WinSocks import *
import os

#setup communication with server. By default sent messages are not displayed in console. 
#To enable verbose mode use e.g.: WSClient('localhost', True)
ws = WSClient('localhost', True)

#########################################################################

directory = os.path.dirname(sys.argv[0])

files = ['PtCo_skyrmion', 'PtCo_skyrmion_noSHE', 'PtCo_skyrmion_SOT']

ws.SendCommand('chdir', [directory])

#obtain polar coordinates movement path for all files above
for file in files:

    ws.SendCommand('dp_load', [file, 1, 2, 1, 2])
    ws.SendCommand('dp_replacerepeats', [1, 1])
    ws.SendCommand('dp_replacerepeats', [2, 2])
    ws.SendCommand('dp_removeoffset', [1, 1])
    ws.SendCommand('dp_removeoffset', [2, 2])
    ws.SendCommand('dp_cartesiantopolar', [1, 2, 3, 4])
    ws.SendCommand('dp_save', [file + '_rtheta', 3, 4])

#from first two files obtain the difference path in polar coordinates
#this is spin transport minus diffusion effect path, comparable to that obtained with SOT
ws.SendCommand('dp_load', [files[0], 1, 2, 1, 2])
ws.SendCommand('dp_replacerepeats', [1, 1])
ws.SendCommand('dp_replacerepeats', [2, 2])
ws.SendCommand('dp_removeoffset', [1, 1])
ws.SendCommand('dp_removeoffset', [2, 2])

ws.SendCommand('dp_load', [files[1], 1, 2, 3, 4])
ws.SendCommand('dp_replacerepeats', [3, 3])
ws.SendCommand('dp_replacerepeats', [4, 4])
ws.SendCommand('dp_removeoffset', [3, 3])
ws.SendCommand('dp_removeoffset', [4, 4])

ws.SendCommand('dp_subdp', [1, 3, 5])
ws.SendCommand('dp_subdp', [2, 4, 6])
ws.SendCommand('dp_cartesiantopolar', [5, 6, 7, 8])
ws.SendCommand('dp_save', [files[0] + '_minus_diffusion_rtheta', 7, 8])



    
    


