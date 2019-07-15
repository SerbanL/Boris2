from WinSocks import *
import os
import math

#setup communication with server. By default sent messages are not displayed in console. 
#To enable verbose mode use e.g.: WSClient('localhost', True)
ws = WSClient('localhost', True)

#########################################################################

directory = os.path.dirname(sys.argv[0])

simfile = 'circle_staircase'
outfile = 'circle_staircase_data.txt'

ws.SendCommand('chdir', [directory])
ws.SendCommand('loadsim', [simfile])

for angle in range(0, 360):

    ws.SendCommand('setangle', [90, angle])
    ws.SendCommand('editstagevalue', [0, 1e6 * math.cos(math.radians(angle)), 1e6 * math.sin(math.radians(angle)), 0])
    ws.SendCommand('reset')

    ws.Run()

    e_demag = ws.SendCommand('showdata', ['e_demag'])
    e_rough = ws.SendCommand('showdata', ['e_rough'])

    ws.SaveDataToFile(outfile, [angle, e_demag, e_rough])

    

    
    


