import os
import sys
from WinSocks import WSClient

#setup communication with server. By default sent messages are not displayed in console. 
#To enable verbose mode use e.g.: WSClient('localhost', True)
ws = WSClient('localhost')

########################################

#the working directory : same as this script file, typically expecting simulation file to be in same directory as this script file
directory = os.path.dirname(sys.argv[0]) + "\\"
ws.chdir(directory)

########################################

#Simulation script here


