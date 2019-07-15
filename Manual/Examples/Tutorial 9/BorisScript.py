from WinSocks import *

#setup communication with server. By default sent messages are not displayed in console. 
#To enable verbose mode use e.g.: WSClient('localhost', True)
ws = WSClient('localhost')

######################################
#
#fields = ws.SendCommand(command)
#
#This sends the command together with any parameters as a single string to the console, and receives any returned parameters as a list of floats - if only one float returned then this is returned as a stand-alone float
#
#You can also send a list of commands using e.g.:
#
#ws.SendCommands([command1, command2])
#
#In this case no parameters are returned.
#
######################################
#
#fields = ws.SendCommand(command, values)
#
#This sends the command string with a list of parameters (parameters could be numbers - usually the case - and/or strings), and receives any returned parameters as a list of floats - if only one float returned then this is returned as a stand-alone float
#
######################################
#
#ws.Run()
#
#This sends the run command to run the simulation, and waits for it to finish.
#
#Alternatively you can use:
#
#running = ws.IsSimulationRunning()
#
#Call this to poll the simRunning flag: typically used in a while loop. Return value is an integer (0 or 1).
#You don't need to worry about flooding the server with requests as there's a built-in polling interval of 500 ms. To use a custom polling interval:
#
#running = ws.IsSimulationRunning(pollTime_ms)
#
######################################
#
#fields = ws.PollCommand(command)
#
#This is similar to SendCommand but is specifically designed to be used for commands which return parameters and has a built-in polling interval of 500 ms. To use a custom polling interval:
#
#fields = ws.PollCommand(command, pollTime_ms)
#
#NOTE:
#
#fields will contain a list of returned parameters as floats or fields = ['']. The latter is returned if command could not be send (because not enough time has elapsed since last command sent).
#Therefore you need to check for this, e.g. the following loop stops when mxh falls below 0.1:
#
#while True:
#    mxh = ws.PollCommand('showdata mxh')
#    if mxh[0] is not '' and mxh[0] < 0.1:
#        break
#
#These type of values can also be checked using the IsInInterval function:
#
#inInterval = ws.IsInInterval(command, interval, component)
#
#command : the console command to send, which should return values
#
#interval = [val1, val2] is the 2-element list specifying an interval against which the returned parameter is checked.
#Special cases val1 = '-inf' and val2 = '+inf' are allowed.
#If multiple parameters are returned then component specifies which parameter should be checked.
#e.g. the following checks when x component of average magnetisation switches from -ve to +ve:
#
#while ws.IsInInterval('showdata avm', ['-inf', 0], 0):
#   pass
#
######################################

#example: set magnetization to the left, a large field to the right, and wait until magnetization switches.

ws.SendCommands(['default',
                 'mesh 128 128 2',
                 'setangle 90 180',
                 'setfield 1e5 90 0',
                 'run'])

while ws.IsInInterval('showdata <M>', ['-inf', 0], 0):
    pass

#the above is the same as:
#
#while True:
#    avm = ws.PollCommand('showdata', ['avm'])
#    if avm[0] is not '' and avm[0] > 0:
#        break;   

ws.SendCommand('stop')
