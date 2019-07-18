import time
import socket
import sys

maxLenMessage = 4096
busyWaitms = 250
defaultPort = 1542

def StF(string, retVal = 0):

        #sometimes string to float conversion fails, so safer to use this. Specify return value on fail, default is 0

        try:
            fVal = float(string)
            return fVal

        except:
            return retVal

def StF_NZ(string):
    
    #sometimes string to float conversion fails, so safer to use this. Do not return zero: if zero return 1 instead.

    try:
        fVal = float(string)
        if fVal == 0:
            return 1
        else:
            return fVal

    except:
        return 1

def Get(list, element):

    try:
        return list[element]
    except:
        return 0

###############################################################################################################################

class WSClient:
    """Handles windows socket client communication with Boris server"""

    # Create a TCP/IP socket
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    simRunningPollTimer_ms = 0
    simRunningPollInterval_ms = 500

    def __init__(self, serverip, verbose = False):
        self.verbose = verbose
        server_address = (serverip, defaultPort)
        print('connecting to %s:%s' % server_address)
        self.sock.connect(server_address)

    def WaitForResponse(self):

        data = self.sock.recv(maxLenMessage)
        print('RX : %s' % data)

    def Run(self):

        self.SendCommand("run")
        self.WaitForResponse()

    def SendCommand(self, command, values = None):
        """Send a single command with any given parameters from the values list to console, and return values as a list if any"""

        while True:

            message = command

            if values is not None:
            
                for n in range(len(values)):
                    message += ' ' + str(values[n])

            try:
                # Send data
                if self.verbose == True: 
                    self.sock.sendall('>' + message)
                else: 
                    self.sock.sendall('*' + message)

                print('TX : %s' % message)

            except:
                pass

             # Look for the response
            data = self.sock.recv(maxLenMessage)

            print('RX : %s' % data)

            #note, the returned data always starts with a tab
            fields = data.split('\t')

            if len(fields) == 1:
                break;

            if fields[1] != '!':
                break
            else:
                print('Something went wrong! Trying again. Did you enter the correct format? This is what was received: %s' % fields[2])

        #list of floats where there should be floats instead of strings
        fFields = []

        for n in range(1, len(fields)):
            fFields.extend([StF(fields[n])])

        if len(fFields) == 1:
            return fFields[0]
        else:
            return fFields

    def GetString(self, command):
        """Send command to console and return any information sent back as a string"""

        try:
            # Send data
            if self.verbose == True: 
                self.sock.sendall('>' + command)
            else: 
                self.sock.sendall('*' + command)

            print('TX : %s' % command)

        except:
            pass

         # Look for the response
        string = self.sock.recv(maxLenMessage)

        print('RX : %s' % string)

        #note, the returned string always starts with a tab
        fields = string.split('\t')

        return fields[1]

    def SendCommands(self, commands):
        """send multiple commands but no return values sent back to caller"""

        for n in range(len(commands)):
            self.SendCommand(commands[n])

    def IsSimulationRunning(self, simRunningPollInterval_ms = 500):
        """Poll the isrunning flag: keep sending the isrunning command and check return value. The command is sent every simRunningPollInterval_ms ms at most."""

        self.simRunningPollInterval_ms = simRunningPollInterval_ms

        simRunning = 1

        if abs(time.clock()*1000 - self.simRunningPollTimer_ms) > self.simRunningPollInterval_ms:
            self.simRunningPollTimer_ms = time.clock()*1000
            simRunning = self.SendCommand('isrunning')
        else:
            #this method will typically be used in a while loop, so make sure user won't flood cpu with function calls
            time.sleep(simRunningPollInterval_ms/10000)

        return simRunning

    def PollCommand(self, command, simRunningPollInterval_ms = 500):
        """Poll the given command and return values. The command is sent every simRunningPollInterval_ms ms at most."""

        self.simRunningPollInterval_ms = simRunningPollInterval_ms

        fFields = ['']

        if abs(time.clock()*1000 - self.simRunningPollTimer_ms) > self.simRunningPollInterval_ms:
            self.simRunningPollTimer_ms = time.clock()*1000
            fFields = self.SendCommand(command)
        else:
            #this method will typically be used in a while loop, so make sure user won't flood cpu with function calls
            time.sleep(simRunningPollInterval_ms/10000)

        return fFields
    
    def IsInInterval(self, command, interval, component, simRunningPollInterval_ms = 500):
        """Send command and get returned value. Check if the indicated component number of returned values list is contained in the given interval. The command is sent every simRunningPollInterval_ms ms at most."""

        if len(interval) != 2:
            return 1

        if (interval[0] == '-inf' and interval[1] == '+inf') or interval[0] == '+inf' or interval[1] == '-inf':
            return 1

        fFields = self.PollCommand(command, simRunningPollInterval_ms)

        fValue = 0.0

        if fFields is float:
            fValue = fFields
        elif fFields[0] is not '' and component >= 0 and component < len(fFields):
            fValue = fFields[component]
        else:
            return 1

        if interval[0] == '-inf' and fValue <= interval[1]:
            return 1

        if interval[1] == '+inf' and fValue >= interval[0]:
            return 1

        #not in interval so signal this - default is to assume value is in interval
        return 0

    def SaveDataToFile(self, fileName, dataList):

        command = 'savecomment ' + fileName + ' '

        for n in range(len(dataList)):
            
            command += str(dataList[n])
            if n < len(dataList) - 1:
                command += '\t'

        self.SendCommand(command)
     


