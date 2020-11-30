import os
from NetSocks import NSClient

#On Linux specify path to BorisLin in boris_path for auto instance startup; on Windows 'C:/Program Files (x86)/Boris' assumed
ns = NSClient('localhost', serverport = 1542)

########################################

directory = os.path.dirname(os.path.realpath(__file__)) + "/"
ns.default()
ns.chdir(directory)

########################################

ns.Run()



