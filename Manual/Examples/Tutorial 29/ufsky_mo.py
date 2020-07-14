"""
This script is part of Boris Computational Spintronics v2.8

@author: Serban Lepadatu, 2020
"""

import os
import sys
from NetSocks import NSClient

#setup communication with server. By default sent messages are not displayed in console. 
#To enable verbose mode use e.g.: NSClient('localhost', True)
ns = NSClient('localhost')

########################################

#the working directory : same as this script file, typically expecting simulation file to be in same directory as this script file
directory = os.path.dirname(sys.argv[0]) + "/"
#make sure Boris is reset to default in case a ready set simulation file is not loaded below (if so this can be commented out)
ns.default()
#set working directory same as this script file
ns.chdir(directory)

########################################

#negative helicity
ns.loadsim('ufsky_fm_mo')

ns.setparam('Co/Pt', 'Hmo', -40e6)

for i in range(30):
    ns.reset()
    ns.Run()
    ns.savemeshimage('ufsky_mo_nve_pulse_%d' % (i))
    
#positive helicity
ns.loadsim('ufsky_fm_mo')

ns.setparam('Co/Pt', 'Hmo', +40e6)

for i in range(30):
    ns.reset()
    ns.Run()
    ns.savemeshimage('ufsky_mo_pve_pulse_%d' % (i))

