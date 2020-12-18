"""
This script is part of Boris Computational Spintronics v3.0

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient

#setup communication with server
ns = NSClient()
ns.configure(True)

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

