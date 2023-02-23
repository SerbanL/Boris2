from NetSocks import NSClient

ns = NSClient(); ns.configure(True)

####################################

Py = ns.Ferromagnet([160e-9, 80e-9, 10e-9], [5e-9])

Py.setfield(3e4, 90, 10)
Py.setangle(90, 180)
ns.editstagestop(0, 'time', 10e-9)
ns.Run()
