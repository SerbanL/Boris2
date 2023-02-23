from NetSocks import NSClient, Shape

ns = NSClient(); ns.configure(True)

########################################

l, w, t = 800e-9, 800e-9, 100e-9

FM = ns.Ferromagnet([l, w, t], [5e-9])
ns.delrect()

tor = Shape.torus([l, w, t])
disk = Shape.disk([l/2, w/2, t])

#set disk centred
FM.shape_set(disk.move([l/2, w/2, t/2]))
#set torus centred
FM.shape_set(tor.move([l/2, w/2, t/2]))

#set magnetization angles of defined shapes
FM.shape_setangle(disk, [180, 0])
FM.shape_setangle(tor, [90, 90])

#set material parameter values of regions defined by shapes
Ms = FM.param.Ms.setparam()
FM.param.Ms.shape_setparam(disk, 600e3 / Ms)
FM.param.Ms.shape_setparam(tor, 1200e3 / Ms)




