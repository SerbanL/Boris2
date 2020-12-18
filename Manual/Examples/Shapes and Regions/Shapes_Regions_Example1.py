from NetSocks import NSClient
from NetSocks import Shape

ns = NSClient()
ns.configure(True)

########################################

for loop in range(100):

    l, w, t = 800e-9, 800e-9, 100e-9
    ns.meshrect([l, w, t])
    ns.delrect()
    
    tor = Shape.torus([l, w, t])
    disk = Shape.disk([l/2, w/2, t])
    
    #set disk centred
    ns.shape_set(disk.move([l/2, w/2, t/2]))
    #set torus centred
    ns.shape_set(tor.move([l/2, w/2, t/2]))
    
    #set magnetization angles of defined shapes
    ns.shape_setangle(disk, [180, 0])
    ns.shape_setangle(tor, [90, 90])
    
    #set material parameter values of regions defined by shapes
    ns.setparam(ns.meshfocus(), 'Ms', 1)
    ns.shape_setparam('Ms', disk, 600e3)
    ns.shape_setparam('Ms', tor, 1200e3)




