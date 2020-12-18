from NetSocks import NSClient
from NetSocks import Shape

ns = NSClient()
ns.configure(True)

########################################

l, w, t = 800e-9, 800e-9, 800e-9
ns.meshrect([l, w, t])
ns.delrect()

disk = Shape.disk([300e-9, 300e-9, 300e-9])
triangle = Shape.triangle([150e-9, 150e-9, 300e-9])
shape = disk - triangle.move([60e-9, 60e-9, 0])

shape.rotate([10, 10, 10])
shape.scale([0.3, 0.3, 0.3])
shape.setrepetitions([4, 4, 4], [200e-9, 200e-9, 200e-9])

ns.shape_set(shape.setposition([100e-9, 100e-9, 100e-9]))





