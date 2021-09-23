from NetSocks import NSClient
from NetSocks import Shape

ns = NSClient(); ns.configure(True)

########################################

l, w, t = 800e-9, 800e-9, 200e-9
ns.meshrect([l, w, t])
ns.cellsize([4e-9, 4e-9, 4e-9])
ns.delrect()

#define hollow half-torus
tor1 = Shape.torus([l, w, t])
tor2 = Shape.torus([l - t/2 + t*0.75/2, w - t/2 + t*0.75/2, t*0.75])
htor = tor1 - tor2 - Shape.rect([l, w/2, t]).move([0, -w/4, 0])

#define masks for left and right sides of half-torus
left = htor - Shape.rect([l/2, w/2, t]).move([l/4, w/4, 0])
right = htor - Shape.rect([l/2, w/2, t]).move([-l/4, w/4, 0])

#set shape and magnetization directions for left and right sides
ns.shape_set(htor.move([l/2, w/2, t/2]))
ns.shape_setangle(left.move([l/2, w/2, t/2]), [90, 90])
ns.shape_setangle(right.move([l/2, w/2, t/2]), [90, 270])

#rectangular base
base = Shape.rect([l * 0.6, w*0.5, 10e-9])
ns.shape_set(base.move([l/2, w/4, t/3]))

#repeated tetrahedra
d = 40e-9
tetra = Shape.tetrahedron([d, d, d])
tetra.setrepetitions([1, w / (4*d), 1], [0, d * 2, 0])
ns.shape_set(tetra.move([l * 0.2 + d, d, t/3 + d/2]))
ns.shape_setangle(tetra, [90, 0])

#repeated pyramids
pyramid = Shape.pyramid([d, d, d])
pyramid.setrepetitions([1, w / (4*d), 1], [0, d * 2, 0])
ns.shape_set(pyramid.move([l * 0.2 + 3*d, d, t/3 + d/2]))
ns.shape_setangle(pyramid, [0, 0])

#repeated rotated triangles
triangle = Shape.triangle([d, d, d])
triangle.setrepetitions([1, w / (4*d), 1], [0, d * 2, 0])
triangle.rotate([0, 90, 45])
ns.shape_set(triangle.move([l * 0.2 + 5*d, d, t/3 + d/2]))
ns.shape_setangle(triangle, [180, 0])

#repeated rotated excentric tubes
tube = Shape.disk([d, d, d]) - Shape.disk([d/2, d/2, d]).move([d/5, 0, 0])
tube.setrepetitions([1, w / (4*d), 1], [0, d * 2, 0])
tube.rotate([-30, 45, 45])
ns.shape_set(tube.move([l * 0.2 + 7*d, d, t/3 + d/2]))
ns.shape_setangle(tube, [90, 90])

#repeated ellipsoids
ellipsoid = Shape.ellipsoid([d/2, d*1.5, d*0.75])
ellipsoid.setrepetitions([1, w / (4*d), 1], [0, d * 2, 0])
ns.shape_set(ellipsoid.move([l * 0.2 + 9*d, d, t/3 + d/2]))
ns.shape_setangle(ellipsoid, [90, 270])

#repeated cones
cone = Shape.cone([d, d, d])
cone.setrepetitions([1, w / (4*d), 1], [0, d * 2, 0])
ns.shape_set(cone.move([l * 0.2 + 11*d, d, t/3 + d/2]))
ns.shape_setangle(cone, [45, 45])

ns.displaydetail(2e-9)





