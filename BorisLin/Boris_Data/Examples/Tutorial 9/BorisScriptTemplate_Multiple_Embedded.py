from NetSocks import NSClient

ns = NSClient(embedded = True)

####################################

#example simulation method
def simulate(H, Ms, angle):
    ns.configure(True)
    
    FM = ns.Ferromagnet([200e-9, 100e-9, 20e-9], [5e-9])
    FM.param.Ms = Ms
    FM.setfield(H, 90, angle)
    FM.addpinneddata('Ha')
    ns.cuda(1)
    ns.Relax(['iter', 10000])

####################################

angle = 45

#pragma parallel for
for (H, Ms) in zip([10e3, 20e3, 30e3, 40e3], [600e3, 700e3, 800e3, 900e3]):
    simulate(H, Ms, angle)
