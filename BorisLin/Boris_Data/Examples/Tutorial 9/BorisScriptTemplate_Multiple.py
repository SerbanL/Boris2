from NetSocks import NSMultiClient

nsm = NSMultiClient(scriptserverports = range(1000, 1002), cudaDevices = range(0, 2))
#reset to default, but turn off verbosity so messages do not clash due to parallel execution
nsm.configure(True, False)

####################################

#example simulation method
#this must take an NSClient object as first argument, and any other number of arguments after
def simulate(ns, H, Ms, angle):

    FM = ns.Ferromagnet([200e-9, 100e-9, 20e-9], [5e-9])
    FM.param.Ms = Ms
    FM.setfield(H, 90, angle)
    FM.addpinneddata('Ha')
    ns.cuda(1)
    ns.Relax(['iter', 10000])

####################################

#example parameters passed to simulate function
H = [10e3, 20e3, 30e3, 40e3]
Ms = [600e3, 700e3, 800e3, 900e3]
angle = 45

#execute simulations (4 in total) in parallel using all configured devices (2 devices, thus 2 simulations each)
nsm.Run(simulate, H, Ms, [angle]*len(H))