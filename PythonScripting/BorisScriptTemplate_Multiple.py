from NetSocks import NSClient
import threading

#edit as needed
ports = [1000, 1001]
nsc = [NSClient('localhost', port) for port in ports]
for ns in nsc: ns.configure(True)

########################################

def Simulate(ns, idx):
    print("Simulation number : ", idx)
    #simulation configuration goes here
    ns.Run()
    
########################################
    
sims = [threading.Thread(target = Simulate, args=(ns,idx)) for (ns, idx) in zip(nsc, range(len(nsc)))]
for sim in sims: sim.start()