from NetSocks import NSClientConfig
import threading

#example script
from BorisScriptTemplate import simulate

#edit as needed
nscfglist = [NSClientConfig(scriptserverport = 1000, cudaDevice = 0), NSClientConfig(scriptserverport = 1001, cudaDevice = 1)]

#example parameters passed to simulate function
Mslist = [750e3, 775e3]
    
##############################################################
 
sims = [threading.Thread(target = simulate, args=(Ms, nscfg)) for (Ms, nscfg) in zip(Mslist, nscfglist)]
for sim in sims: sim.start()