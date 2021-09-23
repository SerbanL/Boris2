from NetSocks import NSClient, NSClientConfig

def simulate(Ms, nscfg = NSClientConfig()):

    ns = NSClient(nscfg); ns.configure(True, False)
    
    ns.setparam(ns.meshfocus(), 'Ms', Ms)
    
    ns.cuda(1)
    ns.Run()
