from NetSocks import NSClient

ns = NSClient(embedded = True); ns.configure(True)

####################################

Py = ns.Ferromagnet([160e-9, 80e-9, 10e-9], [5e-9])

def check_switching(Mx_previous, Mx_previous_save):    
    M = Py.showdata('<M>')
    
    if M[0] * Mx_previous < 0.0:
        print("Switching detected. Stopping simulation.")
        return 0, M[0], Mx_previous_save
    
    #trigger a save when Mx has changed by more than 5 kA/m
    if M[0] - Mx_previous_save > 5e3:
        Mx_previous_save = M[0]
        ns.savedata()
    
    return 1, M[0], Mx_previous_save

Py.setfield(4e4, 90, 0)
Py.setangle(90, 170)
ns.setsavedata('Mswitching.txt', ['time'], ['<M>', Py])
Mx = Py.showdata('<M>')[0]
ns.Run(check_switching, Mx, Mx)
