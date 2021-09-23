from NetSocks import NSClient, customize_plots

ns = NSClient(); ns.configure(True); customize_plots()

import numpy as np
import matplotlib.pyplot as plt

######################################################
# Setup mesh, modules, and parameters - USER EDIT

#number of iterations per Monte Carlo chunk
iters_chunk = 2000

#number of chunks to relax before time-out
iters_chunk_timeout = 10

#threshold normalized magnetization change for accepting state as relax
dm_threshold = 0.01

#find magnetization value for given external conditions:
#1. run iters_chunk iterations and check if new magnetization length along field (averaged over iters_chunk ensembles)
#   is close enough to previous value (dm_threshold)
#2. if accepted then store value and proceed to next field step, else run another iteration up to iters_chunk_timeout
def __run_chunk(M_previous, H_dir):
    
    ns.reset()
    ns.Run()
    
    ns.dp_load('temporary.txt', [0, 1, 2, 0, 1, 2])
    ns.dp_dotprod(0, H_dir, 3)
    
    M = ns.dp_mean(3)[0]
    
    return abs((M - M_previous) / M_previous) < dm_threshold, M, (M + M_previous) / 2

######################################################

def micromagnetic_hysteresis_MMC(
        outputFile,
        meshRect_nm, cellsize_nm,
        optionalCode_start, optionalCode_steps,
        materialName = 'FM',
        Ms = 800e3, A = 13e-12, K1 = 10e3, Kea1 = [1, 0, 0],
        Temperature = 0.0, T_Curie = 680,
        Hmax = 1e5, Hsteps = 100, H_angles_deg = [90, 0], cuda = True, updown_sweeps = False):
    
    """
    Simulate a simple hysteresis loop, saved in outputFile, for given mesh dimensions (meshRect_nm), with given cellsize_nm.
    Use micromagnetic Monte Carlo solver.
    The user can inject custom configuration code through optionalCode_start function before simulation starts, which has the form def optionalCode_start(ns).
    Further custom code can be injected after each field step is solved, with optionalCode_steps, which has the form optionalCode_steps(ns, H)
    Pass in required micromagnetic parameters, unless materialName is given; in this case set material from database.
    The field is swept along H_angles_deg (theta, phi), in given number of steps between -Hmax and Hmax.
    If updown_sweeps is False then only the up sweep is simulated, and the down sweep completed using dp_completehysteresis, else both sweeps are simulated.
    Return: M vs H loop as H, M lists.
    """
    
    ns = NSClient(); ns.configure(True); customize_plots()
    
    ####### Setup
    
    if materialName == 'FM':
        ns.setmesh(materialName, np.array(meshRect_nm) * 1e-9)
    else: 
        ns.setmaterial(materialName, np.array(meshRect_nm) * 1e-9)
        
    ns.cellsize(np.array(cellsize_nm) * 1e-9)
        
    if materialName == 'FM':
        ns.setparam(materialName, 'Ms', Ms)
        ns.setparam(materialName, 'A', A)
        ns.setparam(materialName, 'K1', K1)
        ns.setparam(materialName, 'ea1', Kea1)
        if K1 != 0.0: ns.addmodule(materialName, 'aniuni')

    ns.setstage('MonteCarlo')
    ns.editstagestop(0, 'iter', iters_chunk)
    
    ns.editdatasave(0, 'iter', 1)
    ns.setdata('<M>')
    ns.savedatafile('temporary.txt')
    
    ns.curietemperature(T_Curie)
    ns.temperature(Temperature)
    
    if cuda: ns.cuda(1)
    
    ####### User options
    
    optionalCode_start(ns)
    
    #always start from this direction
    ns.setangle(H_angles_deg[0] + 180, H_angles_deg[1])
    
    ####### Run
    
    ns.dp_newfile(outputFile)
    
    theta = np.radians(H_angles_deg[0]); phi = np.radians(H_angles_deg[1])
    H_dir = [np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]
    Mvec = ns.showdata('<M>')
    M_previous = abs(Mvec[0] * H_dir[0] + Mvec[1] * H_dir[1] + Mvec[2] * H_dir[2])
    
    for H in np.arange(-Hmax, Hmax + 1.0, 2 * Hmax / Hsteps):

        ns.setfield(H, H_angles_deg)
        for chunk_iter in range(iters_chunk_timeout):
            run_next, M_previous, M = __run_chunk(M_previous, H_dir)
            if run_next: break
    
        ns.SaveDataToFile(outputFile, [H, M])
        optionalCode_steps(ns, H)
        
    if updown_sweeps:
        for H in np.arange(Hmax, -(Hmax + 1.0), -2 * Hmax / Hsteps):
    
            ns.setfield(H, H_angles_deg)
            for chunk_iter in range(iters_chunk_timeout):
                run_next, M_previous, M = __run_chunk(M_previous, H_dir)
                if run_next: break
        
            ns.SaveDataToFile(outputFile, [H, M])
            optionalCode_steps(ns, H)
    
    ####### Plot

    theta = np.radians(H_angles_deg[0]); phi = np.radians(H_angles_deg[1])
    H_dir = [np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]
 
    if not updown_sweeps: 
        ns.dp_load(outputFile, [0, 1, 0, 1])
        ns.dp_completehysteresis(0, 1)
        ns.dp_save(outputFile, [0, 1])

    data = ns.Get_Data_Columns(outputFile, [0, 1])
    
    plt.axes(xlabel = 'H (kA/m)', ylabel = 'M (kA/m)')
    plt.plot(np.array(data[0]) / 1e3, np.array(data[1]) / 1e3, '.--', label = 'MMC')
    plt.legend()
    plt.show()
    
    return data[0], data[1]
######################################################