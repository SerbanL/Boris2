from NetSocks import NSClient, customize_plots

import numpy as np
import matplotlib.pyplot as plt

######################################################

def micromagnetic_hysteresis_SDesc(
        outputFile,
        meshRect_nm, cellsize_nm,
        optionalCode_start,
        materialName = 'FM',
        Ms = 800e3, A = 13e-12, K1 = 10e3, Kea1 = [1, 0, 0],
        Temperature = 0.0, T_Curie = 680,
        Hmax = 1e5, Hsteps = 100, H_angles_deg = [90, 0], mxh = 1e-5, cuda = True, updown_sweeps = False):
    
    """
    Simulate a simple hysteresis loop, saved in outputFile, for given mesh dimensions (meshRect_nm), with given cellsize_nm.
    Use SDesc solver.
    The user can inject custom configuration code through optionalCode_start function before simulation starts, which has the form def optionalCode_start(ns).
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
    
    ns.setstage('Hpolar_seq')
    ns.editstagevalue(0, [-Hmax, H_angles_deg[0], H_angles_deg[1], Hmax, H_angles_deg[0], H_angles_deg[1], Hsteps])
    if updown_sweeps:
        ns.addstage('Hpolar_seq')
        ns.editstagevalue(1, [Hmax, H_angles_deg[0], H_angles_deg[1], -Hmax, H_angles_deg[0], H_angles_deg[1], Hsteps])
        
    ns.editstagestop(-1, 'mxh_iter', [mxh, 10000])

    ns.editdatasave(-1, 'step', 1)
    ns.setdata('Ha')
    ns.adddata('<M>')
    ns.savedatafile('rawdata_' + outputFile)
    
    ns.setode('LLGStatic', 'SDesc')
    
    ns.curietemperature(T_Curie)
    ns.temperature(Temperature)
    
    if cuda: ns.cuda(1)
    
    ####### User options
    
    optionalCode_start(ns)
    
    #always start from this direction
    ns.setangle(H_angles_deg[0] + 180, H_angles_deg[1])
    
    ####### Run
    
    ns.Run()
    
    ####### Plot

    theta = np.radians(H_angles_deg[0]); phi = np.radians(H_angles_deg[1])
    H_dir = [np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]

    ns.dp_load('rawdata_' + outputFile, [0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5])
    ns.dp_dotprod(0, H_dir, 6)
    ns.dp_dotprod(3, H_dir, 7)
    
    if not updown_sweeps: ns.dp_completehysteresis(6, 7)
    
    ns.dp_save(outputFile, [6, 7])

    data = ns.Get_Data_Columns(outputFile, [0, 1])
    
    plt.axes(xlabel = 'H (kA/m)', ylabel = 'M (kA/m)')
    plt.plot(np.array(data[0]) / 1e3, np.array(data[1]) / 1e3, '.--', label = 'SDesc')
    plt.legend()
    plt.show()
    
    return data[0], data[1]
    
######################################################