"""
This script is part of Boris Computational Spintronics Article

Generates data in Figure 7(b).

@author: Serban Lepadatu, 2020
"""

import os
import sys
from NetSocks import NSClient
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.optimize import curve_fit

#setup communication with server
ns = NSClient('localhost')

########################################

#plots customizations; see https://matplotlib.org/3.2.1/tutorials/introductory/customizing.html

#axes
mpl.rcParams['axes.linewidth'] = 2.0
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['font.family'] = 'Arial'

#ticks
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15

mpl.rcParams['xtick.major.size'] = 6.0
mpl.rcParams['xtick.major.width'] = 2.0
mpl.rcParams['xtick.minor.size'] = 6.0
mpl.rcParams['xtick.minor.width'] = 2.0

mpl.rcParams['ytick.major.size'] = 6.0
mpl.rcParams['ytick.major.width'] = 2.0
mpl.rcParams['ytick.minor.size'] = 6.0
mpl.rcParams['ytick.minor.width'] = 2.0

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

mpl.rcParams['xtick.major.pad'] = 5.0
mpl.rcParams['xtick.minor.pad'] = 5.0

#legend
mpl.rcParams['patch.linewidth'] = 2.0
mpl.rcParams['legend.fontsize'] = 14.0
mpl.rcParams['legend.edgecolor'] = 'black'
mpl.rcParams['legend.title_fontsize'] = 15
mpl.rcParams['legend.labelspacing'] = 0.1
mpl.rcParams['legend.borderpad'] = 0.3

#math type
mpl.rcParams['mathtext.default'] = 'regular'

mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['savefig.bbox'] = 'tight'

########################################

#the working directory : same as this script file
directory = os.path.dirname(sys.argv[0]) + "/"
#restore program to default state
ns.default()
ns.chdir(directory)

########################################

#fixed layer polarization direction
p = [-1, 0, 0]

ns.loadsim('cppgmr')

ns.meshfocus('free_layer')
meshrect = ns.meshrect()

#free layer thickness
d = meshrect[5] - meshrect[2]

muB_e = 5.788381608E-05

#vary free layer angle in the plane
phi = 90

########################################

def Tsx(theta, Jc, PL, sigL, PR, sigR, r):
    
    mx = np.cos(np.radians(theta)) * np.sin(np.radians(phi))
    my = np.sin(np.radians(theta)) * np.sin(np.radians(phi))
    mz = np.cos(np.radians(phi))
    m = [mx, my, mz]
    
    qp = (PL * sigL**2 * np.sqrt((sigR**2 + 1)/(sigL**2 + 1)) + PR * sigR**2 * np.sqrt((sigL**2 - 1)/(sigR**2 - 1))) / 2
    qm = (PL * sigL**2 * np.sqrt((sigR**2 + 1)/(sigL**2 + 1)) - PR * sigR**2 * np.sqrt((sigL**2 - 1)/(sigR**2 - 1))) / 2
    A = np.sqrt((sigL**2 + 1)*(sigR**2 + 1))
    B = np.sqrt((sigL**2 - 1)*(sigR**2 - 1))
        
    dp = m[0]*p[0] + m[1]*p[1] + m[2]*p[2]
    neta = qp / (A + B * dp) + qm / (A - B * dp)
    eq_const = muB_e * Jc * neta / d
    m_x_p = [m[1]*p[2]-m[2]*p[1], m[2]*p[0]-m[0]*p[2], m[0]*p[1]-m[1]*p[0]]
    m_x_m_x_p = [m[1]*m_x_p[2] - m[2]*m_x_p[1], m[2]*m_x_p[0] - m[0]*m_x_p[2], m[0]*m_x_p[1] - m[1]*m_x_p[0]]
        
    return eq_const * (m_x_m_x_p[0] + r * m_x_p[0])

########################################

def Tsy(theta, Jc, PL, sigL, PR, sigR, r):
    
    mx = np.cos(np.radians(theta)) * np.sin(np.radians(phi))
    my = np.sin(np.radians(theta)) * np.sin(np.radians(phi))
    mz = np.cos(np.radians(phi))
    m = [mx, my, mz]
    
    qp = (PL * sigL**2 * np.sqrt((sigR**2 + 1)/(sigL**2 + 1)) + PR * sigR**2 * np.sqrt((sigL**2 - 1)/(sigR**2 - 1))) / 2
    qm = (PL * sigL**2 * np.sqrt((sigR**2 + 1)/(sigL**2 + 1)) - PR * sigR**2 * np.sqrt((sigL**2 - 1)/(sigR**2 - 1))) / 2
    A = np.sqrt((sigL**2 + 1)*(sigR**2 + 1))
    B = np.sqrt((sigL**2 - 1)*(sigR**2 - 1))
        
    dp = m[0]*p[0] + m[1]*p[1] + m[2]*p[2]
    neta = qp / (A + B * dp) + qm / (A - B * dp)
    eq_const = muB_e * Jc * neta / d
    m_x_p = [m[1]*p[2]-m[2]*p[1], m[2]*p[0]-m[0]*p[2], m[0]*p[1]-m[1]*p[0]]
    m_x_m_x_p = [m[1]*m_x_p[2] - m[2]*m_x_p[1], m[2]*m_x_p[0] - m[0]*m_x_p[2], m[0]*m_x_p[1] - m[1]*m_x_p[0]]
        
    return eq_const * (m_x_m_x_p[1] + r * m_x_p[1])

########################################

def Tsz(theta, Jc, PL, sigL, PR, sigR, r):
    
    mx = np.cos(np.radians(theta)) * np.sin(np.radians(phi))
    my = np.sin(np.radians(theta)) * np.sin(np.radians(phi))
    mz = np.cos(np.radians(phi))
    m = [mx, my, mz]
    
    qp = (PL * sigL**2 * np.sqrt((sigR**2 + 1)/(sigL**2 + 1)) + PR * sigR**2 * np.sqrt((sigL**2 - 1)/(sigR**2 - 1))) / 2
    qm = (PL * sigL**2 * np.sqrt((sigR**2 + 1)/(sigL**2 + 1)) - PR * sigR**2 * np.sqrt((sigL**2 - 1)/(sigR**2 - 1))) / 2
    A = np.sqrt((sigL**2 + 1)*(sigR**2 + 1))
    B = np.sqrt((sigL**2 - 1)*(sigR**2 - 1))
        
    dp = m[0]*p[0] + m[1]*p[1] + m[2]*p[2]
    neta = qp / (A + B * dp) + qm / (A - B * dp)
    eq_const = muB_e * Jc * neta / d
    m_x_p = [m[1]*p[2]-m[2]*p[1], m[2]*p[0]-m[0]*p[2], m[0]*p[1]-m[1]*p[0]]
    m_x_m_x_p = [m[1]*m_x_p[2] - m[2]*m_x_p[1], m[2]*m_x_p[0] - m[0]*m_x_p[2], m[0]*m_x_p[1] - m[1]*m_x_p[0]]
        
    return eq_const * (m_x_m_x_p[2] + r * m_x_p[2])

########################################

def Ts_mag(theta, Jc, PL, sigL, PR, sigR, r):
    
    return np.sqrt(Tsx(theta, Jc, PL, sigL, PR, sigR, r)*Tsx(theta, Jc, PL, sigL, PR, sigR, r) 
                   + Tsy(theta, Jc, PL, sigL, PR, sigR, r)*Tsy(theta, Jc, PL, sigL, PR, sigR, r)
                   + Tsz(theta, Jc, PL, sigL, PR, sigR, r)*Tsz(theta, Jc, PL, sigL, PR, sigR, r))

########################################

def get_torques(I):
    
    theta_range = []
    Tsi_x_range = []
    Tsi_y_range = []
    Tsi_z_range = []
    Tsi_mag_range = []
    
    Jc_list = []
    
    ns.cuda(0)
    ns.setangle(phi, 0.0, 'free_layer')
    ns.setcurrent(I)
    ns.computefields()
    
    for theta in range(0, 181, 5):
        
        ns.setangle(phi, theta, 'free_layer')
        ns.computefields()
        Tsi = ns.averagemeshrect()
        
        Jc_val = ns.showdata('<Jc>', 'free_layer')
        Jc_list.append(Jc_val[2])
        
        theta_range.append(theta)
        Tsi_x_range.append(Tsi[0])
        Tsi_y_range.append(Tsi[1])
        Tsi_z_range.append(Tsi[2])
        Tsi_mag_range.append(np.sqrt(Tsi[0]*Tsi[0] + Tsi[1]*Tsi[1] + Tsi[2]*Tsi[2]))
        
    Jc = sum(Jc_list) / len(Jc_list)
    
    """
    #This doesn't work well with curve_fit
    eps = 1e-6
    popt, pcov = curve_fit(lambda theta, PL, sigL, PR, sigR, r: Ts_mag(theta, Jc, PL, sigL, PR, sigR, r), theta_range, Tsi_mag_range, 
                           bounds = ((0.0, 1.0 + eps, 0.0, 1.0 + eps, 0.0), (np.inf, np.inf, np.inf, np.inf, np.inf)))
    """
    
    eps = 1e-6
    
    #Fit Tsx
    popt, pcov = curve_fit(lambda theta, PL, sigL, PR, sigR: Tsx(theta, Jc, PL, sigL, PR, sigR, 0.0), theta_range, Tsi_x_range, 
                           bounds = ((0.0, 1.0 + eps, 0.0, 1.0 + eps), (np.inf, np.inf, np.inf, np.inf)))
    
    PLx = popt[0]
    sigLx = popt[1]
    PRx = popt[2]
    sigRx = popt[3]
    
    print('fitting Tsx')
    print('PL = %f, sigL = %f, PR = %f, sigR = %f' % (PLx, sigLx, PRx, sigRx))
    
    #Fit Tsy
    popt, pcov = curve_fit(lambda theta, PL, sigL, PR, sigR: Tsy(theta, Jc, PL, sigL, PR, sigR, 0.0), theta_range, Tsi_y_range, 
                           bounds = ((0.0, 1.0 + eps, 0.0, 1.0 + eps), (np.inf, np.inf, np.inf, np.inf)))
    
    PLy = popt[0]
    sigLy = popt[1]
    PRy = popt[2]
    sigRy = popt[3]
    
    print('fitting Tsy')
    print('PL = %f, sigL = %f, PR = %f, sigR = %f' % (PLy, sigLy, PRy, sigRy))
    
    PL = (PLx + PLy) / 2
    PR = (PRx + PRy) / 2
    sigL = (sigLx + sigLy) / 2
    sigR = (sigRx + sigRy) / 2
    
    qp = (PL * sigL**2 * np.sqrt((sigR**2 + 1)/(sigL**2 + 1)) + PR * sigR**2 * np.sqrt((sigL**2 - 1)/(sigR**2 - 1))) / 2
    qm = (PL * sigL**2 * np.sqrt((sigR**2 + 1)/(sigL**2 + 1)) - PR * sigR**2 * np.sqrt((sigL**2 - 1)/(sigR**2 - 1))) / 2
    A = np.sqrt((sigL**2 + 1)*(sigR**2 + 1))
    B = np.sqrt((sigL**2 - 1)*(sigR**2 - 1))
    
    print('qp = %f, qm = %f, A = %f, B = %f' % (qp, qm, A, B))
    
    #fit r separately from Tsz
    popt, pcov = curve_fit(lambda theta, r: Tsz(theta, Jc, PL, sigL, PR, sigR, r), theta_range, Tsi_z_range)
    r = popt[0]

    print('fitting Tsz')
    print('r = %f' % (r))
        
    Tsi_x_eq_range = [Tsx(theta, Jc, PL, sigL, PR, sigR, r) for theta in theta_range]
    Tsi_y_eq_range = [Tsy(theta, Jc, PL, sigL, PR, sigR, r) for theta in theta_range]
    Tsi_z_eq_range = [Tsz(theta, Jc, PL, sigL, PR, sigR, r) for theta in theta_range]
    
    plt.axes(xlabel = '$\Theta$ (deg.)', ylabel = 'Spin Torque (A/ms)')
    plt.plot(theta_range, Tsi_x_range, 'o', color = 'r', label = 'T$_{s,x}$')
    plt.plot(theta_range, Tsi_y_range, 's', color = 'b', label = 'T$_{s,y}$')
    plt.plot(theta_range, Tsi_z_range, 'D', color = 'g', label = 'T$_{s,z}$')
    plt.plot(theta_range, Tsi_x_eq_range, '--', color = 'black', label = 'Slonczewski')
    plt.plot(theta_range, Tsi_y_eq_range, '--', color = 'black')
    plt.plot(theta_range, Tsi_z_eq_range, '--', color = 'black')
    plt.legend(loc = 'upper right')
    plt.savefig('spin_torque_fit.png')
    plt.show()
    
########################################

ns.meshfocus('free_layer')
ns.display('Tsi', 'free_layer')

get_torques(16.3541e-3)

