"""
This script is part of Boris Computational Spintronics Article.

Generates data in Figure 9(a).

@author: Serban Lepadatu, 2020
"""

import os
import sys
from NetSocks import NSClient
import numpy as np
import scipy as sp
import scipy.fftpack
import matplotlib.pyplot as plt
import matplotlib as mpl

#setup communication with server. By default sent messages are not displayed in console. 
#To enable verbose mode use e.g.: NSClient('localhost', True)
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
mpl.rcParams['legend.fontsize'] = 12.0
mpl.rcParams['legend.edgecolor'] = 'black'
mpl.rcParams['legend.title_fontsize'] = 15
mpl.rcParams['legend.labelspacing'] = 0.1
mpl.rcParams['legend.borderpad'] = 0.3

#math type
mpl.rcParams['mathtext.default'] = 'regular'

mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['savefig.bbox'] = 'tight'

########################################

#the working directory : same as this script file, typically expecting simulation file to be in same directory as this script file
directory = os.path.dirname(sys.argv[0]) + "/"
#make sure Boris is reset to default in case a ready set simulation file is not loaded below (if so this can be commented out)
ns.default()
#set working directory same as this script file
ns.chdir(directory)

########################################

def func_lorentz(x, y0, S, A, w, x0):
    return y0 + S * (w + A*(x - x0)) / (4*(x-x0)**2 + w**2)

########################################  

def setup_af_fmr_simulation():
    
    meshdim = [80e-9, 80e-9, 2e-9]
    cellsize = [5e-9, 5e-9, 2e-9]
    
    output_file = 'afmr_peak.txt'
    
    #set thin-film Co/Pt with demag, exchange, Zeeman (enabled by default), and uniaxial anisotropy
    ns.setafmesh('AF', meshdim)
    ns.cellsize(cellsize)
    ns.addmodule('AF', 'aniuni')
    
    ns.setparam('AF', 'ea1', [0,1,0])
    
    ns.setparam('AF', 'Anh', [0.0, 0.0])
    
    #thin film geometry : set pbc conditions
    ns.addmodule('AF', 'demag')
    ns.pbc('AF', 'x', 10)
    ns.pbc('AF', 'y', 10)
    ns.setdtspeedup(1e-9)
    
    #setup sinc pulse using a formula
    ns.setstage('Hequation')
    
    ns.editstagevalue(0, 'He * sinc(2*PI*fc*(t-t0)), H0, 0')

    ns.editstagestop(0, 'time', total_time)
    ns.editdatasave(0, 'time', time_step)
    
    #define the equation constants
    #Excitation field (A/m)
    ns.equationconstants('He', He)
    #f cutoff (Hz)
    ns.equationconstants('fc', fc)
    #sinc pulse center (s)
    ns.equationconstants('t0', total_time / 2)
    
    #save magnetisation only
    ns.setdata('<M>', 'AF')
    ns.savedatafile(output_file)
    
    ns.iterupdate(0)
    
    #RK4 is good for this type of work, just set a low enough time-step
    ns.setode('LLG', 'RK4')
    #make sure the time-step divides the sampling time
    ns.setdt(250e-15)
    
    ns.cuda(1)

########################################

#H0: bias field strength (A/m)
#fc: frequency cutoff (Hz)
def simulate_fmr_peak(H0, show_result = False):
    
    ns.setangle(90, 90)
    
    #Bias field (A/m)
    ns.equationconstants('H0', H0)
    
    ns.reset()
    ns.Run()
    
    #Analyse
    
    #get 2D list as position along horizontal, time along vertical
    Mav = ns.Get_Data_Columns(ns.savedatafile())
    Ms_AFM = ns.setparam('AF', 'Ms_AFM')
    Mx = [row[0]/Ms_AFM[0] for row in Mav]
    
    #2D fft
    fourier_data = np.fft.fftshift(np.abs(sp.fftpack.fft(Mx))**2)
    
    #get value ranges
    freq_len = len(fourier_data)
    freq = sp.fftpack.fftfreq(freq_len, time_step)
    
    #extract result from fourier data in a plottable form
    result = fourier_data[int(freq_len/2):freq_len]
    
    #use Lorentz fitting procedure built into Boris - curve_fit from scipy.optimize fails miserably even if you set bounds or starting guess; 
    #routine in Boris simply works! (and doesn't need any bounds, and it determines starting guesses by itself)
    xy = []
    xy.append(freq[0:len(result)])
    xy.append(result)
    ns.Save_Data_Columns('fmr_temp.txt', xy)
    ns.dp_load('fmr_temp.txt', [0, 1, 0, 1])
    #parameters returned as S, A, x0, w, y0, sigS, siaA, sigx0, sigw, sigy0
    lorentz_params = ns.dp_fitlorentz2(0, 1)
    lorentz = [func_lorentz(f, lorentz_params[4], lorentz_params[0], lorentz_params[1], lorentz_params[3], lorentz_params[2]) for f in freq[0:len(result)]]
    
    if show_result:
    
        print('f0 (GHz) = %.2f, FWHM (GHz) = %.2f' % (lorentz_params[2] / 1e9, lorentz_params[3] / 1e9))    
        plt.axes(xlabel = 'f (Hz)', ylabel = 'FMR Signal (a.u.)', title = 'Frequency-swept FMR')
        plt.plot(freq[0:len(result)], result, '.')
        plt.plot(freq[0:len(result)], lorentz)
        plt.show()
        
    return lorentz_params[2], lorentz_params[3]

########################################    

#Setup    

#cutoff frequency (Hz)
fc = 500e9
#time step for saving magnetisation (s) : determined by Nyquist criterion from cutoff frequency
time_step = (0.5 / fc)
#total time to simulate; increasing this gives you more frequency points in the transform
total_time = 10000e-12
#Field and excitation field strength
He = 1000

setup_af_fmr_simulation()

#useful constants
#saturation magentisation (A/m)
Ms_AFM = ns.setparam('AF', 'Ms_AFM')
#gyromagnetic ratio modulus (times mu0)
grel_AFM = ns.setparam('AF', 'grel_AFM')
gamma = grel_AFM[0] * 2.212761569e5
#Gilbert damping
damping_AFM = ns.setparam('AF', 'damping_AFM')
mu0 = 4*np.pi*1e-7

########################################

#Simulate

H0 = 0.0

#antiferromagnetic exchange (homogeneous)
Ah_range = np.arange(-1e6, -32e6, -2.5e6)

#uniaxial anisotropy
K1_AFM_range = [50e3, 100e3, 150e3, 200e3]

markers = ['o', 's', '^', 'D']

plt.axes(xlabel = 'A$_h$ (J/m$^3$)', ylabel = 'f (GHz)')

for idx in range(len(K1_AFM_range)):
    
    K1_AFM = K1_AFM_range[idx]
    ns.setparam('AF', 'K1_AFM', [K1_AFM, K1_AFM])

    f0_sim_range = []
    f0_Kittel_range = []
    
    for Ah in Ah_range:
    
        ns.setparam('AF', 'Ah', [Ah, Ah])
        
        f0, df = simulate_fmr_peak(H0)
        
        H_Weiss_exchange = 4 * np.abs(Ah) / (mu0 * Ms_AFM[0])
        H_anisotropy = 2 * K1_AFM / (mu0 * Ms_AFM[0])
        
        f0_Kittel = gamma * (H0 + np.sqrt(H_anisotropy * (2 * H_Weiss_exchange + H_anisotropy))) / (2 * np.pi)
        print('Kittel formula, f0 (GHz) = %.2f' % (f0_Kittel / 1e9))
        
        f0_sim_range.append(f0/1e9)
        f0_Kittel_range.append(f0_Kittel/1e9)
        
    plt.plot(Ah_range, f0_sim_range, markers[idx], label = '%d' % (K1_AFM / 1e3))
    if idx == len(K1_AFM_range) - 1: plt.plot(Ah_range, f0_Kittel_range, '--', color = 'black', label = 'Kittel')
    else: plt.plot(Ah_range, f0_Kittel_range, '--', color = 'black')

########################################
    
plt.legend(title = 'K$_1$ (kJ/m$^3$)')
plt.savefig('AFMResonance.png')
plt.show()
    

    

    


