"""
This script is part of Boris Computational Spintronics v3.0

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.fftpack
import scipy.signal

#setup communication with server
ns = NSClient()
ns.configure(True)

########################################

def func_lorentz(x, y0, S, A, w, x0):
    return y0 + S * (w + A*(x - x0)) / (4*(x-x0)**2 + w**2)

########################################  

def setup_af_fmr_simulation():
    
    meshdim = [80e-9, 80e-9, 2e-9]
    cellsize = [5e-9, 5e-9, 2e-9]
    
    output_file = 'af_fmr_peak.txt'
    
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
    
    #ns.cuda(1)

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
#anisotropy constant (for Co/Pt the easy axis is in the z direction)
K1_AFM = ns.setparam('AF', 'K1_AFM')
#gyromagnetic ratio modulus (times mu0)
grel_AFM = ns.setparam('AF', 'grel_AFM')
gamma = grel_AFM[0] * 2.212761569e5
#Gilbert damping
damping_AFM = ns.setparam('AF', 'damping_AFM')
mu0 = 4*np.pi*1e-7

Ah = ns.setparam('AF', 'Ah')

########################################

#Simulate

H0 = 0.0

f0, df = simulate_fmr_peak(H0, True)

H_Weiss_exchange = 4 * np.abs(Ah[0]) / (mu0 * Ms_AFM[0])
H_anisotropy = 2 * K1_AFM[0] / (mu0 * Ms_AFM[0])

f0_Kittel = gamma * (H0 + np.sqrt(H_anisotropy * (2 * H_Weiss_exchange + H_anisotropy))) / (2 * np.pi)
print('Kittel formula, f0 (GHz) = %.2f' % (f0_Kittel / 1e9))

########################################
    

    

    


