"""
This script is part of Boris Computational Spintronics v3.0

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import scipy.fftpack

#setup communication with server
ns = NSClient()
ns.configure(True)

########################################

def func_lorentz(x, y0, S, w, x0):
    return y0 + S * w / (4*(x-x0)**2 + w**2)

########################################  

#Kittel formula for out-of-plane FMR with uniaxial anisotropy along bias field
#w0 = gamma * (H0 + (k - 1)*Ms), where k = 2K1/mu0*Ms^2
#this function return f0
def kittel_oop_uniaxial(H0, K1_fit):
    return (gamma * H0 / (2*np.pi) + gamma * Ms * (2*K1_fit/(mu0 * Ms**2) - 1) / (2*np.pi))

#Kittel formula for in-plane FMR with uniaxial anisotropy along bias field
#w0 = gamma * sqrt((H0 + k*Ms) * (H0 + (1 + k)*Ms)), where k = 2K1/mu0*Ms^2
#this function return f0
def kittel_ip_uniaxial(H0, K1_fit):
    k = 2*K1_fit/(mu0 * Ms**2)
    return (gamma / (2*np.pi)) * np.sqrt((H0 + k*Ms)*(H0 + (1 + k)*Ms))

########################################

def setup_simulation_uniaxial(out_of_plane = True):
    
    """Setup FMR simulation for loaded thin-film material and uniaxial anisotropy along bias field"""
    
    meshdim = [80e-9, 80e-9, 2e-9]
    cellsize = [5e-9, 5e-9, 2e-9]
    
    output_file = 'fmr_peak.txt'
    
    #set thin-film Co/Pt with demag, exchange, Zeeman (enabled by default), and uniaxial anisotropy
    ns.setmaterial(material_name, meshdim)
    ns.cellsize(cellsize)
    #ns.addmodule(material_name, 'demag_N')
    ns.addmodule(material_name, 'aniuni')
    
    if out_of_plane:
        ns.setparam(material_name, 'ea1', [0,0,1])
    else:
        ns.setparam(material_name, 'ea1', [0,1,0])
    
    #thin film geometry : set pbc conditions
    ns.pbc(material_name, 'x', 10)
    ns.pbc(material_name, 'y', 10)
    
    #setup sinc pulse using a formula
    ns.setstage('Hequation')
    
    if out_of_plane:
        ns.editstagevalue(0, 'He * sinc(2*PI*fc*(t-t0)), 0, H0')
    else:
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
    ns.setdata('<M>', material_name)
    ns.savedatafile(output_file)
    
    ns.iterupdate(0)
    
    #this type of simulations needs to be run using a fixed time-step method so the magnetisation data is saved precisely at the sampling points
    #RK4 is good for this type of work, just set a low enough time-step
    ns.setode('LLG', 'RK4')
    #make sure the time-step divides the sampling time
    ns.setdt(500e-15)
    
    #ns.cuda(1)

########################################

#H0: bias field strength (A/m)
#fc: frequency cutoff (Hz)
def simulate_fmr_peak(H0, out_of_plane = True, show_result = False):
    
    """Simulate FMR peak for given bias field strength and return centre frequency and FWHM"""
    
    if out_of_plane:
        ns.setangle(0, 0)
    else:
        ns.setangle(90, 90)
    
    #Bias field (A/m)
    ns.equationconstants('H0', H0)
    
    ns.reset()
    ns.Run()
    
    #Analyse
    
    #get 2D list as position along horizontal, time along vertical
    Mav = ns.Get_Data_Columns(ns.savedatafile())
    Ms = ns.setparam(material_name, 'Ms')
    Mx = [row[0]/Ms for row in Mav]
    
    #2D fft
    fourier_data = np.fft.fftshift(np.abs(sp.fftpack.fft(Mx))**2)
    
    #get value ranges
    freq_len = len(fourier_data)
    freq = sp.fftpack.fftfreq(freq_len, time_step)
    
    #extract result from fourier data in a plottable form
    result = fourier_data[int(freq_len/2):freq_len]
    
    #use Lorentz fitting procedure built into Boris - curve_fit from scipy.optimize fails miserably even if you set bounds; routine in Boris simply works! (and doesn't need any bounds)
    xy = []
    xy.append(freq[0:len(result)])
    xy.append(result)
    ns.Save_Data_Columns('fmr_temp.txt', xy)
    ns.dp_load('fmr_temp.txt', [0, 1, 0, 1])
    #parameters returned as S, x0, w, y0, sigS, sigx0, sigw, sigy0
    lorentz_params = ns.dp_fitlorentz(0, 1)
    lorentz = [func_lorentz(f, lorentz_params[3], lorentz_params[0], lorentz_params[2], lorentz_params[1]) for f in freq[0:len(result)]]
    
    if show_result:
    
        print('f0 (GHz) = %.2f, FWHM (GHz) = %.2f' % (lorentz_params[1] / 1e9, lorentz_params[2] / 1e9))    
        gdamping = lorentz_params[2] / (2 * lorentz_params[1])
        print('damping = %.4f' % gdamping)
        plt.axes(xlabel = 'f (Hz)', ylabel = 'FMR Signal (a.u.)', title = 'Frequency-swept FMR')
        plt.plot(freq[0:len(result)], result, '.')
        plt.plot(freq[0:len(result)], lorentz)
        plt.show()
        
    return lorentz_params[1], lorentz_params[2]

########################################    

#Setup    

#set to False to simulate in-plane FMR
out_of_plane_fmr = True

material_name = 'Co/Pt'
#cutoff frequency (Hz)
fc = 200e9
#time step for saving magnetisation (s) : determined by Nyquist criterion from cutoff frequency
time_step = (0.5 / fc)
#total time to simulate; increasing this gives you more frequency points in the transform
total_time = 10000e-12
#Field and excitation field strength
He = 1000

setup_simulation_uniaxial(out_of_plane_fmr)

#useful constants
#saturation magentisation (A/m)
Ms = ns.setparam(material_name, 'Ms')
#anisotropy constant (for Co/Pt the easy axis is in the z direction)
K1 = ns.setparam(material_name, 'K1')
#gyromagnetic ratio modulus (times mu0)
gamma = 2.212761569e5 * ns.setparam(material_name, 'grel')
#Gilbert damping
damping = ns.setparam(material_name, 'damping')
mu0 = 4*np.pi*1e-7

########################################

#Simulate

H0range = np.arange(1e5, 1e6, 1e5)
f0range = []
dfrange = []

#simulate FMR peaks for given field range
for H0 in H0range:

    f0, df = simulate_fmr_peak(H0, out_of_plane_fmr, True)
    f0range.append(f0)
    dfrange.append(df)

########################################

#Analyse

if out_of_plane_fmr:
    popt, pcov = curve_fit(kittel_oop_uniaxial, H0range, f0range)
else:
    popt, pcov = curve_fit(kittel_ip_uniaxial, H0range, f0range)
    
K1_fit = popt[0]
perr = np.sqrt(np.diag(pcov))
std_K1_fit = perr[0]
print('K1 fitted (J/m^3) = %f +/- %0.4f' % (K1_fit, std_K1_fit))

plt.axes(xlabel = 'H (A/m)', ylabel = 'f (Hz)', title = 'FMR Kittel relation')
plt.plot(H0range, f0range, 'o')

if out_of_plane_fmr:
    plt.plot(H0range, [kittel_oop_uniaxial(H0, K1_fit) for H0 in H0range], '--')
else:
    plt.plot(H0range, [kittel_ip_uniaxial(H0, K1_fit) for H0 in H0range], '--')
    
plt.show()

#df = 2*damping*f0 + gamma * dH0, where dH0 is the inhomogeneous broadening
plt.axes(xlabel = 'f (Hz)', ylabel = 'delta f (Hz)', title = 'Frequency-swept FMR damping')
plt.plot(f0range, dfrange, 'o')

popt, pcov = curve_fit(lambda x, gradient, intercept: gradient*x + intercept, f0range, dfrange)
fitted_damping = popt[0] / 2
dH0 = popt[1] / gamma
perr = np.sqrt(np.diag(pcov))
std_fitted_damping = perr[0]

print('damping fitted = %.4f +/- %.4f' % (fitted_damping, std_fitted_damping))

plt.plot(f0range, [fitted_damping*2*f + dH0*gamma for f in f0range], '--')
plt.show()

    

    

    


