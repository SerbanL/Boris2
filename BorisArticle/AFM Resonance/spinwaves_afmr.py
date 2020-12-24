"""
This script is part of Boris Computational Spintronics Article

Generates data in Figure 10(b).

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient
import numpy as np
import scipy as sp
import scipy.fftpack
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.signal

ns = NSClient()
ns.configure(True)

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
mpl.rcParams['legend.fontsize'] = 10.5
mpl.rcParams['legend.edgecolor'] = 'black'
mpl.rcParams['legend.title_fontsize'] = 15
mpl.rcParams['legend.labelspacing'] = 0.1
mpl.rcParams['legend.borderpad'] = 0.3

#math type
mpl.rcParams['mathtext.default'] = 'regular'

mpl.rcParams['savefig.dpi'] = 600
mpl.rcParams['savefig.bbox'] = 'tight'

########################################

#Simulation script here

#Setup

#fmr cutoff frequency (Hz)
fc = 1000e9
#time step for saving magnetisation (s) : determined by Nyquist criterion from cutoff frequency
time_step = (0.5 / fc)

#k vector cutoff (rad / m)
kc = 2*np.pi*0.1255e9

#the mesh step used to obtain M values along the k vector direction; again this must be set from Nyquist criterion
#thus mesh_step = 0.5 / (kc / 2PI), which gives 4 nm
mesh_step = 4e-9

meshdim = [1000e-9, 50e-9, 1e-9]
cellsize = [1e-9, 2e-9, 1e-9]

#sinc pulse time centre (s) and total time to simulate (s)
#increase t0 if you want better resolution
t0 = 200e-12
total_time = 2 * t0

#bias field and excitation strength (A/m)
H0 = 0
He = 40e3

output_file = 'spinwaves_afmr.txt'

#define useful constants
A = 0.5e-11
Ah = -20e6
Anh = -10e-12
Ms = 800e3
K1 = 50e3
gamma = 2.212761569e5
damping = 0.002
mu0 = 4*np.pi*1e-7

ns.setafmesh('AF', meshdim)

#setup mesh dimensions, pbc (10 images along x only), set magnetisation along +x direction
#ns.meshrect(meshdim)
ns.cellsize(cellsize)
ns.pbc('AF', 'x', 10)
ns.setangle(90, 0)

#setup sinc pulse using a formula
ns.setstage('Hequation')
ns.editstagevalue(0, 'H0, He * sinc(kc*(x-Lx/2))*sinc(kc*(y-Ly/2))*sinc(2*PI*fc*(t-t0)), 0')

#define the equation constants
#Bias field (A/m)
ns.equationconstants('H0', H0)
#Excitation field (A/m)
ns.equationconstants('He', He)
#k cutoff (rad/m)
ns.equationconstants('kc', kc)
#f cutoff (Hz) -> hence time step of 1ps due to Nyquist criterion
ns.equationconstants('fc', fc)
#time center (s)
ns.equationconstants('t0', t0)

#set all other parameters given above
ns.setparam('AF', 'damping_AFM', [damping, damping])
ns.setparam('AF', 'Ms_AFM', [Ms, Ms])
ns.setparam('AF', 'K1_AFM', [K1, K1])
ns.setparam('AF', 'Ah', [Ah, Ah])
ns.setparam('AF', 'Anh', [Anh, Anh])
ns.setparam('AF', 'A_AFM', [A, A])
ns.setparam('AF', 'grel_AFM', [gamma / 2.212761569e5, gamma / 2.212761569e5])

#make sure output file is wiped clean
ns.dp_newfile(output_file)

#set fixed time-step RK4 method with a small enough time step (5 fs is sufficient - small value due to 1 nm cellsize)
ns.setode('LLG', 'RK4')
ns.setdt(50e-15)

########################################

#Run

time = 0.0
ns.cuda(1)

#simulate, saving data every time_step
while time < total_time:
    
    ns.editstagestop(0, 'time', time + time_step)
    ns.Run()
    
    #get magnetisation profile along length through center
    ns.dp_getexactprofile([cellsize[0]/2, meshdim[1]/2 + cellsize[1]/2, 0], [meshdim[0] - cellsize[0]/2, meshdim[1]/2 + cellsize[1]/2, 0], mesh_step, 0)
    #save only the y component of magnetisation at time_step intervals
    ns.dp_div(2, Ms)
    ns.dp_saveappendasrow(output_file, 2)
    
    time += time_step

########################################

#Analyse

#get 2D list as position along horizontal, time along vertical
pos_time = ns.Get_Data_Columns(output_file)

#2D fft (use log function to improve image contrast)
fourier_data = np.fft.fftshift(np.log(np.abs(sp.fftpack.fft2(pos_time))))

#get value ranges
freq_len = len(fourier_data)
k_len = len(fourier_data[0])
freq = sp.fftpack.fftfreq(freq_len, time_step)
kvector = sp.fftpack.fftfreq(k_len, mesh_step)

#maximum k and f values
k_max = 2*np.pi*kvector[int(0.5 * len(kvector))] / 1e9
f_min_GHz = np.abs(freq[0]) / 1e9
f_max_GHz = np.abs(freq[int(0.5 * len(freq))]) / 1e9
f_points = int(0.5 * freq_len)

#extract result from fourier data in a plottable form
result = [fourier_data[i]/1e9 for i in range(int(0.5 * freq_len),freq_len)]

#plot spin wave dispersion
plt.axes(xlabel = 'k (rad/nm)', ylabel = 'f (GHz)')
plt.imshow(result, origin='lower', interpolation='bilinear', 
           extent = [-k_max, k_max, f_min_GHz, f_max_GHz], 
           aspect ="auto")

#plot analytical lines
krange = np.arange(-k_max, k_max, 0.01)

f = []

for k in krange:
    
    He = 4 * np.abs(Ah) / (mu0 * Ms)
    Ha = 2 * K1 / (mu0 * Ms)
    cA = 2*A/(mu0 * Ms)
    cAnh = np.abs(Anh)/(mu0 * Ms)
    fval = gamma * np.sqrt((2*He + Ha + (cA + cAnh)*(k*1e9)**2)*(Ha + (cA + cAnh)*(k*1e9)**2)) / (1e9 * 2*np.pi)
    f.append(fval)

plt.plot(krange, f, 'r:')

#restrict plotting limits
plt.ylim(f_min_GHz, f_max_GHz)
plt.xlim(-k_max, k_max)

plt.savefig('spinwaves_afmr.png', dpi = 600)
plt.show()
    

    

    


