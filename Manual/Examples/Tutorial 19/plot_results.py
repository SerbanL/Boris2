from NetSocks import NSClient, customize_plots
import matplotlib.pyplot as plt
import numpy as np

#setup communication with server
ns = NSClient(); ns.configure(True)
customize_plots()

########################################

data = ns.Get_Data_Columns('ishefmr_data.txt', [0, 1, 4, 7, 8])
time_ns = np.array(data[0]) / 1e-9
Hrf = np.array(data[1]) / 1e3
Mx = np.array(data[2]) / 1e3
Vm = np.array(data[3]) / 1e-9
Vp = np.array(data[4]) / 1e-9

fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize = (12, 4))
fig.subplots_adjust(hspace = 0, wspace = 0.25)

axes[0].set(xlabel = 'Time (ns)', ylabel = 'H, M (kA/m)')

axes[0].plot(time_ns, Hrf, label = 'H$_{rf}$')
axes[0].plot(time_ns, Mx, label = 'M$_{x}$')
axes[0].legend()
axes[0].annotate('(a)', xy=(0, 1), xytext=(10, -10), va='top', xycoords='axes fraction', textcoords='offset points', size = 16)

axes[1].set(xlabel = 'Time (ns)', ylabel = 'Voltage (nV)')
axes[1].plot(time_ns, Vm, label = 'V (-y)')
axes[1].plot(time_ns, Vp, label = 'V (+y)')
axes[1].legend()
axes[1].annotate('(b)', xy=(0, 1), xytext=(10, -10), va='top', xycoords='axes fraction', textcoords='offset points', size = 16)

plt.savefig('Figure19.2.png')
plt.show()
