from NetSocks import NSClient, customize_plots
import matplotlib.pyplot as plt
import numpy as np

#setup communication with server
ns = NSClient(); ns.configure(True)
customize_plots()

########################################

loop_saf = ns.Get_Data_Columns('saf_hysteresis_loop.txt', [0, 1])
loop_syf = ns.Get_Data_Columns('syf_hysteresis_loop.txt', [0, 1])

plt.axes(xlabel = 'H (kA/m)', ylabel = 'M (kA/m)')

plt.plot(np.array(loop_saf[0]) / 1e3, np.array(loop_saf[1]) / 1e3, label = 'SAF')
plt.plot(np.array(loop_syf[0]) / 1e3, np.array(loop_syf[1]) / 1e3, label = 'SyF')
plt.legend()
plt.savefig('Figure12.2.png')
plt.show()
