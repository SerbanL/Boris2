from NetSocks import NSClient, customize_plots
import matplotlib.pyplot as plt
import numpy as np

#setup communication with server
ns = NSClient(); ns.configure(True)
customize_plots()

########################################

#u = |J|*P*muB/(e*Ms*(1+beta^2))
P = ns.setparam('permalloy', 'P')
Ms = ns.setparam('permalloy', 'Ms')
beta = ns.setparam('permalloy', 'beta')
alpha = ns.setparam('permalloy', 'damping')
muB = 9.274009994e-24
e = 1.60217662e-19

data_noOe = ns.Get_Data_Columns('dwvelocity_vs_Jc_noOe.txt', [0, 1, 2])
u_data_noOe = [-J*P*muB/(e*Ms*(1 + beta**2)) for J in data_noOe[0]]

data_Oe = ns.Get_Data_Columns('dwvelocity_vs_Jc_withOe.txt', [0, 1, 2])
u_data_Oe = [-J*P*muB/(e*Ms*(1 + beta**2)) for J in data_Oe[0]]

fig, axes = plt.subplots(nrows = 2, ncols = 1, figsize = (6, 9))
fig.subplots_adjust(hspace = 0.25, wspace = 0)

axes[0].set(xlabel = 'u (m/s)', ylabel = 'DW Velocity (m/s)')
axes[0].plot(u_data_noOe, data_noOe[1], 'o', label = 'no Oe')
axes[0].plot(u_data_Oe, data_Oe[1], 's', label = 'with Oe')
axes[0].plot(u_data_noOe, [(beta/alpha) * u for u in u_data_noOe], '-', color = 'black', label = 'Analytical')
axes[0].annotate('(a)', xy=(0, 1), xytext=(10, -10), va='top', xycoords='axes fraction', textcoords='offset points', size = 16)
axes[0].legend(loc = 'lower right')

umag_noOe, umag_Oe = [], []
vdiff_noOe, vdiff_Oe = [], []
err_noOe, err_Oe = [], []

for i in range(int(len(u_data_noOe)/2)):
    
    vdiff_noOe.append(np.abs(data_noOe[1][i+int(len(u_data_noOe)/2)]) - np.abs(data_noOe[1][i]))
    umag_noOe.append(np.abs(u_data_noOe[i]))
    err_noOe.append(np.sqrt(data_noOe[2][i]**2 + data_noOe[2][int(len(u_data_noOe)/2)]**2))
    
    vdiff_Oe.append(np.abs(data_Oe[1][i+int(len(u_data_Oe)/2)]) - np.abs(data_Oe[1][i]))
    umag_Oe.append(np.abs(u_data_Oe[i]))
    err_Oe.append(np.sqrt(data_Oe[2][i]**2 + data_Oe[2][int(len(u_data_Oe)/2)]**2))

axes[1].set(xlabel = 'u (m/s)', ylabel = 'DW Velocity asymmetry (m/s)')
axes[1].errorbar(umag_noOe, vdiff_noOe, yerr = err_noOe, fmt = 'o', label = 'no Oe')
axes[1].errorbar(umag_Oe, vdiff_Oe, yerr = err_Oe, fmt = 's', label = 'with Oe')
axes[1].annotate('(b)', xy=(0, 1), xytext=(10, -10), va='top', xycoords='axes fraction', textcoords='offset points', size = 16)
axes[1].legend(loc = 'lower right')

plt.savefig('Figure11.2.png')
plt.show()

