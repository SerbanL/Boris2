from NetSocks import NSClient, customize_plots
import numpy as np
import matplotlib.pyplot as plt

customize_plots()
ns = NSClient(); ns.configure(True)

##############################################
#
# Exercise 8.3a

ns.loadsim('amr_longitudinal')
ns.Run()

##############################################
#
# Exercise 8.3b

ns.loadsim('amr_transverse')
ns.Run()

##### Plotting

data = ns.Get_Data_Columns('longitudinal_mr.txt', [0, 1, 2, 3])
Hl = [Hx * np.cos(np.radians(5)) + Hy * np.sin(np.radians(5)) for (Hx, Hy) in zip(data[0], data[1])]
Rl = data[3]

data = ns.Get_Data_Columns('transverse_mr.txt', [0, 1, 2, 3])
Ht = [Hx * np.cos(np.radians(85)) + Hy * np.sin(np.radians(85)) for (Hx, Hy) in zip(data[0], data[1])]
Rt = data[3]

plt.axes(xlabel = 'H (kA/m)', ylabel = 'R (Ohms)')
plt.plot(np.array(Hl) / 1e3, Rl, '-', color = 'blue', label = 'longitudinal')
plt.plot(np.array(Ht) / 1e3, Rt, '-', color = 'red', label = 'transverse')
plt.legend()
plt.savefig('Figure8.2.png')
plt.show()