from NetSocks import NSClient, customize_plots
import numpy as np
import matplotlib.pyplot as plt

customize_plots()
ns = NSClient(); ns.configure(True)

##############################################
#
# Exercise 2.1

ns.loadsim('hysteresis')
ns.Run()

data = ns.Get_Data_Columns('hysteresis.txt', [4, 7])

plt.axes(xlabel = 'H (kA/m)', ylabel = 'M (kA/m)')
plt.plot(np.array(data[0]) / 1e3, np.array(data[1]) / 1e3, 'o-', color = 'black')
plt.savefig('Figure2.2.png')
plt.show()

##############################################
#
# Exercise 2.2a

ns.loadsim('muMAG4 field1')
ns.Run()

data = ns.Get_Data_Columns('mumag4_field1.txt', [3, 7, 8, 9])

plt.axes(xlabel = 'Time (ns)', ylabel = 'M (kA/m)')
plt.plot(np.array(data[0]) / 1e-9, np.array(data[1]) / 1e3, '-', color = 'red', label = 'M$_x$')
plt.plot(np.array(data[0]) / 1e-9, np.array(data[2]) / 1e3, '-', color = 'blue', label = 'M$_y$')
plt.plot(np.array(data[0]) / 1e-9, np.array(data[3]) / 1e3, '-', color = 'green', label = 'M$_z$')
plt.legend()
plt.savefig('Figure2.4a.png')
plt.show()

##############################################
#
# Exercise 2.2b

ns.loadsim('muMAG4 field2')
ns.Run()

data = ns.Get_Data_Columns('mumag4_field2.txt', [3, 7, 8, 9])

plt.axes(xlabel = 'Time (ns)', ylabel = 'M (kA/m)')
plt.plot(np.array(data[0]) / 1e-9, np.array(data[1]) / 1e3, '-', color = 'red', label = 'M$_x$')
plt.plot(np.array(data[0]) / 1e-9, np.array(data[2]) / 1e3, '-', color = 'blue', label = 'M$_y$')
plt.plot(np.array(data[0]) / 1e-9, np.array(data[3]) / 1e3, '-', color = 'green', label = 'M$_z$')
plt.legend()
plt.savefig('Figure2.5a.png')
plt.show()