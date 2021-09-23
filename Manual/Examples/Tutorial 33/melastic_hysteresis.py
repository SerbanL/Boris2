"""
This script is part of Boris Computational Spintronics v3.0

@author: Serban Lepadatu, 2020
"""

from NetSocks import NSClient, customize_plots
import matplotlib.pyplot as plt

#setup communication with server
ns = NSClient(); ns.configure(True); customize_plots()

########################################

#stress magnitude in Pa
stress_mag = 4e8

########################################

#1. No stress

ns.loadsim('melastic')

stress = 0.0
ns.setstress(stress, 90, 0)
output_file_nostress = 'melastic_stress_%.f.txt' % stress
ns.savedatafile(output_file_nostress)
ns.Run()

plt.axes(xlabel = 'H (A/m)', ylabel = 'Mx  (A/m)')
hyster_nostress = ns.Get_Data_Columns(output_file_nostress, [0, 3])
plt.plot([H/1e3 for H in hyster_nostress[0]], [M/1e3 for M in hyster_nostress[1]], label = 'no stress')
plt.show()

########################################

#2. Compressive stress (hence reinforces easy axis along field)

ns.loadsim('melastic')

stress = -stress_mag
ns.setstress(stress, 90, 0)
output_file_compressive = 'melastic_stress_%.f.txt' % stress
ns.savedatafile(output_file_compressive)
ns.Run()

plt.axes(xlabel = 'H (kA/m)', ylabel = 'Mx  (kA/m)')
hyster_compressive = ns.Get_Data_Columns(output_file_compressive, [0, 3])
plt.plot([H/1e3 for H in hyster_compressive[0]], [M/1e3 for M in hyster_compressive[1]], label = 'compressive stress: %.1f (MPa)' % (stress_mag/1e6))
plt.show()

########################################

#3. Extensive stress (hence weakens easy axis along field)

ns.loadsim('melastic')

stress = stress_mag
ns.setstress(stress, 90, 0)
output_file_extensive = 'melastic_stress_%.f.txt' % stress
ns.savedatafile(output_file_extensive)
ns.Run()

plt.axes(xlabel = 'H (A/m)', ylabel = 'Mx  (A/m)')
hyster_extensive = ns.Get_Data_Columns(output_file_extensive, [0, 3])
plt.plot([H/1e3 for H in hyster_extensive[0]], [M/1e3 for M in hyster_extensive[1]], label = 'extensive stress: %.1f (MPa)' % (stress_mag/1e6))
plt.show()

########################################

plt.axes(xlabel = 'H (A/m)', ylabel = 'Mx  (A/m)')

hyster_nostress = ns.Get_Data_Columns(output_file_nostress, [0, 3])
plt.plot([H/1e3 for H in hyster_nostress[0]], [M/1e3 for M in hyster_nostress[1]], label = 'no stress')

hyster_compressive = ns.Get_Data_Columns(output_file_compressive, [0, 3])
plt.plot([H/1e3 for H in hyster_compressive[0]], [M/1e3 for M in hyster_compressive[1]], label = 'compressive stress: %.1f (MPa)' % (stress_mag/1e6))

hyster_extensive = ns.Get_Data_Columns(output_file_extensive, [0, 3])
plt.plot([H/1e3 for H in hyster_extensive[0]], [M/1e3 for M in hyster_extensive[1]], label = 'extensive stress: %.1f (MPa)' % (stress_mag/1e6))

plt.legend(loc = 'lower left')
plt.savefig('melastic.png', dpi = 600)
plt.show()
