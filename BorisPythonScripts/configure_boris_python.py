import os, sys
#get location of site-packages directory
site_packages = next(p for p in sys.path if 'site-packages' in p)

#make a path file which will contain paths to all "standard" Boris python modules
#these can then be imported using a single import command as needed
path_file = open(os.path.join(site_packages, 'boris_python_scripts.pth'), "w")

BorisPythonScripts_dir = 'BorisPythonScripts/'

dirlist = [str(s[0])+'/' for p in ['~/Documents/Boris Data/' + BorisPythonScripts_dir, '~/Documents/Boris_Data/' + BorisPythonScripts_dir] for s in os.walk(os.path.expanduser(p))]

for directory in dirlist:
    path_file.write(os.path.expanduser(directory + '\n'))
    
path_file.write(os.path.expanduser('~/Documents/Boris Data/') + '\n')
path_file.write(os.path.expanduser('~/Documents/Boris_Data/') + '\n')
    
path_file.close()



