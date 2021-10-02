#NetSocks Module Updated on : 29/09/2021
#Boris version : 3.4
import sys
import os
import subprocess
import platform
import socket
import matplotlib.pyplot as plt
import numpy as np
from itertools import product
import struct
import time
from copy import deepcopy

import matplotlib as mpl

########################################

#plots customizations; see https://matplotlib.org/3.2.1/tutorials/introductory/customizing.html

def customize_plots(labelsize = 20, ticklabelsize = 15, legendfontsize = 13, savefigdpi = 600):

    #axes
    mpl.rcParams['axes.linewidth'] = 2.0
    mpl.rcParams['axes.labelsize'] = labelsize
    mpl.rcParams['font.family'] = 'Arial'
    
    #ticks
    mpl.rcParams['xtick.labelsize'] = ticklabelsize
    mpl.rcParams['ytick.labelsize'] = ticklabelsize
    
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
    mpl.rcParams['legend.fontsize'] = legendfontsize
    mpl.rcParams['legend.edgecolor'] = 'black'
    mpl.rcParams['legend.title_fontsize'] = 15
    mpl.rcParams['legend.labelspacing'] = 0.1
    mpl.rcParams['legend.borderpad'] = 0.3
    
    #math type
    mpl.rcParams['mathtext.default'] = 'regular'
    
    mpl.rcParams['savefig.dpi'] = savefigdpi
    mpl.rcParams['savefig.bbox'] = 'tight'
    
    #font
    mpl.rcParams['font.size'] = 12

###############################################################################################################################

class ElementaryShape:
    
    #################### DATA
    
    #name of elementary shape: disk, rect, triangle, ellipsoid, pyramid, tetrahedron, cone, torus
    name = ""
    
    #dimensions in metres
    dimensions = np.array([0.0, 0.0, 0.0])
    
    #shape centre coordinates position, relative to mesh
    position = np.array([0.0, 0.0, 0.0])
    
    #rotation in degrees as psi (around y), theta (around x), phi (around z)
    rotation = np.array([0.0, 0.0, 0.0])
    
    #number of x, y, z repetitions for generating arrays
    repetitions = np.array([1, 1, 1])
    
    #x, y, z displacement in metres used when generating arrays
    displacement = np.array([0.0, 0.0, 0.0])
    
    #method used to draw shape: add, sub, xor, and
    method = "add"
    
    #################### CTOR
    
    def __init__(
            self,
            name = "",
            dimensions = np.array([0.0, 0.0, 0.0]), 
            position = np.array([0.0, 0.0, 0.0]),
            rotation = np.array([0.0, 0.0, 0.0]),
            repetitions = np.array([1, 1, 1]),
            displacement = np.array([0.0, 0.0, 0.0]),
            method = "add"): 
        
        self.name = deepcopy(name)
        self.dimensions = np.array(deepcopy(dimensions))
        self.position = np.array(deepcopy(position))
        self.rotation = np.array(deepcopy(rotation))
        self.repetitions = np.array(deepcopy(repetitions))
        self.displacement = np.array(deepcopy(displacement))
        self.method = np.array(deepcopy(method))
        
    #################### DIMENSIONS
        
    def setdimensions(self, dimensions):
        """set dimensions of elementary shape"""
        self.dimensions = np.array(deepcopy(dimensions))
        
    def scale(self, scalefactors):
        """scale dimensions of elementary shape"""
        self.dimensions *= np.array(scalefactors)
            
    #################### POSITION
            
    def setposition(self, position):
        """set position of elementary shape"""
        self.position = np.array(deepcopy(position))
        
    def move(self, positionshift):
        """translate position of elementary shape"""
        self.position += np.array(positionshift)
        
    #################### ROTATION
            
    def setrotation(self, rotation):
        """set rotation of elementary shape"""
        self.rotation = np.array(deepcopy(rotation))
        
    def rotate(self, rotation):
        """rotate elementary shape"""
        self.rotation += np.array(rotation)
        
    #################### REPETITIONS
            
    def setrepetitions(self, repetitions, displacement):
        """set number of repetitions and displacement of elementary shape for generating array"""
        self.repetitions = np.array(deepcopy(repetitions))
        self.displacement = np.array(deepcopy(displacement))
        
    #################### METHOD
    
    def setaddshape(self):
        self.method = "add"
        
    def setsubshape(self):
        self.method = "sub"
        
    def is_additive(self):
        return self.method == "add"
    
    def is_subtractive(self):
        return self.method == "sub"
    
    #################### CONVERSION
    
    def tostring(self):
        lst = [self.dimensions, self.position, self.rotation, self.repetitions, self.displacement]
        text_lst = [self.name] + [" ".join(map(str, elem)) for elem in lst] + [self.method]
        return " ".join(map(str, text_lst))

###############################################################################################################################

class Shape:
    
    #################### DATA
    
    #List of elementary shapes
    shapes = []
    
    #################### CTOR
    
    def __init__(self, shape = ElementaryShape()):
        self.shapes = [deepcopy(shape)]
        
    #################### OPERATORS
        
    #Add two shapes, returning new copy
    def __add__(self, shape_right):
        
        newshape = Shape()
        newshape.shapes = deepcopy(self.shapes) + deepcopy(shape_right.shapes)
        return newshape
    
    #Subtract two shapes, returning new copy
    def __sub__(self, shape_right):
        
        newshape_left = deepcopy(self)
        newshape_right = deepcopy(shape_right)
        for shape in newshape_right.shapes:
            if shape.is_additive(): shape.setsubshape()
            elif shape.is_subtractive(): shape.setaddshape()
        return newshape_left + newshape_right
      
    #################### AUXILIARY
    
    def rotate_object_yxz(self, r, psi_theta_phi_deg):
        
        psi, theta, phi = psi_theta_phi_deg[0] * np.pi / 180, psi_theta_phi_deg[1] * np.pi / 180, psi_theta_phi_deg[2] * np.pi / 180
        rr = np.array(r)
        
        rr[0] = (np.cos(psi) * np.cos(phi) + np.sin(psi) * np.sin(theta) * np.sin(phi)) * r[0] + (np.cos(phi) * np.sin(psi) * np.sin(theta) - np.cos(psi) * np.sin(phi)) * r[1] + (np.cos(theta) * np.sin(psi)) * r[2]
        rr[1] = (np.cos(theta) * np.sin(phi)) * r[0] + (np.cos(theta) * np.cos(phi)) * r[1] - np.sin(theta) * r[2]
        rr[2] = (np.cos(psi) * np.sin(theta) * np.sin(phi) - np.cos(phi) * np.sin(psi)) * r[0] + (np.cos(psi) * np.cos(phi) * np.sin(theta) + np.sin(psi) * np.sin(phi)) * r[1] + (np.cos(psi) * np.cos(theta)) * r[2]
        return rr
    
    #################### DIMENSIONS
    
    def setdimensions(self, dimensions):
        """set dimensions of first elementary shape contained, and scale everything else in proportion"""
        scalefactors = np.array(dimensions) / self.shapes[0].dimensions
        self.scale(scalefactors)
        return self
        
    def scale(self, scalefactors):
        """scale dimensions of shape"""
        current_position = self.shapes[0].position
        for shape in self.shapes: 
            shape.scale(scalefactors)
            shape.setposition(current_position + (shape.position - current_position) * scalefactors)
        return self
    
    #################### POSITION    
    
    def setposition(self, position):
        """set position of shape, defined by the position of the first elementary shape contained"""
        current_position = self.shapes[0].position
        self.shapes[0].setposition(position)
        for shape in self.shapes[1:]: shape.move(position - current_position)
        return self
        
    def move(self, positionshift):
        """translate position of shape"""
        for shape in self.shapes: shape.move(positionshift)
        return self
        
    #################### ROTATION
    
    def setrotation(self, rotation):
        """set rotation of shape, around shape position as defined by the first elementary shape contained"""
        current_position = self.shapes[0].position
        for shape in self.shapes:
            shape.setposition(current_position + self.rotate_object_yxz(shape.position - current_position, rotation))
            shape.setrotation(rotation)
        return self
        
    def rotate(self, rotation):
        """rotate shape around shape position as defined by the first elementary shape contained"""
        current_position = self.shapes[0].position
        for shape in self.shapes:
            shape.setposition(current_position + self.rotate_object_yxz(shape.position - current_position, rotation))
            shape.rotate(rotation)
        return self
    
    #################### REPETITIONS
    
    def setrepetitions(self, repetitions, displacement):
        """set number of repetitions and displacement of shape for generating array"""
        for shape in self.shapes: shape.setrepetitions(repetitions, displacement)
        return self
        
    #################### CONVERSION
    
    def tostring(self):
        shapes_text = [shape.tostring() for shape in self.shapes]
        return " ".join(shapes_text)
        
    #################### ELEMENTARY SHAPES GENERATORS
        
    #define an elementary disk shape
    def disk(dimensions = np.array([0.0, 0.0, 0.0]), position = np.array([0.0, 0.0, 0.0])):
        return Shape(ElementaryShape("disk", dimensions, position))
    
    #define an elementary rectangle shape
    def rect(dimensions = np.array([0.0, 0.0, 0.0]), position = np.array([0.0, 0.0, 0.0])):
        return Shape(ElementaryShape("rect", dimensions, position))
    
    #define an elementary triangle shape
    def triangle(dimensions = np.array([0.0, 0.0, 0.0]), position = np.array([0.0, 0.0, 0.0])):
        return Shape(ElementaryShape("triangle", dimensions, position))
    
    #define an elementary ellipsoid shape
    def ellipsoid(dimensions = np.array([0.0, 0.0, 0.0]), position = np.array([0.0, 0.0, 0.0])):
        return Shape(ElementaryShape("ellipsoid", dimensions, position))
    
    #define an elementary pyramid shape
    def pyramid(dimensions = np.array([0.0, 0.0, 0.0]), position = np.array([0.0, 0.0, 0.0])):
        return Shape(ElementaryShape("pyramid", dimensions, position))
    
    #define an elementary tetrahedron shape
    def tetrahedron(dimensions = np.array([0.0, 0.0, 0.0]), position = np.array([0.0, 0.0, 0.0])):
        return Shape(ElementaryShape("tetrahedron", dimensions, position))
    
    #define an elementary cone shape
    def cone(dimensions = np.array([0.0, 0.0, 0.0]), position = np.array([0.0, 0.0, 0.0])):
        return Shape(ElementaryShape("cone", dimensions, position))
    
    #define an elementary torus shape
    def torus(dimensions = np.array([0.0, 0.0, 0.0]), position = np.array([0.0, 0.0, 0.0])):
        return Shape(ElementaryShape("torus", dimensions, position))

        
###############################################################################################################################

class NSClientConfig:
    
    scriptserverip = 'localhost'
    scriptserverport = 1542
    scriptserverpwd = ''
    
    cudaDevice = 0
    
    #verbosity of Boris console
    verbose = False
    
    boris_path = ''
    
    boris_exe = ''
    
    window = ''
    
    verbose = False
    
    def __init__(self, 
                 scriptserverip = 'localhost', scriptserverport = 1542, scriptserverpwd = '', 
                 cudaDevice = -1, 
                 boris_path = '', boris_exe = '',
                 window = 'back',
                 verbose = False):
        
        self.scriptserverip = scriptserverip
        self.scriptserverport = scriptserverport
        self.scriptserverpwd = scriptserverpwd
        self.cudaDevice = cudaDevice
        self.boris_path = boris_path
        self.boris_exe = boris_exe
        self.window = window
        self.verbose = verbose

###############################################################################################################################

class NSClient:

    #################### DATA
    
    maxLenMessage = 4096
    timeout_ms = 30000000
    
    scriptserverip = 'localhost'
    scriptserverport = 1542
    scriptserverpwd = ''
    
    cudaDevice = 0
    
    #verbosity of Boris console
    verbose = False
    
    #verbosity of Python script (adjusted by configure)
    script_verbose = True

    #on Windows assume this is where Boris.exe is (should be if installed with installer)
    #on Linux don't attempt to define a default : user will have to provide path if they want automatic startup
    win_default_boris_path = 'C:/Program Files (x86)/Boris'

    #################### CTOR / DTOR

    def __initialize_nsclient(self, 
                              scriptserverip = 'localhost', scriptserverport = 1542, scriptserverpwd = '', 
                              cudaDevice = -1, 
                              boris_path = '', boris_exe = '',
                              window = 'back',
                              verbose = False):
    
        self.scriptserverip = scriptserverip
        if scriptserverport >= 0: self.scriptserverport = scriptserverport
        self.verbose = verbose
        self.scriptserverpwd = scriptserverpwd
        
        #if boris path not specified use default if using Windows
        #Linux users will need to provide path
        if boris_path == '':
            
            if platform.system() == 'Windows':
                boris_path = self.win_default_boris_path
                
        if boris_exe == '':
            
            if platform.system() == 'Windows':
                boris_exe = 'Boris.exe'
            
            else:
                boris_exe = './BorisLin'
            
        #start a new Boris instance if none exists listening on serverport, if we have a path to Boris.exe            
        if len(boris_path):
            
            if scriptserverport >= 0:
            
                with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            
                    try:
                        sock.connect((self.scriptserverip, self.scriptserverport))
                    except:
                        if scriptserverip == 'localhost':
                            print("No server found on port %d. Starting new instance." % self.scriptserverport)
                            os.chdir(boris_path)
                            subprocess.Popen([boris_exe, str(self.scriptserverport), str(cudaDevice), window, scriptserverpwd])
                            os.chdir(os.path.dirname(sys.argv[0]) + "/")
                            
                            #now make sure server is running and ready to accept input
                            for tryidx in range(10):
                            
                                with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock2:
                                    
                                    try:
                                        sock2.connect((self.scriptserverip, self.scriptserverport))
                                        print("New instance started.")
                                        break
                                    except:
                                        time.sleep(0.1)
                            else:
                                print("Couldn't start new instance with required server port - start it manually.");
                                    
                                    
                        else:
                            print("No server found on port %d. Make sure remote host has a Boris instance running for given port and is accessible." % self.scriptserverport)
                            
            else:
                
                with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            
                    self.scriptserverport = 1542
                    
                    while True:
                    
                        try:
                            sock.connect((self.scriptserverip, self.scriptserverport))
                            self.scriptserverport += 1
                        except:
                            if scriptserverip == 'localhost':
                                print("No server found on port %d. Starting new instance." % self.scriptserverport)
                                os.chdir(boris_path)
                                subprocess.Popen([boris_exe, str(self.scriptserverport), str(cudaDevice), window, scriptserverpwd])
                                os.chdir(os.path.dirname(sys.argv[0]) + "/")
                                
                                #now make sure server is running and ready to accept input
                                for tryidx in range(10):
                                
                                    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock2:
                                        
                                        try:
                                            sock2.connect((self.scriptserverip, self.scriptserverport))
                                            print("New instance started.")
                                            break
                                        except:
                                            time.sleep(0.1)
                                else:
                                    print("Couldn't start new instance with required server port - start it manually.");
                                        
                                        
                            else:
                                print("No server found on port %d. Make sure remote host has a Boris instance running for given port and is accessible." % self.scriptserverport)
                                
                            break
                        
        if cudaDevice > 0: 
            self.cudaDevice = cudaDevice
            self.selectcudadevice(cudaDevice)
    

    def __init__(self, 
                 scriptserverip = 'localhost', scriptserverport = 1542, scriptserverpwd = '', 
                 cudaDevice = -1, 
                 boris_path = '', boris_exe = '',
                 window = 'back',
                 verbose = False):
        
        self.__initialize_nsclient(scriptserverip , scriptserverport, scriptserverpwd, 
                 cudaDevice, 
                 boris_path, boris_exe,
                 window,
                 verbose)

    def __init__(self, nscfg = NSClientConfig()):
        
        self.__initialize_nsclient(
                nscfg.scriptserverip, nscfg.scriptserverport, nscfg.scriptserverpwd, 
                nscfg.cudaDevice, 
                nscfg.boris_path, nscfg.boris_exe,
                nscfg.window,
                nscfg.verbose)

    #################### AUXILIARY

    #use this on return parameters : returned parameters can either be numbers or a word (text without spaces)
    #the words by themselves cannot be converted to numbers
    #thus try to convert to a number first, if not must be a word
    def convert_returned_parameter(self, string):
        try:    return float(string)
        except: return string
    
    #SendCommand also serves as auxiliary method
    
    #check if we can convert the row text to list of numbers; expecting tab-spaced data
    def can_convert(self, row_text, separator = '\t'):
        
        for entry in row_text.rstrip().split(separator):
            try:    float(entry)
            except: return False

        return True
    
    #Set directory same as script directory if running on localhost only. Reset to default state unless specified otherwise
    def configure(self, reset_to_default = True, script_verbose = True):
        """Set directory same as script directory, and reset to default state if called with True; also set script verbosity"""
        
        self.script_verbose = script_verbose
        
        if self.scriptserverip == 'localhost':
            directory = os.path.dirname(sys.argv[0]) + "/"
            print("Working directory is: ", directory)
            if reset_to_default: self.default()
            self.chdir(directory)
            
        elif reset_to_default: self.default()
            
    
    #################### PLOTTING HELPERS
    
    #load columns from tab-spaced data file, e.g. as outputted by a Boris simulation
    def Get_Data_Columns(self, fileName, column_indexes = '', separator = '\t'):
        """Get indexed columns from tab-spaced data file as a list"""
        #Get data locally
        f = open(fileName, 'r')
        rows = [[float(number) for number in row.rstrip().split(separator)] for row in f.readlines() if self.can_convert(row, separator)]
        f.close()
        
        if isinstance(column_indexes, list):
            output_columns = []
            
            for index in column_indexes:
                output_columns.append([row[index] for row in rows])
            
            return output_columns
        
        else: 
            if not isinstance(column_indexes, str):
                return [row[column_indexes] for row in rows]
            else:
                return rows
        
    #Save data columns in file name, tab separated
    def Save_Data_Columns(self, fileName, data_columns):
        """Save data columns to tab-spaced data file"""
        f = open(fileName, 'w')
        
        max_len = 0
        
        for column in data_columns:
            max_len = max_len if len(column) < max_len else len(column)
            
        for row_idx in range(max_len):
            
            line = ''
            
            for col_idx in range(len(data_columns)):
                
                if row_idx < len(data_columns[col_idx]): 
                    line += str(data_columns[col_idx][row_idx])
                else: 
                    line += ' '
                
                if col_idx < len(data_columns) - 1:
                    line += '\t'
            
            f.write(line + '\n')
        
        f.close()
        
        
    #Simple plot of y vs x : helps if you just want to see a simple simulation output plot
    #you can get data from simulation output file and plot it in just 2 lines of code in your Python script (can be done in 1 line with 2 calls to Get_Data_Columns instead)
    def Plot_Data(self, x, y, xlabel = '', ylabel = '', title = '', label_ = '', imageFile = ''):
        """Simple plotting helper"""
        plt.axes(xlabel = xlabel, ylabel = ylabel, title = title)
        plt.grid()
        plt.plot(x, y, label = label_)
        if len(label_): plt.legend()
        
        if len(imageFile): plt.savefig(imageFile + '.png')
        
        plt.show()
        
    def PlotPolar_Data(self, r, theta_deg, xlabel = '', ylabel = '', title = '', label_ = '', imageFile = ''):
        """Simple plotting helper: polar"""
        plt.axes(xlabel = xlabel, ylabel = ylabel, title = title)
        plt.grid()
        plt.polar([np.radians(t) for t in theta_deg], r, label = label_)
        if len(label_): plt.legend()
        
        if len(imageFile): plt.savefig(imageFile + '.png')
        
        plt.show()

    #################### OVF2 HELPERS

    #write an OVF2 file for a mesh with given rectangle (m), number of nodes and values in vec list ordered by x, then y, finally z.
    #rect_m must be a list with 6 elements : [xmin, ymin, zmin, xmax, ymax, zmax]
    #nodes must be a list with 3 integers : [xnodes, ynodes, znodes]; the cellsize is determined from rect_m and nodes
    #vec can be a scalar quantity (list of floats), or a vector quantity (list of 3-element lists)
    def Write_OVF2(self, fileName, vec, nodes, rect_m):
        """Write an OVF2 file from numpy array with given integer number of cells (nodes) and mesh rectangle (m)"""
        lines = []
        
        if not isinstance(vec, np.ndarray): vec = np.asarray(vec)
        is_vectorial = (vec.shape[1] == 3)
        
        lines.append("# OOMMF OVF 2.0")
        lines.append("# Segment count: 1")
        lines.append("#")
        lines.append("# Begin: Segment")
        lines.append("# Begin: Header")
        lines.append("#")
        lines.append("# Title: " + fileName + "")
        lines.append("#")
        lines.append("# meshunit: m")
        lines.append("#")
        lines.append("# meshtype: rectangular")
        lines.append("#")
        lines.append("# xmin: " + str(rect_m[0]))
        lines.append("# ymin: " + str(rect_m[1]))
        lines.append("# zmin: " + str(rect_m[2]))
        lines.append("# xmax: " + str(rect_m[3]))
        lines.append("# ymax: " + str(rect_m[4]))
        lines.append("# zmax: " + str(rect_m[5]))
        lines.append("#")
        lines.append("# xnodes: " + str(nodes[0]))
        lines.append("# ynodes: " + str(nodes[1]))
        lines.append("# znodes: " + str(nodes[2]))
        lines.append("#")
        lines.append("# xstepsize: " + str((rect_m[3] - rect_m[0]) / nodes[0]))
        lines.append("# ystepsize: " + str((rect_m[4] - rect_m[1]) / nodes[1]))
        lines.append("# zstepsize: " + str((rect_m[5] - rect_m[2]) / nodes[2]))
        lines.append("#")
        if is_vectorial: lines.append("# valuedim: 3")
        else: lines.append("# valuedim: 1")
        lines.append("#")
        lines.append("# End: Header")
        lines.append("#")
        lines.append("# Begin: data binary 8")
               
        with open(fileName, 'wb') as f:
            
            #write header
            for line in lines: f.write(bytes(line + '\n', 'utf-8'))
            
            #insert check value and flatten
            vec = np.insert(vec, 0, 123456789012345.0)
            #write as 8-byte floats (including the check value)
            vec.tofile(f)
            
            #write termination
            f.write(bytes("# End: data " + "binary 8" + '\n', 'utf-8'))
            f.write(bytes("# End: Segment" + '\n', 'utf-8'))

            return True
        
        return False
    
    def Read_OVF2(self, fileName):
        """Return a numpy array, integer number of cells, and mesh rectangle, from OVF2 file"""
        vec = np.ndarray(0)
        
        meshtype_rectangular = False
        meshunit = 1.0
        valuedim = 0
        data_bytes = 0
        meshRect = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        n = [0, 0, 0]
        h = [0.0, 0.0, 0.0]
        
        #return 0 if something went wrong : abort loading OVF file
    	#return 1 if everything fine, continue
    	#return 2 if start of data header found
        def scan_line(line):
            
            nonlocal meshtype_rectangular
            nonlocal meshunit
            nonlocal valuedim
            nonlocal data_bytes
            nonlocal meshRect
            nonlocal n
            nonlocal h
            
            if line.find('# meshtype: ') != -1:
                value = line[len('# meshtype: '):]
                if value != "rectangular": return 0
                meshtype_rectangular = True
                
            elif line.find('# meshunit: ') != -1:
                value = line[len('# meshunit: '):]
                if value == "m": meshunit = 1.0
                elif value == "nm": meshunit = 1e-9
                else: return 0
            
            elif line.find('# valuedim: ') != -1:
                value = line[len('# valuedim: '):]
                if value == "1": valuedim = 1
                elif value == "3": valuedim = 3
                else: return 0
                
            elif line.find('# xmin: ') != -1:
                value = line[len('# xmin: '):]
                meshRect[0] = float(value) * meshunit
                
            elif line.find('# ymin: ') != -1:
                value = line[len('# ymin: '):]
                meshRect[1] = float(value) * meshunit
                
            elif line.find('# zmin: ') != -1:
                value = line[len('# zmin: '):]
                meshRect[2] = float(value) * meshunit
                
            elif line.find('# xmax: ') != -1:
                value = line[len('# xmax: '):]
                meshRect[3] = float(value) * meshunit
                
            elif line.find('# ymax: ') != -1:
                value = line[len('# ymax: '):]
                meshRect[4] = float(value) * meshunit
                
            elif line.find('# zmax: ') != -1:
                value = line[len('# zmax: '):]
                meshRect[5] = float(value) * meshunit
                
            elif line.find('# xnodes: ') != -1:
                value = line[len('# xnodes: '):]
                n[0] = int(value)
                
            elif line.find('# ynodes: ') != -1:
                value = line[len('# ynodes: '):]
                n[1] = int(value)
                
            elif line.find('# znodes: ') != -1:
                value = line[len('# znodes: '):]
                n[2] = int(value)
                
            elif line.find('# xstepsize: ') != -1:
                value = line[len('# xstepsize: '):]
                h[0] = float(value) * meshunit
                
            elif line.find('# ystepsize: ') != -1:
                value = line[len('# ystepsize: '):]
                h[1] = float(value) * meshunit
                
            elif line.find('# zstepsize: ') != -1:
                value = line[len('# zstepsize: '):]
                h[2] = float(value) * meshunit
                
            elif line.lower().find('# begin: data ') != -1:
                value = line.lower()[len('# begin: data '):]
                                  
                if value == 'binary 4': 
                    data_bytes = 4
                    return 2
                
                if value == 'binary 8': 
                    data_bytes = 8
                    return 2

                if value == 'text': 
                    data_bytes = 1
                    return 2
            
            return 1
        
        with open(fileName, 'rb') as f:
            
            lineidx = 0
            for rawline in f:
            
                line = str(rawline.rstrip()).strip("b'")
                
                if lineidx == 0: 
                    if line != "# OOMMF OVF 2.0": return vec, n, meshRect
                else:
                    check = scan_line(line)
                    if check == 0: return vec, n, meshRect
                    elif check == 2:
                        
                        if not meshtype_rectangular or meshRect == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] or n == [0, 0, 0] or h == [0.0, 0.0, 0.0] or data_bytes == 0:
                            return vec, n, meshRect
                        
                        if data_bytes == 4:
                            value = struct.unpack('f', f.read(4))
                            if value[0] != 1234567.0: return vec, n, meshRect
                        elif data_bytes == 8:
                            value = struct.unpack('d', f.read(8))
                            if value[0] != 123456789012345.0: return vec, n, meshRect
                        elif data_bytes == 1:
                            pass
                        else:
                            return vec, n, meshRect
                        
                        if valuedim == 1:
                            vec = np.zeros((n[0]*n[1]*n[2], 1))
                            
                            #k, j, i order is important, since itertools product iterates the outer index first, and the file has i iterated first
                            for k, j, i in product(range(n[2]), range(n[1]), range(n[0])):
                                if data_bytes == 4:
                                    value = struct.unpack('f', f.read(4))
                                    vec[i + j*n[0] + k*n[0]*n[1]] = value[0]
                                elif data_bytes == 8:
                                    value = struct.unpack('d', f.read(8))
                                    vec[i + j*n[0] + k*n[0]*n[1]] = value[0] 
                                elif data_bytes == 1:
                                    rawline = f.readline()
                                    line = str(rawline.rstrip()).strip("b'")
                                    vec[i + j*n[0] + k*n[0]*n[1]] = float(line)
                            
                        elif valuedim == 3:
                            vec = np.zeros((n[0]*n[1]*n[2], 3))

                            for k, j, i in product(range(n[2]), range(n[1]), range(n[0])):
                                if data_bytes == 4:
                                    value = struct.unpack('3f', f.read(4*3))
                                    vec[i + j*n[0] + k*n[0]*n[1]] = [value[0], value[0], value[2]]
                                elif data_bytes == 8:
                                    value = struct.unpack('3d', f.read(8*3))
                                    vec[i + j*n[0] + k*n[0]*n[1]] = [value[0], value[1], value[2]] 
                                elif data_bytes == 1:
                                    rawline = f.readline()
                                    line = str(rawline.rstrip()).strip("b'").replace("\t", " ").split(" ")
                                    vec[i + j*n[0] + k*n[0]*n[1]] = [float(line[0]), float(line[1]), float(line[2])]
                                
                        return vec, n, meshRect
                
                lineidx += 1
        
        return vec, n, meshRect
        
    #################### SPECIAL COMMANDS

    #Send run command and wait for simulation to finish : blocking call
    def Run(self):

        """Run simulation and wait for it to finish: blocking call"""
        
        #Now wait for response using a blocking socket
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        
            sock.connect((self.scriptserverip, self.scriptserverport))
            
            # Look for the response
            try:
                if self.script_verbose: print('TX : run')
                sock.sendall(bytes(self.scriptserverpwd + '*' + "run", 'utf-8'))
                #receive response to run command
                str(sock.recv(self.maxLenMessage), 'utf-8')
                #now wait for "stopped" signal
                data = str(sock.recv(self.maxLenMessage), 'utf-8')
                if self.script_verbose: print('RX : %s' % data)
            except:
                print("SendCommand (receive): timed out.")
        
    #Save in given filename a new row containing parameters in dataList as tab-spaced characters
    #This uses the savecomment command to save in the local Boris data directory as currently configured in Boris
    def SaveDataToFile(self, fileName, dataList):
        """Save a single row of tab-spaced data, appending to file"""
        command = 'savecomment ' + fileName + ' '

        for entry in dataList: command += str(entry) + '\t'
        
        #send command with parameters (remove tab if last character - not needed)
        self.SendCommand(command.rstrip('\t'))

    #################### CONSOLE COMMANDS <-> METHODS

    #Structure of command methods (use this to generate them programatically after grabbing commands list with their USAGE from Boris):
    #Thus if you add new commands in Boris you can just run a separate script to update this module automatically so you don't have to keep track of changes
    #def name(self, param1 = '', param2 = '', ...):
    #   return self.SendCommand("name", [param1, param2, ...]) 
    
    def _2dmulticonvolution(self, status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("2dmulticonvolution", [status])
    	self.SendCommand("buffercommand", ["2dmulticonvolution", status])
    
    def addafmesh(self, name = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("addafmesh", [name, rectangle])
    	self.SendCommand("buffercommand", ["addafmesh", name, rectangle])
    
    def addameshcubic(self, name = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("addameshcubic", [name, rectangle])
    	self.SendCommand("buffercommand", ["addameshcubic", name, rectangle])
    
    def addconductor(self, name = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("addconductor", [name, rectangle])
    	self.SendCommand("buffercommand", ["addconductor", name, rectangle])
    
    def adddata(self, meshname = '', dataname = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("adddata", [meshname, dataname, rectangle])
    	self.SendCommand("buffercommand", ["adddata", meshname, dataname, rectangle])
    
    def adddipole(self, name = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("adddipole", [name, rectangle])
    	self.SendCommand("buffercommand", ["adddipole", name, rectangle])
    
    def addelectrode(self, electrode_rect = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("addelectrode", [electrode_rect])
    	self.SendCommand("buffercommand", ["addelectrode", electrode_rect])
    
    def addinsulator(self, name = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("addinsulator", [name, rectangle])
    	self.SendCommand("buffercommand", ["addinsulator", name, rectangle])
    
    def addmaterial(self, name = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("addmaterial", [name, rectangle])
    	self.SendCommand("buffercommand", ["addmaterial", name, rectangle])
    
    def addmdbentry(self, meshname = '', materialname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("addmdbentry", [meshname, materialname])
    	self.SendCommand("buffercommand", ["addmdbentry", meshname, materialname])
    
    def addmesh(self, name = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("addmesh", [name, rectangle])
    	self.SendCommand("buffercommand", ["addmesh", name, rectangle])
    
    def addmodule(self, meshname = '', handle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("addmodule", [meshname, handle])
    	self.SendCommand("buffercommand", ["addmodule", meshname, handle])
    
    def addpinneddata(self, meshname = '', dataname = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("addpinneddata", [meshname, dataname, rectangle])
    	self.SendCommand("buffercommand", ["addpinneddata", meshname, dataname, rectangle])
    
    def addrect(self, meshname = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("addrect", [meshname, rectangle])
    	self.SendCommand("buffercommand", ["addrect", meshname, rectangle])
    
    def addstage(self, meshname = '', stagetype = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("addstage", [meshname, stagetype])
    	self.SendCommand("buffercommand", ["addstage", meshname, stagetype])
    
    def adjustcamdistance(self, dZ = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("adjustcamdistance", [dZ])
    	self.SendCommand("buffercommand", ["adjustcamdistance", dZ])
    
    def ambient(self, meshname = '', ambient_temperature = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("ambient", [meshname, ambient_temperature])
    	self.SendCommand("buffercommand", ["ambient", meshname, ambient_temperature])
    
    def astepctrl(self, err_fail = '', dT_incr = '', dT_min = '', dT_max = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("astepctrl", [err_fail, dT_incr, dT_min, dT_max])
    	self.SendCommand("buffercommand", ["astepctrl", err_fail, dT_incr, dT_min, dT_max])
    
    def atomicmoment(self, meshname = '', ub_multiple = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("atomicmoment", [meshname, ub_multiple])
    	self.SendCommand("buffercommand", ["atomicmoment", meshname, ub_multiple])
    
    def averagemeshrect(self, meshname = '', rectangle = '', dp_index = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("averagemeshrect", [meshname, rectangle, dp_index])
    	self.SendCommand("buffercommand", ["averagemeshrect", meshname, rectangle, dp_index])
    
    def benchtime(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("benchtime")
    	self.SendCommand("buffercommand", ["benchtime"])
    
    def blochpreparemovingmesh(self, meshname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("blochpreparemovingmesh", [meshname])
    	self.SendCommand("buffercommand", ["blochpreparemovingmesh", meshname])
    
    def buffercommand(self, command = '', params = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("buffercommand", [command, params])
    	self.SendCommand("buffercommand", ["buffercommand", command, params])
    
    def cellsize(self, meshname = '', value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("cellsize", [meshname, value])
    	self.SendCommand("buffercommand", ["cellsize", meshname, value])
    
    def center(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("center")
    	self.SendCommand("buffercommand", ["center"])
    
    def chdir(self, directory = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("chdir", [directory])
    	self.SendCommand("buffercommand", ["chdir", directory])
    
    def checkupdates(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("checkupdates")
    	self.SendCommand("buffercommand", ["checkupdates"])
    
    def clearcommbuffer(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("clearcommbuffer")
    	self.SendCommand("buffercommand", ["clearcommbuffer"])
    
    def clearelectrodes(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("clearelectrodes")
    	self.SendCommand("buffercommand", ["clearelectrodes"])
    
    def clearequationconstants(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("clearequationconstants")
    	self.SendCommand("buffercommand", ["clearequationconstants"])
    
    def clearmovingmesh(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("clearmovingmesh")
    	self.SendCommand("buffercommand", ["clearmovingmesh"])
    
    def clearparamstemp(self, meshname = '', paramname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("clearparamstemp", [meshname, paramname])
    	self.SendCommand("buffercommand", ["clearparamstemp", meshname, paramname])
    
    def clearparamsvar(self, meshname = '', paramname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("clearparamsvar", [meshname, paramname])
    	self.SendCommand("buffercommand", ["clearparamsvar", meshname, paramname])
    
    def clearroughness(self, meshname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("clearroughness", [meshname])
    	self.SendCommand("buffercommand", ["clearroughness", meshname])
    
    def clearscreen(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("clearscreen")
    	self.SendCommand("buffercommand", ["clearscreen"])
    
    def computefields(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("computefields")
    	self.SendCommand("buffercommand", ["computefields"])
    
    def copymeshdata(self, meshname_from = '', meshname_to = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("copymeshdata", [meshname_from, meshname_to])
    	self.SendCommand("buffercommand", ["copymeshdata", meshname_from, meshname_to])
    
    def copyparams(self, meshname_from = '', meshname_to = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("copyparams", [meshname_from, meshname_to])
    	self.SendCommand("buffercommand", ["copyparams", meshname_from, meshname_to])
    
    def coupletodipoles(self, status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("coupletodipoles", [status])
    	self.SendCommand("buffercommand", ["coupletodipoles", status])
    
    def crosstie(self, meshname = '', direction = '', radius = '', thickness = '', centre = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("crosstie", [meshname, direction, radius, thickness, centre])
    	self.SendCommand("buffercommand", ["crosstie", meshname, direction, radius, thickness, centre])
    
    def cuda(self, status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("cuda", [status])
    	self.SendCommand("buffercommand", ["cuda", status])
    
    def curietemperature(self, meshname = '', curie_temperature = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("curietemperature", [meshname, curie_temperature])
    	self.SendCommand("buffercommand", ["curietemperature", meshname, curie_temperature])
    
    def data(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("data")
    	self.SendCommand("buffercommand", ["data"])
    
    def dataprecision(self, precision = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dataprecision", [precision])
    	self.SendCommand("buffercommand", ["dataprecision", precision])
    
    def default(self):
       	self.SendCommand("default")
        self.selectcudadevice(self.cudaDevice)
    
    def deldata(self, index = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("deldata", [index])
    	self.SendCommand("buffercommand", ["deldata", index])
    
    def delelectrode(self, index = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("delelectrode", [index])
    	self.SendCommand("buffercommand", ["delelectrode", index])
    
    def delequationconstant(self, name = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("delequationconstant", [name])
    	self.SendCommand("buffercommand", ["delequationconstant", name])
    
    def delmdbentry(self, materialname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("delmdbentry", [materialname])
    	self.SendCommand("buffercommand", ["delmdbentry", materialname])
    
    def delmesh(self, name = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("delmesh", [name])
    	self.SendCommand("buffercommand", ["delmesh", name])
    
    def delmodule(self, meshname = '', handle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("delmodule", [meshname, handle])
    	self.SendCommand("buffercommand", ["delmodule", meshname, handle])
    
    def delpinneddata(self, index = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("delpinneddata", [index])
    	self.SendCommand("buffercommand", ["delpinneddata", index])
    
    def delrect(self, meshname = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("delrect", [meshname, rectangle])
    	self.SendCommand("buffercommand", ["delrect", meshname, rectangle])
    
    def delstage(self, index = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("delstage", [index])
    	self.SendCommand("buffercommand", ["delstage", index])
    
    def designateground(self, electrode_index = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("designateground", [electrode_index])
    	self.SendCommand("buffercommand", ["designateground", electrode_index])
    
    def disabletransportsolver(self, status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("disabletransportsolver", [status])
    	self.SendCommand("buffercommand", ["disabletransportsolver", status])
    
    def diskbufferlines(self, lines = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("diskbufferlines", [lines])
    	self.SendCommand("buffercommand", ["diskbufferlines", lines])
    
    def display(self, meshname = '', name = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("display", [meshname, name])
    	self.SendCommand("buffercommand", ["display", meshname, name])
    
    def displaybackground(self, meshname = '', name = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("displaybackground", [meshname, name])
    	self.SendCommand("buffercommand", ["displaybackground", meshname, name])
    
    def displaydetail(self, size = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("displaydetail", [size])
    	self.SendCommand("buffercommand", ["displaydetail", size])
    
    def displaymodule(self, meshname = '', modulename = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("displaymodule", [meshname, modulename])
    	self.SendCommand("buffercommand", ["displaymodule", meshname, modulename])
    
    def displayrenderthresholds(self, thresh1 = '', thresh2 = '', thresh3 = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("displayrenderthresholds", [thresh1, thresh2, thresh3])
    	self.SendCommand("buffercommand", ["displayrenderthresholds", thresh1, thresh2, thresh3])
    
    def displaythresholds(self, minimum = '', maximum = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("displaythresholds", [minimum, maximum])
    	self.SendCommand("buffercommand", ["displaythresholds", minimum, maximum])
    
    def displaythresholdtrigger(self, trigtype = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("displaythresholdtrigger", [trigtype])
    	self.SendCommand("buffercommand", ["displaythresholdtrigger", trigtype])
    
    def displaytransparency(self, foreground = '', background = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("displaytransparency", [foreground, background])
    	self.SendCommand("buffercommand", ["displaytransparency", foreground, background])
    
    def dmcellsize(self, meshname = '', value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dmcellsize", [meshname, value])
    	self.SendCommand("buffercommand", ["dmcellsize", meshname, value])
    
    def dp_add(self, dp_source = '', value = '', dp_dest = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_add", [dp_source, value, dp_dest])
    	self.SendCommand("buffercommand", ["dp_add", dp_source, value, dp_dest])
    
    def dp_adddp(self, dp_x1 = '', dp_x2 = '', dp_dest = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_adddp", [dp_x1, dp_x2, dp_dest])
    	self.SendCommand("buffercommand", ["dp_adddp", dp_x1, dp_x2, dp_dest])
    
    def dp_anghistogram(self, meshname = '', dp_x = '', dp_y = '', cx = '', cy = '', cz = '', nx = '', ny = '', nz = '', numbins = '', min = '', max = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_anghistogram", [meshname, dp_x, dp_y, cx, cy, cz, nx, ny, nz, numbins, min, max])
    	self.SendCommand("buffercommand", ["dp_anghistogram", meshname, dp_x, dp_y, cx, cy, cz, nx, ny, nz, numbins, min, max])
    
    def dp_append(self, dp_original = '', dp_new = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_append", [dp_original, dp_new])
    	self.SendCommand("buffercommand", ["dp_append", dp_original, dp_new])
    
    def dp_calcsot(self, hm_mesh = '', fm_mesh = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_calcsot", [hm_mesh, fm_mesh])
    	self.SendCommand("buffercommand", ["dp_calcsot", hm_mesh, fm_mesh])
    
    def dp_calctopochargedensity(self, meshname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_calctopochargedensity", [meshname])
    	self.SendCommand("buffercommand", ["dp_calctopochargedensity", meshname])
    
    def dp_cartesiantopolar(self, dp_in_x = '', dp_in_y = '', dp_out_r = '', dp_out_theta = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_cartesiantopolar", [dp_in_x, dp_in_y, dp_out_r, dp_out_theta])
    	self.SendCommand("buffercommand", ["dp_cartesiantopolar", dp_in_x, dp_in_y, dp_out_r, dp_out_theta])
    
    def dp_chunkedstd(self, dp_index = '', chunk = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_chunkedstd", [dp_index, chunk])
    	self.SendCommand("buffercommand", ["dp_chunkedstd", dp_index, chunk])
    
    def dp_clear(self, indexes = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_clear", [indexes])
    	self.SendCommand("buffercommand", ["dp_clear", indexes])
    
    def dp_clearall(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_clearall")
    	self.SendCommand("buffercommand", ["dp_clearall"])
    
    def dp_coercivity(self, dp_index_x = '', dp_index_y = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_coercivity", [dp_index_x, dp_index_y])
    	self.SendCommand("buffercommand", ["dp_coercivity", dp_index_x, dp_index_y])
    
    def dp_completehysteresis(self, dp_index_x = '', dp_index_y = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_completehysteresis", [dp_index_x, dp_index_y])
    	self.SendCommand("buffercommand", ["dp_completehysteresis", dp_index_x, dp_index_y])
    
    def dp_countskyrmions(self, meshname = '', x = '', y = '', radius = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_countskyrmions", [meshname, x, y, radius])
    	self.SendCommand("buffercommand", ["dp_countskyrmions", meshname, x, y, radius])
    
    def dp_crossingsfrequency(self, dp_in_x = '', dp_in_y = '', dp_level = '', dp_freq_up = '', dp_freq_dn = '', steps = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_crossingsfrequency", [dp_in_x, dp_in_y, dp_level, dp_freq_up, dp_freq_dn, steps])
    	self.SendCommand("buffercommand", ["dp_crossingsfrequency", dp_in_x, dp_in_y, dp_level, dp_freq_up, dp_freq_dn, steps])
    
    def dp_crossingshistogram(self, dp_in_x = '', dp_in_y = '', dp_level = '', dp_counts = '', steps = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_crossingshistogram", [dp_in_x, dp_in_y, dp_level, dp_counts, steps])
    	self.SendCommand("buffercommand", ["dp_crossingshistogram", dp_in_x, dp_in_y, dp_level, dp_counts, steps])
    
    def dp_div(self, dp_source = '', value = '', dp_dest = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_div", [dp_source, value, dp_dest])
    	self.SendCommand("buffercommand", ["dp_div", dp_source, value, dp_dest])
    
    def dp_divdp(self, dp_x1 = '', dp_x2 = '', dp_dest = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_divdp", [dp_x1, dp_x2, dp_dest])
    	self.SendCommand("buffercommand", ["dp_divdp", dp_x1, dp_x2, dp_dest])
    
    def dp_dotprod(self, dp_vector = '', ux = '', uy = '', uz = '', dp_out = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_dotprod", [dp_vector, ux, uy, uz, dp_out])
    	self.SendCommand("buffercommand", ["dp_dotprod", dp_vector, ux, uy, uz, dp_out])
    
    def dp_dotproddp(self, dp_x1 = '', dp_x2 = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_dotproddp", [dp_x1, dp_x2])
    	self.SendCommand("buffercommand", ["dp_dotproddp", dp_x1, dp_x2])
    
    def dp_dumptdep(self, meshname = '', paramname = '', max_temperature = '', dp_index = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_dumptdep", [meshname, paramname, max_temperature, dp_index])
    	self.SendCommand("buffercommand", ["dp_dumptdep", meshname, paramname, max_temperature, dp_index])
    
    def dp_erase(self, dp_index = '', start_index = '', length = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_erase", [dp_index, start_index, length])
    	self.SendCommand("buffercommand", ["dp_erase", dp_index, start_index, length])
    
    def dp_extract(self, dp_in = '', dp_out = '', start_index = '', length = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_extract", [dp_in, dp_out, start_index, length])
    	self.SendCommand("buffercommand", ["dp_extract", dp_in, dp_out, start_index, length])
    
    def dp_fitadiabatic(self, meshname = '', abs_err = '', Rsq = '', T_ratio = '', stencil = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_fitadiabatic", [meshname, abs_err, Rsq, T_ratio, stencil])
    	self.SendCommand("buffercommand", ["dp_fitadiabatic", meshname, abs_err, Rsq, T_ratio, stencil])
    
    def dp_fitdw(self, dp_x = '', dp_y = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_fitdw", [dp_x, dp_y])
    	self.SendCommand("buffercommand", ["dp_fitdw", dp_x, dp_y])
    
    def dp_fitlorentz(self, dp_x = '', dp_y = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_fitlorentz", [dp_x, dp_y])
    	self.SendCommand("buffercommand", ["dp_fitlorentz", dp_x, dp_y])
    
    def dp_fitlorentz2(self, dp_x = '', dp_y = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_fitlorentz2", [dp_x, dp_y])
    	self.SendCommand("buffercommand", ["dp_fitlorentz2", dp_x, dp_y])
    
    def dp_fitnonadiabatic(self, meshname = '', abs_err = '', Rsq = '', T_ratio = '', stencil = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_fitnonadiabatic", [meshname, abs_err, Rsq, T_ratio, stencil])
    	self.SendCommand("buffercommand", ["dp_fitnonadiabatic", meshname, abs_err, Rsq, T_ratio, stencil])
    
    def dp_fitskyrmion(self, dp_x = '', dp_y = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_fitskyrmion", [dp_x, dp_y])
    	self.SendCommand("buffercommand", ["dp_fitskyrmion", dp_x, dp_y])
    
    def dp_fitsot(self, meshname = '', hm_mesh = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_fitsot", [meshname, hm_mesh, rectangle])
    	self.SendCommand("buffercommand", ["dp_fitsot", meshname, hm_mesh, rectangle])
    
    def dp_fitsotstt(self, meshname = '', hm_mesh = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_fitsotstt", [meshname, hm_mesh, rectangle])
    	self.SendCommand("buffercommand", ["dp_fitsotstt", meshname, hm_mesh, rectangle])
    
    def dp_fitstt(self, meshname = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_fitstt", [meshname, rectangle])
    	self.SendCommand("buffercommand", ["dp_fitstt", meshname, rectangle])
    
    def dp_get(self, dp_arr = '', index = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_get", [dp_arr, index])
    	self.SendCommand("buffercommand", ["dp_get", dp_arr, index])
    
    def dp_getampli(self, dp_source = '', pointsPeriod = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_getampli", [dp_source, pointsPeriod])
    	self.SendCommand("buffercommand", ["dp_getampli", dp_source, pointsPeriod])
    
    def dp_getaveragedprofile(self, meshname = '', start = '', end = '', step = '', dp_index = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_getaveragedprofile", [meshname, start, end, step, dp_index])
    	self.SendCommand("buffercommand", ["dp_getaveragedprofile", meshname, start, end, step, dp_index])
    
    def dp_getexactprofile(self, meshname = '', start = '', end = '', step = '', dp_index = '', stencil = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_getexactprofile", [meshname, start, end, step, dp_index, stencil])
    	self.SendCommand("buffercommand", ["dp_getexactprofile", meshname, start, end, step, dp_index, stencil])
    
    def dp_getprofile(self, start = '', end = '', dp_index = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_getprofile", [start, end, dp_index])
    	self.SendCommand("buffercommand", ["dp_getprofile", start, end, dp_index])
    
    def dp_histogram(self, meshname = '', dp_x = '', dp_y = '', cx = '', cy = '', cz = '', numbins = '', min = '', max = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_histogram", [meshname, dp_x, dp_y, cx, cy, cz, numbins, min, max])
    	self.SendCommand("buffercommand", ["dp_histogram", meshname, dp_x, dp_y, cx, cy, cz, numbins, min, max])
    
    def dp_histogram2(self, meshname = '', dp_x = '', dp_y = '', numbins = '', min = '', max = '', M2 = '', deltaM2 = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_histogram2", [meshname, dp_x, dp_y, numbins, min, max, M2, deltaM2])
    	self.SendCommand("buffercommand", ["dp_histogram2", meshname, dp_x, dp_y, numbins, min, max, M2, deltaM2])
    
    def dp_linreg(self, dp_index_x = '', dp_index_y = '', dp_index_z = '', dp_index_out = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_linreg", [dp_index_x, dp_index_y, dp_index_z, dp_index_out])
    	self.SendCommand("buffercommand", ["dp_linreg", dp_index_x, dp_index_y, dp_index_z, dp_index_out])
    
    def dp_load(self, filename = '', file_indexes = '', dp_indexes = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_load", [filename, file_indexes, dp_indexes])
    	self.SendCommand("buffercommand", ["dp_load", filename, file_indexes, dp_indexes])
    
    def dp_mean(self, dp_index = '', exclusion_ratio = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_mean", [dp_index, exclusion_ratio])
    	self.SendCommand("buffercommand", ["dp_mean", dp_index, exclusion_ratio])
    
    def dp_minmax(self, dp_index = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_minmax", [dp_index])
    	self.SendCommand("buffercommand", ["dp_minmax", dp_index])
    
    def dp_monotonic(self, dp_in_x = '', dp_in_y = '', dp_out_x = '', dp_out_y = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_monotonic", [dp_in_x, dp_in_y, dp_out_x, dp_out_y])
    	self.SendCommand("buffercommand", ["dp_monotonic", dp_in_x, dp_in_y, dp_out_x, dp_out_y])
    
    def dp_mul(self, dp_source = '', value = '', dp_dest = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_mul", [dp_source, value, dp_dest])
    	self.SendCommand("buffercommand", ["dp_mul", dp_source, value, dp_dest])
    
    def dp_muldp(self, dp_x1 = '', dp_x2 = '', dp_dest = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_muldp", [dp_x1, dp_x2, dp_dest])
    	self.SendCommand("buffercommand", ["dp_muldp", dp_x1, dp_x2, dp_dest])
    
    def dp_newfile(self, filename = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_newfile", [filename])
    	self.SendCommand("buffercommand", ["dp_newfile", filename])
    
    def dp_peaksfrequency(self, dp_in_x = '', dp_in_y = '', dp_level = '', dp_freq = '', steps = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_peaksfrequency", [dp_in_x, dp_in_y, dp_level, dp_freq, steps])
    	self.SendCommand("buffercommand", ["dp_peaksfrequency", dp_in_x, dp_in_y, dp_level, dp_freq, steps])
    
    def dp_pow(self, dp_source = '', exponent = '', dp_dest = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_pow", [dp_source, exponent, dp_dest])
    	self.SendCommand("buffercommand", ["dp_pow", dp_source, exponent, dp_dest])
    
    def dp_rarefy(self, dp_in = '', dp_out = '', skip = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_rarefy", [dp_in, dp_out, skip])
    	self.SendCommand("buffercommand", ["dp_rarefy", dp_in, dp_out, skip])
    
    def dp_remanence(self, dp_index_x = '', dp_index_y = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_remanence", [dp_index_x, dp_index_y])
    	self.SendCommand("buffercommand", ["dp_remanence", dp_index_x, dp_index_y])
    
    def dp_removeoffset(self, dp_index = '', dp_index_out = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_removeoffset", [dp_index, dp_index_out])
    	self.SendCommand("buffercommand", ["dp_removeoffset", dp_index, dp_index_out])
    
    def dp_replacerepeats(self, dp_index = '', dp_index_out = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_replacerepeats", [dp_index, dp_index_out])
    	self.SendCommand("buffercommand", ["dp_replacerepeats", dp_index, dp_index_out])
    
    def dp_save(self, filename = '', dp_indexes = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_save", [filename, dp_indexes])
    	self.SendCommand("buffercommand", ["dp_save", filename, dp_indexes])
    
    def dp_saveappend(self, filename = '', dp_indexes = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_saveappend", [filename, dp_indexes])
    	self.SendCommand("buffercommand", ["dp_saveappend", filename, dp_indexes])
    
    def dp_saveappendasrow(self, filename = '', dp_index = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_saveappendasrow", [filename, dp_index])
    	self.SendCommand("buffercommand", ["dp_saveappendasrow", filename, dp_index])
    
    def dp_saveasrow(self, filename = '', dp_index = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_saveasrow", [filename, dp_index])
    	self.SendCommand("buffercommand", ["dp_saveasrow", filename, dp_index])
    
    def dp_sequence(self, dp_index = '', start_value = '', increment = '', points = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_sequence", [dp_index, start_value, increment, points])
    	self.SendCommand("buffercommand", ["dp_sequence", dp_index, start_value, increment, points])
    
    def dp_set(self, dp_arr = '', index = '', value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_set", [dp_arr, index, value])
    	self.SendCommand("buffercommand", ["dp_set", dp_arr, index, value])
    
    def dp_showsizes(self, dp_arr = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_showsizes", [dp_arr])
    	self.SendCommand("buffercommand", ["dp_showsizes", dp_arr])
    
    def dp_smooth(self, dp_in = '', dp_out = '', window_size = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_smooth", [dp_in, dp_out, window_size])
    	self.SendCommand("buffercommand", ["dp_smooth", dp_in, dp_out, window_size])
    
    def dp_sub(self, dp_source = '', value = '', dp_dest = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_sub", [dp_source, value, dp_dest])
    	self.SendCommand("buffercommand", ["dp_sub", dp_source, value, dp_dest])
    
    def dp_subdp(self, dp_x1 = '', dp_x2 = '', dp_dest = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_subdp", [dp_x1, dp_x2, dp_dest])
    	self.SendCommand("buffercommand", ["dp_subdp", dp_x1, dp_x2, dp_dest])
    
    def dp_sum(self, dp_index = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_sum", [dp_index])
    	self.SendCommand("buffercommand", ["dp_sum", dp_index])
    
    def dp_thavanghistogram(self, meshname = '', dp_x = '', dp_y = '', cx = '', cy = '', cz = '', nx = '', ny = '', nz = '', numbins = '', min = '', max = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_thavanghistogram", [meshname, dp_x, dp_y, cx, cy, cz, nx, ny, nz, numbins, min, max])
    	self.SendCommand("buffercommand", ["dp_thavanghistogram", meshname, dp_x, dp_y, cx, cy, cz, nx, ny, nz, numbins, min, max])
    
    def dp_thavhistogram(self, meshname = '', dp_x = '', dp_y = '', cx = '', cy = '', cz = '', numbins = '', min = '', max = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_thavhistogram", [meshname, dp_x, dp_y, cx, cy, cz, numbins, min, max])
    	self.SendCommand("buffercommand", ["dp_thavhistogram", meshname, dp_x, dp_y, cx, cy, cz, numbins, min, max])
    
    def dp_topocharge(self, meshname = '', x = '', y = '', radius = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dp_topocharge", [meshname, x, y, radius])
    	self.SendCommand("buffercommand", ["dp_topocharge", meshname, x, y, radius])
    
    def dwall(self, meshname = '', longitudinal = '', transverse = '', width = '', position = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dwall", [meshname, longitudinal, transverse, width, position])
    	self.SendCommand("buffercommand", ["dwall", meshname, longitudinal, transverse, width, position])
    
    def dwposcomponent(self, value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("dwposcomponent", [value])
    	self.SendCommand("buffercommand", ["dwposcomponent", value])
    
    def ecellsize(self, meshname = '', value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("ecellsize", [meshname, value])
    	self.SendCommand("buffercommand", ["ecellsize", meshname, value])
    
    def editdata(self, index = '', meshname = '', dataname = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("editdata", [index, meshname, dataname, rectangle])
    	self.SendCommand("buffercommand", ["editdata", index, meshname, dataname, rectangle])
    
    def editdatasave(self, index = '', savetype = '', savevalue = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("editdatasave", [index, savetype, savevalue])
    	self.SendCommand("buffercommand", ["editdatasave", index, savetype, savevalue])
    
    def editstage(self, index = '', meshname = '', stagetype = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("editstage", [index, meshname, stagetype])
    	self.SendCommand("buffercommand", ["editstage", index, meshname, stagetype])
    
    def editstagestop(self, index = '', stoptype = '', stopvalue = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("editstagestop", [index, stoptype, stopvalue])
    	self.SendCommand("buffercommand", ["editstagestop", index, stoptype, stopvalue])
    
    def editstagevalue(self, index = '', value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("editstagevalue", [index, value])
    	self.SendCommand("buffercommand", ["editstagevalue", index, value])
    
    def electrodes(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("electrodes")
    	self.SendCommand("buffercommand", ["electrodes"])
    
    def equationconstants(self, name = '', value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("equationconstants", [name, value])
    	self.SendCommand("buffercommand", ["equationconstants", name, value])
    
    def errorlog(self, status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("errorlog", [status])
    	self.SendCommand("buffercommand", ["errorlog", status])
    
    def escellsize(self, value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("escellsize", [value])
    	self.SendCommand("buffercommand", ["escellsize", value])
    
    def evalspeedup(self, level = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("evalspeedup", [level])
    	self.SendCommand("buffercommand", ["evalspeedup", level])
    
    def exchangecoupledmeshes(self, meshname = '', status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("exchangecoupledmeshes", [meshname, status])
    	self.SendCommand("buffercommand", ["exchangecoupledmeshes", meshname, status])
    
    def excludemulticonvdemag(self, meshname = '', status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("excludemulticonvdemag", [meshname, status])
    	self.SendCommand("buffercommand", ["excludemulticonvdemag", meshname, status])
    
    def flower(self, meshname = '', direction = '', radius = '', thickness = '', centre = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("flower", [meshname, direction, radius, thickness, centre])
    	self.SendCommand("buffercommand", ["flower", meshname, direction, radius, thickness, centre])
    
    def flusherrorlog(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("flusherrorlog")
    	self.SendCommand("buffercommand", ["flusherrorlog"])
    
    def fmscellsize(self, value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("fmscellsize", [value])
    	self.SendCommand("buffercommand", ["fmscellsize", value])
    
    def generate2dgrains(self, meshname = '', spacing = '', seed = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("generate2dgrains", [meshname, spacing, seed])
    	self.SendCommand("buffercommand", ["generate2dgrains", meshname, spacing, seed])
    
    def generate3dgrains(self, meshname = '', spacing = '', seed = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("generate3dgrains", [meshname, spacing, seed])
    	self.SendCommand("buffercommand", ["generate3dgrains", meshname, spacing, seed])
    
    def gpukernels(self, status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("gpukernels", [status])
    	self.SendCommand("buffercommand", ["gpukernels", status])
    
    def imagecropping(self, left = '', bottom = '', right = '', top = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("imagecropping", [left, bottom, right, top])
    	self.SendCommand("buffercommand", ["imagecropping", left, bottom, right, top])
    
    def individualmaskshape(self, status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("individualmaskshape", [status])
    	self.SendCommand("buffercommand", ["individualmaskshape", status])
    
    def insulatingside(self, meshname = '', side_literal = '', status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("insulatingside", [meshname, side_literal, status])
    	self.SendCommand("buffercommand", ["insulatingside", meshname, side_literal, status])
    
    def invertmag(self, meshname = '', components = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("invertmag", [meshname, components])
    	self.SendCommand("buffercommand", ["invertmag", meshname, components])
    
    def isrunning(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("isrunning")
    	self.SendCommand("buffercommand", ["isrunning"])
    
    def iterupdate(self, iterations = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("iterupdate", [iterations])
    	self.SendCommand("buffercommand", ["iterupdate", iterations])
    
    def linkdtspeedup(self, flag = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("linkdtspeedup", [flag])
    	self.SendCommand("buffercommand", ["linkdtspeedup", flag])
    
    def linkdtstochastic(self, flag = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("linkdtstochastic", [flag])
    	self.SendCommand("buffercommand", ["linkdtstochastic", flag])
    
    def linkstochastic(self, meshname = '', flag = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("linkstochastic", [meshname, flag])
    	self.SendCommand("buffercommand", ["linkstochastic", meshname, flag])
    
    def loadmaskfile(self, meshname = '', z_depth = '', filename = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("loadmaskfile", [meshname, z_depth, filename])
    	self.SendCommand("buffercommand", ["loadmaskfile", meshname, z_depth, filename])
    
    def loadovf2curr(self, meshname = '', filename = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("loadovf2curr", [meshname, filename])
    	self.SendCommand("buffercommand", ["loadovf2curr", meshname, filename])
    
    def loadovf2disp(self, meshname = '', filename = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("loadovf2disp", [meshname, filename])
    	self.SendCommand("buffercommand", ["loadovf2disp", meshname, filename])
    
    def loadovf2field(self, meshname = '', filename = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("loadovf2field", [meshname, filename])
    	self.SendCommand("buffercommand", ["loadovf2field", meshname, filename])
    
    def loadovf2mag(self, meshname = '', renormalize_value = '', filename = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("loadovf2mag", [meshname, renormalize_value, filename])
    	self.SendCommand("buffercommand", ["loadovf2mag", meshname, renormalize_value, filename])
    
    def loadovf2strain(self, meshname = '', filename_diag = '', filename_odiag = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("loadovf2strain", [meshname, filename_diag, filename_odiag])
    	self.SendCommand("buffercommand", ["loadovf2strain", meshname, filename_diag, filename_odiag])
    
    def loadovf2temp(self, meshname = '', filename = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("loadovf2temp", [meshname, filename])
    	self.SendCommand("buffercommand", ["loadovf2temp", meshname, filename])
    
    def loadsim(self, filename = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("loadsim", [filename])
    	self.SendCommand("buffercommand", ["loadsim", filename])
    
    def makevideo(self, filebase = '', fps = '', quality = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("makevideo", [filebase, fps, quality])
    	self.SendCommand("buffercommand", ["makevideo", filebase, fps, quality])
    
    def manual(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("manual")
    	self.SendCommand("buffercommand", ["manual"])
    
    def matcurietemperature(self, meshname = '', curie_temperature = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("matcurietemperature", [meshname, curie_temperature])
    	self.SendCommand("buffercommand", ["matcurietemperature", meshname, curie_temperature])
    
    def materialsdatabase(self, mdbname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("materialsdatabase", [mdbname])
    	self.SendCommand("buffercommand", ["materialsdatabase", mdbname])
    
    def mccomputefields(self, status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("mccomputefields", [status])
    	self.SendCommand("buffercommand", ["mccomputefields", status])
    
    def mcconeangle(self, min_angle = '', max_angle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("mcconeangle", [min_angle, max_angle])
    	self.SendCommand("buffercommand", ["mcconeangle", min_angle, max_angle])
    
    def mcconstrain(self, meshname = '', value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("mcconstrain", [meshname, value])
    	self.SendCommand("buffercommand", ["mcconstrain", meshname, value])
    
    def mcdisable(self, meshname = '', status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("mcdisable", [meshname, status])
    	self.SendCommand("buffercommand", ["mcdisable", meshname, status])
    
    def mcellsize(self, meshname = '', value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("mcellsize", [meshname, value])
    	self.SendCommand("buffercommand", ["mcellsize", meshname, value])
    
    def mcserial(self, meshname = '', value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("mcserial", [meshname, value])
    	self.SendCommand("buffercommand", ["mcserial", meshname, value])
    
    def memory(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("memory")
    	self.SendCommand("buffercommand", ["memory"])
    
    def mesh(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("mesh")
    	self.SendCommand("buffercommand", ["mesh"])
    
    def meshfocus(self, meshname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("meshfocus", [meshname])
    	self.SendCommand("buffercommand", ["meshfocus", meshname])
    
    def meshfocus2(self, meshname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("meshfocus2", [meshname])
    	self.SendCommand("buffercommand", ["meshfocus2", meshname])
    
    def meshrect(self, meshname = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("meshrect", [meshname, rectangle])
    	self.SendCommand("buffercommand", ["meshrect", meshname, rectangle])
    
    def mirrormag(self, meshname = '', axis = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("mirrormag", [meshname, axis])
    	self.SendCommand("buffercommand", ["mirrormag", meshname, axis])
    
    def modules(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("modules")
    	self.SendCommand("buffercommand", ["modules"])
    
    def movingmesh(self, status_or_meshname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("movingmesh", [status_or_meshname])
    	self.SendCommand("buffercommand", ["movingmesh", status_or_meshname])
    
    def movingmeshasym(self, status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("movingmeshasym", [status])
    	self.SendCommand("buffercommand", ["movingmeshasym", status])
    
    def movingmeshthresh(self, value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("movingmeshthresh", [value])
    	self.SendCommand("buffercommand", ["movingmeshthresh", value])
    
    def multiconvolution(self, status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("multiconvolution", [status])
    	self.SendCommand("buffercommand", ["multiconvolution", status])
    
    def ncommon(self, sizes = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("ncommon", [sizes])
    	self.SendCommand("buffercommand", ["ncommon", sizes])
    
    def ncommonstatus(self, status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("ncommonstatus", [status])
    	self.SendCommand("buffercommand", ["ncommonstatus", status])
    
    def neelpreparemovingmesh(self, meshname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("neelpreparemovingmesh", [meshname])
    	self.SendCommand("buffercommand", ["neelpreparemovingmesh", meshname])
    
    def newinstance(self, port = '', cudaDevice = '', password = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("newinstance", [port, cudaDevice, password])
    	self.SendCommand("buffercommand", ["newinstance", port, cudaDevice, password])
    
    def ode(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("ode")
    	self.SendCommand("buffercommand", ["ode"])
    
    def onion(self, meshname = '', direction = '', radius1 = '', radius2 = '', thickness = '', centre = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("onion", [meshname, direction, radius1, radius2, thickness, centre])
    	self.SendCommand("buffercommand", ["onion", meshname, direction, radius1, radius2, thickness, centre])
    
    def params(self, meshname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("params", [meshname])
    	self.SendCommand("buffercommand", ["params", meshname])
    
    def paramstemp(self, meshname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("paramstemp", [meshname])
    	self.SendCommand("buffercommand", ["paramstemp", meshname])
    
    def paramsvar(self, meshname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("paramsvar", [meshname])
    	self.SendCommand("buffercommand", ["paramsvar", meshname])
    
    def pbc(self, meshname = '', flag = '', images = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("pbc", [meshname, flag, images])
    	self.SendCommand("buffercommand", ["pbc", meshname, flag, images])
    
    def preparemovingmesh(self, meshname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("preparemovingmesh", [meshname])
    	self.SendCommand("buffercommand", ["preparemovingmesh", meshname])
    
    def random(self, meshname = '', seed = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("random", [meshname, seed])
    	self.SendCommand("buffercommand", ["random", meshname, seed])
    
    def randomxy(self, meshname = '', seed = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("randomxy", [meshname, seed])
    	self.SendCommand("buffercommand", ["randomxy", meshname, seed])
    
    def refineroughness(self, meshname = '', value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("refineroughness", [meshname, value])
    	self.SendCommand("buffercommand", ["refineroughness", meshname, value])
    
    def refreshmdb(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("refreshmdb")
    	self.SendCommand("buffercommand", ["refreshmdb"])
    
    def refreshscreen(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("refreshscreen")
    	self.SendCommand("buffercommand", ["refreshscreen"])
    
    def renamemesh(self, old_name = '', new_name = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("renamemesh", [old_name, new_name])
    	self.SendCommand("buffercommand", ["renamemesh", old_name, new_name])
    
    def requestmdbsync(self, materialname = '', email = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("requestmdbsync", [materialname, email])
    	self.SendCommand("buffercommand", ["requestmdbsync", materialname, email])
    
    def reset(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("reset")
    	self.SendCommand("buffercommand", ["reset"])
    
    def resetmesh(self, meshname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("resetmesh", [meshname])
    	self.SendCommand("buffercommand", ["resetmesh", meshname])
    
    def robinalpha(self, meshname = '', robin_alpha = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("robinalpha", [meshname, robin_alpha])
    	self.SendCommand("buffercommand", ["robinalpha", meshname, robin_alpha])
    
    def rotcamaboutaxis(self, dAngle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("rotcamaboutaxis", [dAngle])
    	self.SendCommand("buffercommand", ["rotcamaboutaxis", dAngle])
    
    def rotcamaboutorigin(self, dAzim = '', dPolar = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("rotcamaboutorigin", [dAzim, dPolar])
    	self.SendCommand("buffercommand", ["rotcamaboutorigin", dAzim, dPolar])
    
    def roughenmesh(self, meshname = '', depth = '', side = '', seed = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("roughenmesh", [meshname, depth, side, seed])
    	self.SendCommand("buffercommand", ["roughenmesh", meshname, depth, side, seed])
    
    def runcommbuffer(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("runcommbuffer")
    	self.SendCommand("buffercommand", ["runcommbuffer"])
    
    def savecomment(self, filename = '', comment = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("savecomment", [filename, comment])
    	self.SendCommand("buffercommand", ["savecomment", filename, comment])
    
    def savedatafile(self, filename = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("savedatafile", [filename])
    	self.SendCommand("buffercommand", ["savedatafile", filename])
    
    def savedataflag(self, status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("savedataflag", [status])
    	self.SendCommand("buffercommand", ["savedataflag", status])
    
    def saveimage(self, filename = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("saveimage", [filename])
    	self.SendCommand("buffercommand", ["saveimage", filename])
    
    def saveimagefile(self, filename = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("saveimagefile", [filename])
    	self.SendCommand("buffercommand", ["saveimagefile", filename])
    
    def saveimageflag(self, status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("saveimageflag", [status])
    	self.SendCommand("buffercommand", ["saveimageflag", status])
    
    def savemeshimage(self, filename = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("savemeshimage", [filename])
    	self.SendCommand("buffercommand", ["savemeshimage", filename])
    
    def saveovf2(self, meshname = '', data_type = '', filename = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("saveovf2", [meshname, data_type, filename])
    	self.SendCommand("buffercommand", ["saveovf2", meshname, data_type, filename])
    
    def saveovf2mag(self, meshname = '', n = '', data_type = '', filename = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("saveovf2mag", [meshname, n, data_type, filename])
    	self.SendCommand("buffercommand", ["saveovf2mag", meshname, n, data_type, filename])
    
    def saveovf2param(self, meshname = '', data_type = '', paramname = '', filename = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("saveovf2param", [meshname, data_type, paramname, filename])
    	self.SendCommand("buffercommand", ["saveovf2param", meshname, data_type, paramname, filename])
    
    def savesim(self, filename = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("savesim", [filename])
    	self.SendCommand("buffercommand", ["savesim", filename])
    
    def scalemeshrects(self, status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("scalemeshrects", [status])
    	self.SendCommand("buffercommand", ["scalemeshrects", status])
    
    def scellsize(self, meshname = '', value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("scellsize", [meshname, value])
    	self.SendCommand("buffercommand", ["scellsize", meshname, value])
    
    def scriptserver(self, status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("scriptserver", [status])
    	self.SendCommand("buffercommand", ["scriptserver", status])
    
    def selectcudadevice(self, number = '', bufferCommand = False):
    	if not bufferCommand: 
            self.cudaDevice = number
            return self.SendCommand("selectcudadevice", [number])
    	self.SendCommand("buffercommand", ["selectcudadevice", number])
    
    def serverpassword(self, password = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("serverpassword", [password])
    	self.SendCommand("buffercommand", ["serverpassword", password])
    
    def serverport(self, port = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("serverport", [port])
    	self.SendCommand("buffercommand", ["serverport", port])
    
    def serversleepms(self, time_ms = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("serversleepms", [time_ms])
    	self.SendCommand("buffercommand", ["serversleepms", time_ms])
    
    def setafmesh(self, name = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setafmesh", [name, rectangle])
    	self.SendCommand("buffercommand", ["setafmesh", name, rectangle])
    
    def setameshcubic(self, name = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setameshcubic", [name, rectangle])
    	self.SendCommand("buffercommand", ["setameshcubic", name, rectangle])
    
    def setangle(self, meshname = '', polar = '', azimuthal = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setangle", [meshname, polar, azimuthal])
    	self.SendCommand("buffercommand", ["setangle", meshname, polar, azimuthal])
    
    def setatomode(self, equation = '', evaluation = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setatomode", [equation, evaluation])
    	self.SendCommand("buffercommand", ["setatomode", equation, evaluation])
    
    def setcurrent(self, current = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setcurrent", [current])
    	self.SendCommand("buffercommand", ["setcurrent", current])
    
    def setcurrentdensity(self, meshname = '', Jx = '', Jy = '', Jz = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setcurrentdensity", [meshname, Jx, Jy, Jz])
    	self.SendCommand("buffercommand", ["setcurrentdensity", meshname, Jx, Jy, Jz])
    
    def setdata(self, meshname = '', dataname = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setdata", [meshname, dataname, rectangle])
    	self.SendCommand("buffercommand", ["setdata", meshname, dataname, rectangle])
    
    def setdefaultelectrodes(self, sides = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setdefaultelectrodes", [sides])
    	self.SendCommand("buffercommand", ["setdefaultelectrodes", sides])
    
    def setdisplayedparamsvar(self, meshname = '', paramname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setdisplayedparamsvar", [meshname, paramname])
    	self.SendCommand("buffercommand", ["setdisplayedparamsvar", meshname, paramname])
    
    def setdt(self, value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setdt", [value])
    	self.SendCommand("buffercommand", ["setdt", value])
    
    def setdtspeedup(self, value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setdtspeedup", [value])
    	self.SendCommand("buffercommand", ["setdtspeedup", value])
    
    def setdtstoch(self, value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setdtstoch", [value])
    	self.SendCommand("buffercommand", ["setdtstoch", value])
    
    def setelectrodepotential(self, electrode_index = '', potential = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setelectrodepotential", [electrode_index, potential])
    	self.SendCommand("buffercommand", ["setelectrodepotential", electrode_index, potential])
    
    def setelectroderect(self, electrode_index = '', electrode_rect = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setelectroderect", [electrode_index, electrode_rect])
    	self.SendCommand("buffercommand", ["setelectroderect", electrode_index, electrode_rect])
    
    def setfield(self, meshname = '', magnitude = '', polar = '', azimuthal = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setfield", [meshname, magnitude, polar, azimuthal])
    	self.SendCommand("buffercommand", ["setfield", meshname, magnitude, polar, azimuthal])
    
    def setheatdt(self, value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setheatdt", [value])
    	self.SendCommand("buffercommand", ["setheatdt", value])
    
    def setktens(self, term1 = '', term2 = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setktens", [term1, term2])
    	self.SendCommand("buffercommand", ["setktens", term1, term2])
    
    def setmaterial(self, name = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setmaterial", [name, rectangle])
    	self.SendCommand("buffercommand", ["setmaterial", name, rectangle])
    
    def setmesh(self, name = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setmesh", [name, rectangle])
    	self.SendCommand("buffercommand", ["setmesh", name, rectangle])
    
    def setobjectangle(self, meshname = '', polar = '', azimuthal = '', position = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setobjectangle", [meshname, polar, azimuthal, position])
    	self.SendCommand("buffercommand", ["setobjectangle", meshname, polar, azimuthal, position])
    
    def setode(self, equation = '', evaluation = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setode", [equation, evaluation])
    	self.SendCommand("buffercommand", ["setode", equation, evaluation])
    
    def setodeeval(self, evaluation = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setodeeval", [evaluation])
    	self.SendCommand("buffercommand", ["setodeeval", evaluation])
    
    def setparam(self, meshname = '', paramname = '', value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setparam", [meshname, paramname, value])
    	self.SendCommand("buffercommand", ["setparam", meshname, paramname, value])
    
    def setparamtemparray(self, meshname = '', paramname = '', filename = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setparamtemparray", [meshname, paramname, filename])
    	self.SendCommand("buffercommand", ["setparamtemparray", meshname, paramname, filename])
    
    def setparamtempequation(self, meshname = '', paramname = '', text_equation = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setparamtempequation", [meshname, paramname, text_equation])
    	self.SendCommand("buffercommand", ["setparamtempequation", meshname, paramname, text_equation])
    
    def setparamvar(self, meshname = '', paramname = '', generatorname = '', arguments = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setparamvar", [meshname, paramname, generatorname, arguments])
    	self.SendCommand("buffercommand", ["setparamvar", meshname, paramname, generatorname, arguments])
    
    def setpotential(self, potential = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setpotential", [potential])
    	self.SendCommand("buffercommand", ["setpotential", potential])
    
    def setrect(self, meshname = '', polar = '', azimuthal = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setrect", [meshname, polar, azimuthal, rectangle])
    	self.SendCommand("buffercommand", ["setrect", meshname, polar, azimuthal, rectangle])
    
    def setsordamping(self, damping_v = '', damping_s = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setsordamping", [damping_v, damping_s])
    	self.SendCommand("buffercommand", ["setsordamping", damping_v, damping_s])
    
    def setstage(self, meshname = '', stagetype = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setstage", [meshname, stagetype])
    	self.SendCommand("buffercommand", ["setstage", meshname, stagetype])
    
    def setstress(self, meshname = '', magnitude = '', polar = '', azimuthal = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("setstress", [meshname, magnitude, polar, azimuthal])
    	self.SendCommand("buffercommand", ["setstress", meshname, magnitude, polar, azimuthal])
    
    def shape_cone(self, meshname = '', len_x = '', len_y = '', len_z = '', cpos_x = '', cpos_y = '', cpos_z = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("shape_cone", [meshname, len_x, len_y, len_z, cpos_x, cpos_y, cpos_z])
    	self.SendCommand("buffercommand", ["shape_cone", meshname, len_x, len_y, len_z, cpos_x, cpos_y, cpos_z])
    
    def shape_disk(self, meshname = '', dia_x = '', dia_y = '', cpos_x = '', cpos_y = '', z_start = '', z_end = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("shape_disk", [meshname, dia_x, dia_y, cpos_x, cpos_y, z_start, z_end])
    	self.SendCommand("buffercommand", ["shape_disk", meshname, dia_x, dia_y, cpos_x, cpos_y, z_start, z_end])
    
    def shape_displacement(self, x = '', y = '', z = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("shape_displacement", [x, y, z])
    	self.SendCommand("buffercommand", ["shape_displacement", x, y, z])
    
    def shape_ellipsoid(self, meshname = '', dia_x = '', dia_y = '', dia_z = '', cpos_x = '', cpos_y = '', cpos_z = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("shape_ellipsoid", [meshname, dia_x, dia_y, dia_z, cpos_x, cpos_y, cpos_z])
    	self.SendCommand("buffercommand", ["shape_ellipsoid", meshname, dia_x, dia_y, dia_z, cpos_x, cpos_y, cpos_z])
    
    def shape_get(self, meshname = '', shape = Shape()):
        if isinstance(meshname, str): return self.SendCommand("shape_get", [meshname, shape.tostring()])
        else: return self.SendCommand("shape_get", [meshname.tostring()])
    
    def shape_method(self, method = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("shape_method", [method])
    	self.SendCommand("buffercommand", ["shape_method", method])
    
    def shape_pyramid(self, meshname = '', len_x = '', len_y = '', len_z = '', cpos_x = '', cpos_y = '', cpos_z = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("shape_pyramid", [meshname, len_x, len_y, len_z, cpos_x, cpos_y, cpos_z])
    	self.SendCommand("buffercommand", ["shape_pyramid", meshname, len_x, len_y, len_z, cpos_x, cpos_y, cpos_z])
    
    def shape_rect(self, meshname = '', len_x = '', len_y = '', cpos_x = '', cpos_y = '', z_start = '', z_end = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("shape_rect", [meshname, len_x, len_y, cpos_x, cpos_y, z_start, z_end])
    	self.SendCommand("buffercommand", ["shape_rect", meshname, len_x, len_y, cpos_x, cpos_y, z_start, z_end])
    
    def shape_repetitions(self, x = '', y = '', z = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("shape_repetitions", [x, y, z])
    	self.SendCommand("buffercommand", ["shape_repetitions", x, y, z])
    
    def shape_rotation(self, psi = '', theta = '', phi = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("shape_rotation", [psi, theta, phi])
    	self.SendCommand("buffercommand", ["shape_rotation", psi, theta, phi])
    
    def shape_set(self, meshname = '', shape = Shape()):
        if isinstance(meshname, str): return self.SendCommand("shape_set", [meshname, shape.tostring()])
        else: return self.SendCommand("shape_set", [meshname.tostring()])
    
    def shape_setangle(self, meshname = '', shape = Shape(), theta = '', polar = ''):
        if isinstance(meshname, str): return self.SendCommand("shape_setangle", [meshname, shape.tostring(), theta, polar])
        else: return self.SendCommand("shape_setangle", [meshname.tostring(), shape, theta])
    
    def shape_setparam(self, meshname = '', paramname = '', shape = Shape(), scaling_value = ''):
        if isinstance(paramname, str): return self.SendCommand("shape_setparam", [meshname, paramname, shape.tostring(), scaling_value])
        else: return self.SendCommand("shape_setparam", [meshname, paramname.tostring(), shape])
    
    def shape_tetrahedron(self, meshname = '', len_x = '', len_y = '', len_z = '', cpos_x = '', cpos_y = '', cpos_z = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("shape_tetrahedron", [meshname, len_x, len_y, len_z, cpos_x, cpos_y, cpos_z])
    	self.SendCommand("buffercommand", ["shape_tetrahedron", meshname, len_x, len_y, len_z, cpos_x, cpos_y, cpos_z])
    
    def shape_torus(self, meshname = '', len_x = '', len_y = '', len_z = '', cpos_x = '', cpos_y = '', cpos_z = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("shape_torus", [meshname, len_x, len_y, len_z, cpos_x, cpos_y, cpos_z])
    	self.SendCommand("buffercommand", ["shape_torus", meshname, len_x, len_y, len_z, cpos_x, cpos_y, cpos_z])
    
    def shape_triangle(self, meshname = '', len_x = '', len_y = '', cpos_x = '', cpos_y = '', z_start = '', z_end = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("shape_triangle", [meshname, len_x, len_y, cpos_x, cpos_y, z_start, z_end])
    	self.SendCommand("buffercommand", ["shape_triangle", meshname, len_x, len_y, cpos_x, cpos_y, z_start, z_end])
    
    def shiftcamorigin(self, dX = '', dY = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("shiftcamorigin", [dX, dY])
    	self.SendCommand("buffercommand", ["shiftcamorigin", dX, dY])
    
    def showa(self, meshname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("showa", [meshname])
    	self.SendCommand("buffercommand", ["showa", meshname])
    
    def showdata(self, meshname = '', dataname = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("showdata", [meshname, dataname, rectangle])
    	self.SendCommand("buffercommand", ["showdata", meshname, dataname, rectangle])
    
    def showk(self, meshname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("showk", [meshname])
    	self.SendCommand("buffercommand", ["showk", meshname])
    
    def showlengths(self, meshname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("showlengths", [meshname])
    	self.SendCommand("buffercommand", ["showlengths", meshname])
    
    def showmcells(self, meshname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("showmcells", [meshname])
    	self.SendCommand("buffercommand", ["showmcells", meshname])
    
    def showms(self, meshname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("showms", [meshname])
    	self.SendCommand("buffercommand", ["showms", meshname])
    
    def showtc(self, meshname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("showtc", [meshname])
    	self.SendCommand("buffercommand", ["showtc", meshname])
    
    def skyposdmul(self, meshname = '', multiplier = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("skyposdmul", [meshname, multiplier])
    	self.SendCommand("buffercommand", ["skyposdmul", meshname, multiplier])
    
    def skyrmion(self, meshname = '', core = '', chirality = '', diameter = '', position = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("skyrmion", [meshname, core, chirality, diameter, position])
    	self.SendCommand("buffercommand", ["skyrmion", meshname, core, chirality, diameter, position])
    
    def skyrmionbloch(self, meshname = '', core = '', chirality = '', diameter = '', position = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("skyrmionbloch", [meshname, core, chirality, diameter, position])
    	self.SendCommand("buffercommand", ["skyrmionbloch", meshname, core, chirality, diameter, position])
    
    def skyrmionpreparemovingmesh(self, meshname = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("skyrmionpreparemovingmesh", [meshname])
    	self.SendCommand("buffercommand", ["skyrmionpreparemovingmesh", meshname])
    
    def ssolverconfig(self, s_convergence_error = '', s_iters_timeout = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("ssolverconfig", [s_convergence_error, s_iters_timeout])
    	self.SendCommand("buffercommand", ["ssolverconfig", s_convergence_error, s_iters_timeout])
    
    def stages(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("stages")
    	self.SendCommand("buffercommand", ["stages"])
    
    def startupscriptserver(self, status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("startupscriptserver", [status])
    	self.SendCommand("buffercommand", ["startupscriptserver", status])
    
    def startupupdatecheck(self, status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("startupupdatecheck", [status])
    	self.SendCommand("buffercommand", ["startupupdatecheck", status])
    
    def statictransportsolver(self, status = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("statictransportsolver", [status])
    	self.SendCommand("buffercommand", ["statictransportsolver", status])
    
    def stochastic(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("stochastic")
    	self.SendCommand("buffercommand", ["stochastic"])
    
    def stop(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("stop")
    	self.SendCommand("buffercommand", ["stop"])
    
    def surfroughenjagged(self, meshname = '', depth = '', spacing = '', seed = '', sides = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("surfroughenjagged", [meshname, depth, spacing, seed, sides])
    	self.SendCommand("buffercommand", ["surfroughenjagged", meshname, depth, spacing, seed, sides])
    
    def tau(self, meshname = '', tau_11 = '', tau_22 = '', tau_12 = '', tau_21 = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("tau", [meshname, tau_11, tau_22, tau_12, tau_21])
    	self.SendCommand("buffercommand", ["tau", meshname, tau_11, tau_22, tau_12, tau_21])
    
    def tcellsize(self, meshname = '', value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("tcellsize", [meshname, value])
    	self.SendCommand("buffercommand", ["tcellsize", meshname, value])
    
    def temperature(self, meshname = '', value = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("temperature", [meshname, value])
    	self.SendCommand("buffercommand", ["temperature", meshname, value])
    
    def threads(self, number = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("threads", [number])
    	self.SendCommand("buffercommand", ["threads", number])
    
    def tmodel(self, meshname = '', num_temperatures = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("tmodel", [meshname, num_temperatures])
    	self.SendCommand("buffercommand", ["tmodel", meshname, num_temperatures])
    
    def tsolverconfig(self, convergence_error = '', iters_timeout = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("tsolverconfig", [convergence_error, iters_timeout])
    	self.SendCommand("buffercommand", ["tsolverconfig", convergence_error, iters_timeout])
    
    def updatemdb(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("updatemdb")
    	self.SendCommand("buffercommand", ["updatemdb"])
    
    def updatescreen(self, bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("updatescreen")
    	self.SendCommand("buffercommand", ["updatescreen"])
    
    def vecrep(self, meshname = '', vecreptype = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("vecrep", [meshname, vecreptype])
    	self.SendCommand("buffercommand", ["vecrep", meshname, vecreptype])
    
    def versionupdate(self, action = '', target_version = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("versionupdate", [action, target_version])
    	self.SendCommand("buffercommand", ["versionupdate", action, target_version])
    
    def vortex(self, meshname = '', longitudinal = '', rotation = '', core = '', rectangle = '', bufferCommand = False):
    	if not bufferCommand: return self.SendCommand("vortex", [meshname, longitudinal, rotation, core, rectangle])
    	self.SendCommand("buffercommand", ["vortex", meshname, longitudinal, rotation, core, rectangle])
    


    #################### Command Send / Data Receive #######################

    def SendCommand(self, command, values = None):

        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        
            sock.connect((self.scriptserverip, self.scriptserverport))
            sock.settimeout(self.timeout_ms / 1000)

            while True:
                message = command
    
                #command arguments if specified : space-separated
                if values is not None: 
                    for value in values: 
                        #first check if it's a string
                        if isinstance(value, str):
                            #only add it if not empty
                            if len(value): message += ' ' + value
                        #not a string : could be a list or a number
                        else:
                            #first try to convert entries in a list to add as space-seaprated parameters
                            try:
                                for entry in value: message += ' ' + str(entry)
                            #not a list, so add a number after converting to string
                            except:
                                message += ' ' + str(value)
                        
                try:
                    # Send data
                    if self.verbose == True: sock.sendall(bytes(self.scriptserverpwd + '>' + message, 'utf-8'))
                    else: sock.sendall(bytes(self.scriptserverpwd + '*' + message, 'utf-8'))
                    if self.script_verbose: print('TX : %s' % message)
                except:
                    print("SendCommand (send): timed out.")
                    return
    
                # Look for the response
                try:
                    data = str(sock.recv(self.maxLenMessage), 'utf-8')
                except:
                    print("SendCommand (receive): timed out.")
                    return
        
                #note, the returned data always starts with a tab
                fields = data.split('\t')
    
                if len(fields) >= 2 and fields[1] == 'stopped':
                        #if we received the 'stopped' message this means a simulation was running when we sent this command. 
                        #this caused the simulation to stop thus issuing the 'stopped' message since a client is connected.
                        #since this command is expecting another message to be returned, we need to receive it - issue recv call again
                        try:
                            data = str(sock.recv(self.maxLenMessage), 'utf-8')
                        except:
                            print("SendCommand (receive): timed out.")
                            return
    
                        #note, the returned data always starts with a tab
                        fields = data.split('\t')
    
                if self.script_verbose: print('RX : %s' % data)
    
                #the received message should always have at least 2 fields since the message always starts with a tab
                if len(fields) >= 2:      
    
                    #list of floats where there should be floats instead of strings (skip first entry always as this is empty)
                    return_data = [self.convert_returned_parameter(entry) for entry in fields[1:]]
                    #if the returned data has multiple parameters then return it as a list, else as a single element
                    if len(return_data) == 1: return return_data[0]
                    else: return return_data
    
###############################################################################################################################