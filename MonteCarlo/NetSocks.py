#NetSocks Module Updated on : 15/03/2021
#Boris version : >3.01 (in progress)
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

class NSClient:

    #################### DATA
    
    maxLenMessage = 4096
    timeout_ms = 30000000
    
    serverip = 'localhost'
    serverport = 1542
    serverpwd = ''
    
    #verbosity of Boris console
    verbose = False
    
    #verbosity of Python script (adjusted by configure)
    script_verbose = True

    #on Windows assume this is where Boris.exe is (should be if installed with installer)
    #on Linux don't attempt to define a default : user will have to provide path if they want automatic startup
    win_default_boris_path = 'C:/Program Files (x86)/Boris'

    #################### CTOR / DTOR

    def __init__(self, 
                 serverip = 'localhost', serverport = 1542, serverpwd = '', 
                 cudaDevice = -1, 
                 boris_path = '', boris_exe = '',
                 window = 'back',
                 verbose = False):

        self.serverip = serverip     
        self.serverport = serverport
        self.verbose = verbose
        self.serverpwd = serverpwd
        
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
            
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        
                try:
                    sock.connect((self.serverip, self.serverport))
                except:
                    if serverip == 'localhost':
                        print("No server found on port %d. Starting new instance." % serverport)
                        os.chdir(boris_path)
                        subprocess.Popen([boris_exe, str(serverport), str(cudaDevice), window, serverpwd])
                        os.chdir(os.path.dirname(sys.argv[0]) + "/")
                        
                        #now make sure server is running and ready to accept input
                        for tryidx in range(10):
                        
                            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock2:
                                
                                try:
                                    sock2.connect((self.serverip, self.serverport))
                                    print("New instance started.")
                                    break
                                except:
                                    time.sleep(0.1)
                        else:
                            print("Couldn't start new instance with required server port - start it manually.");
                                
                                
                    else:
                        print("No server found on port %d. Make sure remote host has a Boris instance running for given port and is accessible." % serverport)

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
        
        if self.serverip == 'localhost':
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
        
            sock.connect((self.serverip, self.serverport))
            
            # Look for the response
            try:
                if self.script_verbose: print('TX : run')
                sock.sendall(bytes(self.serverpwd + '*' + "run", 'utf-8'))
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
    
    def _2dmulticonvolution(self, status = ''):
    	return self.SendCommand("2dmulticonvolution", [status])
    
    def addafmesh(self, name = '', rectangle = ''):
    	return self.SendCommand("addafmesh", [name, rectangle])
    
    def addameshcubic(self, name = '', rectangle = ''):
    	return self.SendCommand("addameshcubic", [name, rectangle])
    
    def addconductor(self, name = '', rectangle = ''):
    	return self.SendCommand("addconductor", [name, rectangle])
    
    def adddata(self, dataname = '', meshname = '', rectangle = ''):
    	return self.SendCommand("adddata", [dataname, meshname, rectangle])
    
    def adddiamagnet(self, name = '', rectangle = ''):
    	return self.SendCommand("adddiamagnet", [name, rectangle])
    
    def adddipole(self, name = '', rectangle = ''):
    	return self.SendCommand("adddipole", [name, rectangle])
    
    def addelectrode(self, electrode_rect = ''):
    	return self.SendCommand("addelectrode", [electrode_rect])
    
    def addinsulator(self, name = '', rectangle = ''):
    	return self.SendCommand("addinsulator", [name, rectangle])
    
    def addmaterial(self, name = '', rectangle = ''):
    	return self.SendCommand("addmaterial", [name, rectangle])
    
    def addmdbentry(self, meshname = '', materialname = ''):
    	return self.SendCommand("addmdbentry", [meshname, materialname])
    
    def addmesh(self, name = '', rectangle = ''):
    	return self.SendCommand("addmesh", [name, rectangle])
    
    def addmodule(self, meshname = '', handle = ''):
    	return self.SendCommand("addmodule", [meshname, handle])
    
    def addpinneddata(self, dataname = '', meshname = '', rectangle = ''):
    	return self.SendCommand("addpinneddata", [dataname, meshname, rectangle])
    
    def addrect(self, rectangle = '', meshname = ''):
    	return self.SendCommand("addrect", [rectangle, meshname])
    
    def addstage(self, stagetype = '', meshname = ''):
    	return self.SendCommand("addstage", [stagetype, meshname])
    
    def ambient(self, ambient_temperature = '', meshname = ''):
    	return self.SendCommand("ambient", [ambient_temperature, meshname])
    
    def astepctrl(self, err_fail = '', err_high = '', err_low = '', dT_incr = '', dT_min = '', dT_max = ''):
    	return self.SendCommand("astepctrl", [err_fail, err_high, err_low, dT_incr, dT_min, dT_max])
    
    def atomicmoment(self, ub_multiple = '', meshname = ''):
    	return self.SendCommand("atomicmoment", [ub_multiple, meshname])
    
    def averagemeshrect(self, rectangle = ''):
    	return self.SendCommand("averagemeshrect", [rectangle])
    
    def benchtime(self):
    	return self.SendCommand("benchtime")
    
    def blochpreparemovingmesh(self, meshname = ''):
    	return self.SendCommand("blochpreparemovingmesh", [meshname])
    
    def cellsize(self, value = ''):
    	return self.SendCommand("cellsize", [value])
    
    def center(self):
    	return self.SendCommand("center")
    
    def chdir(self, directory = ''):
    	return self.SendCommand("chdir", [directory])
    
    def checkupdates(self):
    	return self.SendCommand("checkupdates")
    
    def clearelectrodes(self):
    	return self.SendCommand("clearelectrodes")
    
    def clearequationconstants(self):
    	return self.SendCommand("clearequationconstants")
    
    def clearmovingmesh(self):
    	return self.SendCommand("clearmovingmesh")
    
    def clearparamstemp(self, meshname = '', paramname = ''):
    	return self.SendCommand("clearparamstemp", [meshname, paramname])
    
    def clearparamsvar(self, meshname = ''):
    	return self.SendCommand("clearparamsvar", [meshname])
    
    def clearparamvar(self, meshname = '', paramname = ''):
    	return self.SendCommand("clearparamvar", [meshname, paramname])
    
    def clearroughness(self, meshname = ''):
    	return self.SendCommand("clearroughness", [meshname])
    
    def clearscreen(self):
    	return self.SendCommand("clearscreen")
    
    def computefields(self):
    	return self.SendCommand("computefields")
    
    def copymeshdata(self, meshname_from = '', meshname_to = ''):
    	return self.SendCommand("copymeshdata", [meshname_from, meshname_to])
    
    def copyparams(self, meshname_from = '', meshname_to = ''):
    	return self.SendCommand("copyparams", [meshname_from, meshname_to])
    
    def coupletodipoles(self, status = ''):
    	return self.SendCommand("coupletodipoles", [status])
    
    def crosstie(self, direction = '', radius = '', thickness = '', centre = '', meshname = ''):
    	return self.SendCommand("crosstie", [direction, radius, thickness, centre, meshname])
    
    def cuda(self, status = ''):
    	return self.SendCommand("cuda", [status])
    
    def curietemperature(self, curie_temperature = '', meshname = ''):
    	return self.SendCommand("curietemperature", [curie_temperature, meshname])
    
    def data(self):
    	return self.SendCommand("data")
    
    def dataprecision(self, precision = ''):
    	return self.SendCommand("dataprecision", [precision])
    
    def default(self):
    	return self.SendCommand("default")
    
    def deldata(self, index = ''):
    	return self.SendCommand("deldata", [index])
    
    def delelectrode(self, index = ''):
    	return self.SendCommand("delelectrode", [index])
    
    def delequationconstant(self, name = ''):
    	return self.SendCommand("delequationconstant", [name])
    
    def delmdbentry(self, materialname = ''):
    	return self.SendCommand("delmdbentry", [materialname])
    
    def delmesh(self, name = ''):
    	return self.SendCommand("delmesh", [name])
    
    def delmodule(self, meshname = '', handle = ''):
    	return self.SendCommand("delmodule", [meshname, handle])
    
    def delpinneddata(self, index = ''):
    	return self.SendCommand("delpinneddata", [index])
    
    def delrect(self, rectangle = '', meshname = ''):
    	return self.SendCommand("delrect", [rectangle, meshname])
    
    def delstage(self, index = ''):
    	return self.SendCommand("delstage", [index])
    
    def designateground(self, electrode_index = ''):
    	return self.SendCommand("designateground", [electrode_index])
    
    def disabletransportsolver(self, status = ''):
    	return self.SendCommand("disabletransportsolver", [status])
    
    def diskbufferlines(self, lines = ''):
    	return self.SendCommand("diskbufferlines", [lines])
    
    def display(self, name = '', meshname = ''):
    	return self.SendCommand("display", [name, meshname])
    
    def displaybackground(self, name = '', meshname = ''):
    	return self.SendCommand("displaybackground", [name, meshname])
    
    def displaydetail(self, size = ''):
    	return self.SendCommand("displaydetail", [size])
    
    def displaymodule(self, modulename = '', meshname = ''):
    	return self.SendCommand("displaymodule", [modulename, meshname])
    
    def displayrenderthresholds(self, thresh1 = '', thresh2 = '', thresh3 = ''):
    	return self.SendCommand("displayrenderthresholds", [thresh1, thresh2, thresh3])
    
    def displaythresholds(self, minimum = '', maximum = ''):
    	return self.SendCommand("displaythresholds", [minimum, maximum])
    
    def displaythresholdtrigger(self, trigtype = ''):
    	return self.SendCommand("displaythresholdtrigger", [trigtype])
    
    def displaytransparency(self, foreground = '', background = ''):
    	return self.SendCommand("displaytransparency", [foreground, background])
    
    def dmcellsize(self, value = ''):
    	return self.SendCommand("dmcellsize", [value])
    
    def dp_add(self, dp_source = '', value = '', dp_dest = ''):
    	return self.SendCommand("dp_add", [dp_source, value, dp_dest])
    
    def dp_adddp(self, dp_x1 = '', dp_x2 = '', dp_dest = ''):
    	return self.SendCommand("dp_adddp", [dp_x1, dp_x2, dp_dest])
    
    def dp_append(self, dp_original = '', dp_new = ''):
    	return self.SendCommand("dp_append", [dp_original, dp_new])
    
    def dp_calcsot(self, hm_mesh = '', fm_mesh = ''):
    	return self.SendCommand("dp_calcsot", [hm_mesh, fm_mesh])
    
    def dp_calctopochargedensity(self):
    	return self.SendCommand("dp_calctopochargedensity")
    
    def dp_cartesiantopolar(self, dp_in_x = '', dp_in_y = '', dp_out_r = '', dp_out_theta = ''):
    	return self.SendCommand("dp_cartesiantopolar", [dp_in_x, dp_in_y, dp_out_r, dp_out_theta])
    
    def dp_chunkedstd(self, dp_index = '', chunk = ''):
    	return self.SendCommand("dp_chunkedstd", [dp_index, chunk])
    
    def dp_clear(self, indexes = ''):
    	return self.SendCommand("dp_clear", [indexes])
    
    def dp_clearall(self):
    	return self.SendCommand("dp_clearall")
    
    def dp_coercivity(self, dp_index_x = '', dp_index_y = ''):
    	return self.SendCommand("dp_coercivity", [dp_index_x, dp_index_y])
    
    def dp_completehysteresis(self, dp_index_x = '', dp_index_y = ''):
    	return self.SendCommand("dp_completehysteresis", [dp_index_x, dp_index_y])
    
    def dp_countskyrmions(self, x = '', y = '', radius = ''):
    	return self.SendCommand("dp_countskyrmions", [x, y, radius])
    
    def dp_crossingsfrequency(self, dp_in_x = '', dp_in_y = '', dp_level = '', dp_freq_up = '', dp_freq_dn = '', steps = ''):
    	return self.SendCommand("dp_crossingsfrequency", [dp_in_x, dp_in_y, dp_level, dp_freq_up, dp_freq_dn, steps])
    
    def dp_crossingshistogram(self, dp_in_x = '', dp_in_y = '', dp_level = '', dp_counts = '', steps = ''):
    	return self.SendCommand("dp_crossingshistogram", [dp_in_x, dp_in_y, dp_level, dp_counts, steps])
    
    def dp_div(self, dp_source = '', value = '', dp_dest = ''):
    	return self.SendCommand("dp_div", [dp_source, value, dp_dest])
    
    def dp_divdp(self, dp_x1 = '', dp_x2 = '', dp_dest = ''):
    	return self.SendCommand("dp_divdp", [dp_x1, dp_x2, dp_dest])
    
    def dp_dotprod(self, dp_vector = '', ux = '', uy = '', uz = '', dp_out = ''):
    	return self.SendCommand("dp_dotprod", [dp_vector, ux, uy, uz, dp_out])
    
    def dp_dotproddp(self, dp_x1 = '', dp_x2 = ''):
    	return self.SendCommand("dp_dotproddp", [dp_x1, dp_x2])
    
    def dp_dumptdep(self, meshname = '', paramname = '', max_temperature = '', dp_index = ''):
    	return self.SendCommand("dp_dumptdep", [meshname, paramname, max_temperature, dp_index])
    
    def dp_erase(self, dp_index = '', start_index = '', length = ''):
    	return self.SendCommand("dp_erase", [dp_index, start_index, length])
    
    def dp_extract(self, dp_in = '', dp_out = '', start_index = '', length = ''):
    	return self.SendCommand("dp_extract", [dp_in, dp_out, start_index, length])
    
    def dp_fitadiabatic(self, abs_err = '', Rsq = '', T_ratio = '', stencil = ''):
    	return self.SendCommand("dp_fitadiabatic", [abs_err, Rsq, T_ratio, stencil])
    
    def dp_fitdw(self, dp_x = '', dp_y = ''):
    	return self.SendCommand("dp_fitdw", [dp_x, dp_y])
    
    def dp_fitlorentz(self, dp_x = '', dp_y = ''):
    	return self.SendCommand("dp_fitlorentz", [dp_x, dp_y])
    
    def dp_fitlorentz2(self, dp_x = '', dp_y = ''):
    	return self.SendCommand("dp_fitlorentz2", [dp_x, dp_y])
    
    def dp_fitnonadiabatic(self, abs_err = '', Rsq = '', T_ratio = '', stencil = ''):
    	return self.SendCommand("dp_fitnonadiabatic", [abs_err, Rsq, T_ratio, stencil])
    
    def dp_fitskyrmion(self, dp_x = '', dp_y = ''):
    	return self.SendCommand("dp_fitskyrmion", [dp_x, dp_y])
    
    def dp_fitsot(self, hm_mesh = '', rectangle = ''):
    	return self.SendCommand("dp_fitsot", [hm_mesh, rectangle])
    
    def dp_fitsotstt(self, hm_mesh = '', rectangle = ''):
    	return self.SendCommand("dp_fitsotstt", [hm_mesh, rectangle])
    
    def dp_fitstt(self, rectangle = ''):
    	return self.SendCommand("dp_fitstt", [rectangle])
    
    def dp_get(self, dp_arr = '', index = ''):
    	return self.SendCommand("dp_get", [dp_arr, index])
    
    def dp_getampli(self, dp_source = '', pointsPeriod = ''):
    	return self.SendCommand("dp_getampli", [dp_source, pointsPeriod])
    
    def dp_getexactprofile(self, start = '', end = '', step = '', dp_index = '', stencil = ''):
    	return self.SendCommand("dp_getexactprofile", [start, end, step, dp_index, stencil])
    
    def dp_getprofile(self, start = '', end = '', dp_index = ''):
    	return self.SendCommand("dp_getprofile", [start, end, dp_index])
    
    def dp_histogram(self, dp_x = '', dp_y = '', bin = '', min = '', max = ''):
    	return self.SendCommand("dp_histogram", [dp_x, dp_y, bin, min, max])
    
    def dp_histogram2(self, dp_x = '', dp_y = '', bin = '', min = '', max = '', M2 = '', deltaM2 = ''):
    	return self.SendCommand("dp_histogram2", [dp_x, dp_y, bin, min, max, M2, deltaM2])
    
    def dp_linreg(self, dp_index_x = '', dp_index_y = '', dp_index_z = '', dp_index_out = ''):
    	return self.SendCommand("dp_linreg", [dp_index_x, dp_index_y, dp_index_z, dp_index_out])
    
    def dp_load(self, filename = '', file_indexes = '', dp_indexes = ''):
    	return self.SendCommand("dp_load", [filename, file_indexes, dp_indexes])
    
    def dp_mean(self, dp_index = '', exclusion_ratio = ''):
    	return self.SendCommand("dp_mean", [dp_index, exclusion_ratio])
    
    def dp_minmax(self, dp_index = ''):
    	return self.SendCommand("dp_minmax", [dp_index])
    
    def dp_monotonic(self, dp_in_x = '', dp_in_y = '', dp_out_x = '', dp_out_y = ''):
    	return self.SendCommand("dp_monotonic", [dp_in_x, dp_in_y, dp_out_x, dp_out_y])
    
    def dp_mul(self, dp_source = '', value = '', dp_dest = ''):
    	return self.SendCommand("dp_mul", [dp_source, value, dp_dest])
    
    def dp_muldp(self, dp_x1 = '', dp_x2 = '', dp_dest = ''):
    	return self.SendCommand("dp_muldp", [dp_x1, dp_x2, dp_dest])
    
    def dp_newfile(self, filename = ''):
    	return self.SendCommand("dp_newfile", [filename])
    
    def dp_peaksfrequency(self, dp_in_x = '', dp_in_y = '', dp_level = '', dp_freq = '', steps = ''):
    	return self.SendCommand("dp_peaksfrequency", [dp_in_x, dp_in_y, dp_level, dp_freq, steps])
    
    def dp_pow(self, dp_source = '', exponent = '', dp_dest = ''):
    	return self.SendCommand("dp_pow", [dp_source, exponent, dp_dest])
    
    def dp_rarefy(self, dp_in = '', dp_out = '', skip = ''):
    	return self.SendCommand("dp_rarefy", [dp_in, dp_out, skip])
    
    def dp_remanence(self, dp_index_x = '', dp_index_y = ''):
    	return self.SendCommand("dp_remanence", [dp_index_x, dp_index_y])
    
    def dp_removeoffset(self, dp_index = '', dp_index_out = ''):
    	return self.SendCommand("dp_removeoffset", [dp_index, dp_index_out])
    
    def dp_replacerepeats(self, dp_index = '', dp_index_out = ''):
    	return self.SendCommand("dp_replacerepeats", [dp_index, dp_index_out])
    
    def dp_save(self, filename = '', dp_indexes = ''):
    	return self.SendCommand("dp_save", [filename, dp_indexes])
    
    def dp_saveappend(self, filename = '', dp_indexes = ''):
    	return self.SendCommand("dp_saveappend", [filename, dp_indexes])
    
    def dp_saveappendasrow(self, filename = '', dp_index = ''):
    	return self.SendCommand("dp_saveappendasrow", [filename, dp_index])
    
    def dp_saveasrow(self, filename = '', dp_index = ''):
    	return self.SendCommand("dp_saveasrow", [filename, dp_index])
    
    def dp_sequence(self, dp_index = '', start_value = '', increment = '', points = ''):
    	return self.SendCommand("dp_sequence", [dp_index, start_value, increment, points])
    
    def dp_set(self, dp_arr = '', index = '', value = ''):
    	return self.SendCommand("dp_set", [dp_arr, index, value])
    
    def dp_showsizes(self, dp_arr = ''):
    	return self.SendCommand("dp_showsizes", [dp_arr])
    
    def dp_smooth(self, dp_in = '', dp_out = '', window_size = ''):
    	return self.SendCommand("dp_smooth", [dp_in, dp_out, window_size])
    
    def dp_sub(self, dp_source = '', value = '', dp_dest = ''):
    	return self.SendCommand("dp_sub", [dp_source, value, dp_dest])
    
    def dp_subdp(self, dp_x1 = '', dp_x2 = '', dp_dest = ''):
    	return self.SendCommand("dp_subdp", [dp_x1, dp_x2, dp_dest])
    
    def dp_topocharge(self, x = '', y = '', radius = ''):
    	return self.SendCommand("dp_topocharge", [x, y, radius])
    
    def dwall(self, longitudinal = '', transverse = '', width = '', position = '', meshname = ''):
    	return self.SendCommand("dwall", [longitudinal, transverse, width, position, meshname])
    
    def dwposcomponent(self, value = ''):
    	return self.SendCommand("dwposcomponent", [value])
    
    def ecellsize(self, value = ''):
    	return self.SendCommand("ecellsize", [value])
    
    def editdata(self, index = '', dataname = '', meshname = '', rectangle = ''):
    	return self.SendCommand("editdata", [index, dataname, meshname, rectangle])
    
    def editdatasave(self, index = '', savetype = '', savevalue = ''):
    	return self.SendCommand("editdatasave", [index, savetype, savevalue])
    
    def editstage(self, index = '', stagetype = '', meshname = ''):
    	return self.SendCommand("editstage", [index, stagetype, meshname])
    
    def editstagestop(self, index = '', stoptype = '', stopvalue = ''):
    	return self.SendCommand("editstagestop", [index, stoptype, stopvalue])
    
    def editstagevalue(self, index = '', value = ''):
    	return self.SendCommand("editstagevalue", [index, value])
    
    def electrodes(self):
    	return self.SendCommand("electrodes")
    
    def equationconstants(self, name = '', value = ''):
    	return self.SendCommand("equationconstants", [name, value])
    
    def errorlog(self, status = ''):
    	return self.SendCommand("errorlog", [status])
    
    def escellsize(self, value = ''):
    	return self.SendCommand("escellsize", [value])
    
    def evalspeedup(self, level = ''):
    	return self.SendCommand("evalspeedup", [level])
    
    def exchangecoupledmeshes(self, status = '', meshname = ''):
    	return self.SendCommand("exchangecoupledmeshes", [status, meshname])
    
    def excludemulticonvdemag(self, status = '', meshname = ''):
    	return self.SendCommand("excludemulticonvdemag", [status, meshname])
    
    def flower(self, direction = '', radius = '', thickness = '', centre = '', meshname = ''):
    	return self.SendCommand("flower", [direction, radius, thickness, centre, meshname])
    
    def flusherrorlog(self):
    	return self.SendCommand("flusherrorlog")
    
    def fmscellsize(self, value = ''):
    	return self.SendCommand("fmscellsize", [value])
    
    def generate2dgrains(self, spacing = '', seed = ''):
    	return self.SendCommand("generate2dgrains", [spacing, seed])
    
    def generate3dgrains(self, spacing = '', seed = ''):
    	return self.SendCommand("generate3dgrains", [spacing, seed])
    
    def getvalue(self, abspos = ''):
    	return self.SendCommand("getvalue", [abspos])
    
    def imagecropping(self, left = '', bottom = '', right = '', top = ''):
    	return self.SendCommand("imagecropping", [left, bottom, right, top])
    
    def individualmaskshape(self, status = ''):
    	return self.SendCommand("individualmaskshape", [status])
    
    def insulatingside(self, side_literal = '', status = '', meshname = ''):
    	return self.SendCommand("insulatingside", [side_literal, status, meshname])
    
    def invertmag(self, components = '', meshname = ''):
    	return self.SendCommand("invertmag", [components, meshname])
    
    def isrunning(self):
    	return self.SendCommand("isrunning")
    
    def iterupdate(self, iterations = ''):
    	return self.SendCommand("iterupdate", [iterations])
    
    def linkdtspeedup(self, flag = ''):
    	return self.SendCommand("linkdtspeedup", [flag])
    
    def linkdtstochastic(self, flag = ''):
    	return self.SendCommand("linkdtstochastic", [flag])
    
    def linkstochastic(self, flag = '', meshname = ''):
    	return self.SendCommand("linkstochastic", [flag, meshname])
    
    def loadmaskfile(self, z_depth = '', filename = ''):
    	return self.SendCommand("loadmaskfile", [z_depth, filename])
    
    def loadovf2curr(self, filename = ''):
    	return self.SendCommand("loadovf2curr", [filename])
    
    def loadovf2disp(self, filename = ''):
    	return self.SendCommand("loadovf2disp", [filename])
    
    def loadovf2mag(self, renormalize_value = '', filename = ''):
    	return self.SendCommand("loadovf2mag", [renormalize_value, filename])
    
    def loadovf2mesh(self, renormalize_value = '', filename = ''):
    	return self.SendCommand("loadovf2mesh", [renormalize_value, filename])
    
    def loadovf2strain(self, filename_diag = '', filename_odiag = ''):
    	return self.SendCommand("loadovf2strain", [filename_diag, filename_odiag])
    
    def loadovf2temp(self, filename = ''):
    	return self.SendCommand("loadovf2temp", [filename])
    
    def loadsim(self, filename = ''):
    	return self.SendCommand("loadsim", [filename])
    
    def makevideo(self, filebase = '', fps = '', quality = ''):
    	return self.SendCommand("makevideo", [filebase, fps, quality])
    
    def manual(self):
    	return self.SendCommand("manual")
    
    def matcurietemperature(self, curie_temperature = '', meshname = ''):
    	return self.SendCommand("matcurietemperature", [curie_temperature, meshname])
    
    def materialsdatabase(self, mdbname = ''):
    	return self.SendCommand("materialsdatabase", [mdbname])
    
    def mccomputefields(self, status = ''):
    	return self.SendCommand("mccomputefields", [status])
    
    def mcconeangle(self, min_angle = '', max_angle = ''):
    	return self.SendCommand("mcconeangle", [min_angle, max_angle])
    
    def mcconstrain(self, value = '', meshname = ''):
    	return self.SendCommand("mcconstrain", [value, meshname])
    
    def mcellsize(self, value = ''):
    	return self.SendCommand("mcellsize", [value])
    
    def mcserial(self, value = '', meshname = ''):
    	return self.SendCommand("mcserial", [value, meshname])
    
    def memory(self):
    	return self.SendCommand("memory")
    
    def mesh(self):
    	return self.SendCommand("mesh")
    
    def meshfocus(self, meshname = ''):
    	return self.SendCommand("meshfocus", [meshname])
    
    def meshfocus2(self, meshname = ''):
    	return self.SendCommand("meshfocus2", [meshname])
    
    def meshrect(self, rectangle = ''):
    	return self.SendCommand("meshrect", [rectangle])
    
    def mirrormag(self, axis = '', meshname = ''):
    	return self.SendCommand("mirrormag", [axis, meshname])
    
    def modules(self):
    	return self.SendCommand("modules")
    
    def movingmesh(self, status_or_meshname = ''):
    	return self.SendCommand("movingmesh", [status_or_meshname])
    
    def movingmeshasym(self, status = ''):
    	return self.SendCommand("movingmeshasym", [status])
    
    def movingmeshthresh(self, value = ''):
    	return self.SendCommand("movingmeshthresh", [value])
    
    def multiconvolution(self, status = ''):
    	return self.SendCommand("multiconvolution", [status])
    
    def ncommon(self, sizes = ''):
    	return self.SendCommand("ncommon", [sizes])
    
    def ncommonstatus(self, status = ''):
    	return self.SendCommand("ncommonstatus", [status])
    
    def neelpreparemovingmesh(self, meshname = ''):
    	return self.SendCommand("neelpreparemovingmesh", [meshname])
    
    def newinstance(self, port = '', cudaDevice = '', password = ''):
    	return self.SendCommand("newinstance", [port, cudaDevice, password])
    
    def ode(self):
    	return self.SendCommand("ode")
    
    def onion(self, direction = '', radius1 = '', radius2 = '', thickness = '', centre = '', meshname = ''):
    	return self.SendCommand("onion", [direction, radius1, radius2, thickness, centre, meshname])
    
    def params(self, meshname = ''):
    	return self.SendCommand("params", [meshname])
    
    def paramstemp(self, meshname = ''):
    	return self.SendCommand("paramstemp", [meshname])
    
    def paramsvar(self, meshname = ''):
    	return self.SendCommand("paramsvar", [meshname])
    
    def pbc(self, meshname = '', flag = '', images = ''):
    	return self.SendCommand("pbc", [meshname, flag, images])
    
    def preparemovingmesh(self, meshname = ''):
    	return self.SendCommand("preparemovingmesh", [meshname])
    
    def random(self, meshname = '', seed = ''):
    	return self.SendCommand("random", [meshname, seed])
    
    def randomxy(self, meshname = '', seed = ''):
    	return self.SendCommand("randomxy", [meshname, seed])
    
    def refineroughness(self, value = '', meshname = ''):
    	return self.SendCommand("refineroughness", [value, meshname])
    
    def refreshmdb(self):
    	return self.SendCommand("refreshmdb")
    
    def refreshscreen(self):
    	return self.SendCommand("refreshscreen")
    
    def renamemesh(self, old_name = '', new_name = ''):
    	return self.SendCommand("renamemesh", [old_name, new_name])
    
    def requestmdbsync(self, materialname = '', email = ''):
    	return self.SendCommand("requestmdbsync", [materialname, email])
    
    def reset(self):
    	return self.SendCommand("reset")
    
    def resetmesh(self, meshname = ''):
    	return self.SendCommand("resetmesh", [meshname])
    
    def robinalpha(self, robin_alpha = '', meshname = ''):
    	return self.SendCommand("robinalpha", [robin_alpha, meshname])
    
    def roughenmesh(self, depth = '', side = '', seed = ''):
    	return self.SendCommand("roughenmesh", [depth, side, seed])
    
    def savecomment(self, filename = '', comment = ''):
    	return self.SendCommand("savecomment", [filename, comment])
    
    def savedatafile(self, filename = ''):
    	return self.SendCommand("savedatafile", [filename])
    
    def savedataflag(self, status = ''):
    	return self.SendCommand("savedataflag", [status])
    
    def saveimagefile(self, filename = ''):
    	return self.SendCommand("saveimagefile", [filename])
    
    def saveimageflag(self, status = ''):
    	return self.SendCommand("saveimageflag", [status])
    
    def savemeshimage(self, filename = ''):
    	return self.SendCommand("savemeshimage", [filename])
    
    def saveovf2(self, data_type = '', filename = ''):
    	return self.SendCommand("saveovf2", [data_type, filename])
    
    def saveovf2mag(self, n = '', data_type = '', filename = ''):
    	return self.SendCommand("saveovf2mag", [n, data_type, filename])
    
    def saveovf2param(self, data_type = '', meshname = '', paramname = '', filename = ''):
    	return self.SendCommand("saveovf2param", [data_type, meshname, paramname, filename])
    
    def savesim(self, filename = ''):
    	return self.SendCommand("savesim", [filename])
    
    def scalemeshrects(self, status = ''):
    	return self.SendCommand("scalemeshrects", [status])
    
    def scellsize(self, value = ''):
    	return self.SendCommand("scellsize", [value])
    
    def scriptserver(self, status = ''):
    	return self.SendCommand("scriptserver", [status])
    
    def selectcudadevice(self, number = ''):
    	return self.SendCommand("selectcudadevice", [number])
    
    def serverpassword(self, password = ''):
    	return self.SendCommand("serverpassword", [password])
    
    def serverport(self, port = ''):
    	return self.SendCommand("serverport", [port])
    
    def serversleepms(self, time_ms = ''):
    	return self.SendCommand("serversleepms", [time_ms])
    
    def setafmesh(self, name = '', rectangle = ''):
    	return self.SendCommand("setafmesh", [name, rectangle])
    
    def setameshcubic(self, name = '', rectangle = ''):
    	return self.SendCommand("setameshcubic", [name, rectangle])
    
    def setangle(self, polar = '', azimuthal = '', meshname = ''):
    	return self.SendCommand("setangle", [polar, azimuthal, meshname])
    
    def setatomode(self, equation = '', evaluation = ''):
    	return self.SendCommand("setatomode", [equation, evaluation])
    
    def setcurrent(self, current = ''):
    	return self.SendCommand("setcurrent", [current])
    
    def setcurrentdensity(self, Jx = '', Jy = '', Jz = '', meshname = ''):
    	return self.SendCommand("setcurrentdensity", [Jx, Jy, Jz, meshname])
    
    def setdata(self, dataname = '', meshname = '', rectangle = ''):
    	return self.SendCommand("setdata", [dataname, meshname, rectangle])
    
    def setdefaultelectrodes(self, sides = ''):
    	return self.SendCommand("setdefaultelectrodes", [sides])
    
    def setdisplayedparamsvar(self, meshname = '', paramname = ''):
    	return self.SendCommand("setdisplayedparamsvar", [meshname, paramname])
    
    def setdt(self, value = ''):
    	return self.SendCommand("setdt", [value])
    
    def setdtspeedup(self, value = ''):
    	return self.SendCommand("setdtspeedup", [value])
    
    def setdtstoch(self, value = ''):
    	return self.SendCommand("setdtstoch", [value])
    
    def setelectrodepotential(self, electrode_index = '', potential = ''):
    	return self.SendCommand("setelectrodepotential", [electrode_index, potential])
    
    def setelectroderect(self, electrode_index = '', electrode_rect = ''):
    	return self.SendCommand("setelectroderect", [electrode_index, electrode_rect])
    
    def setfield(self, magnitude = '', polar = '', azimuthal = '', meshname = ''):
    	return self.SendCommand("setfield", [magnitude, polar, azimuthal, meshname])
    
    def setheatdt(self, value = ''):
    	return self.SendCommand("setheatdt", [value])
    
    def setktens(self, term1 = '', term2 = ''):
    	return self.SendCommand("setktens", [term1, term2])
    
    def setmaterial(self, name = '', rectangle = ''):
    	return self.SendCommand("setmaterial", [name, rectangle])
    
    def setmesh(self, name = '', rectangle = ''):
    	return self.SendCommand("setmesh", [name, rectangle])
    
    def setobjectangle(self, polar = '', azimuthal = '', position = '', meshname = ''):
    	return self.SendCommand("setobjectangle", [polar, azimuthal, position, meshname])
    
    def setode(self, equation = '', evaluation = ''):
    	return self.SendCommand("setode", [equation, evaluation])
    
    def setodeeval(self, evaluation = ''):
    	return self.SendCommand("setodeeval", [evaluation])
    
    def setparam(self, meshname = '', paramname = '', value = ''):
    	return self.SendCommand("setparam", [meshname, paramname, value])
    
    def setparamtemparray(self, meshname = '', paramname = '', filename = ''):
    	return self.SendCommand("setparamtemparray", [meshname, paramname, filename])
    
    def setparamtempequation(self, meshname = '', paramname = '', text_equation = ''):
    	return self.SendCommand("setparamtempequation", [meshname, paramname, text_equation])
    
    def setparamvar(self, meshname = '', paramname = '', generatorname = '', arguments = ''):
    	return self.SendCommand("setparamvar", [meshname, paramname, generatorname, arguments])
    
    def setpotential(self, potential = ''):
    	return self.SendCommand("setpotential", [potential])
    
    def setrect(self, polar = '', azimuthal = '', rectangle = '', meshname = ''):
    	return self.SendCommand("setrect", [polar, azimuthal, rectangle, meshname])
    
    def setsordamping(self, damping_v = '', damping_s = ''):
    	return self.SendCommand("setsordamping", [damping_v, damping_s])
    
    def setstage(self, stagetype = '', meshname = ''):
    	return self.SendCommand("setstage", [stagetype, meshname])
    
    def setstress(self, magnitude = '', polar = '', azimuthal = '', meshname = ''):
    	return self.SendCommand("setstress", [magnitude, polar, azimuthal, meshname])
    
    def shape_cone(self, len_x = '', len_y = '', len_z = '', cpos_x = '', cpos_y = '', cpos_z = ''):
    	return self.SendCommand("shape_cone", [len_x, len_y, len_z, cpos_x, cpos_y, cpos_z])
    
    def shape_disk(self, dia_x = '', dia_y = '', cpos_x = '', cpos_y = '', z_start = '', z_end = ''):
    	return self.SendCommand("shape_disk", [dia_x, dia_y, cpos_x, cpos_y, z_start, z_end])
    
    def shape_displacement(self, x = '', y = '', z = ''):
    	return self.SendCommand("shape_displacement", [x, y, z])
    
    def shape_ellipsoid(self, dia_x = '', dia_y = '', dia_z = '', cpos_x = '', cpos_y = '', cpos_z = ''):
    	return self.SendCommand("shape_ellipsoid", [dia_x, dia_y, dia_z, cpos_x, cpos_y, cpos_z])
    
    def shape_method(self, method = ''):
    	return self.SendCommand("shape_method", [method])
    
    def shape_pyramid(self, len_x = '', len_y = '', len_z = '', cpos_x = '', cpos_y = '', cpos_z = ''):
    	return self.SendCommand("shape_pyramid", [len_x, len_y, len_z, cpos_x, cpos_y, cpos_z])
    
    def shape_rect(self, len_x = '', len_y = '', cpos_x = '', cpos_y = '', z_start = '', z_end = ''):
    	return self.SendCommand("shape_rect", [len_x, len_y, cpos_x, cpos_y, z_start, z_end])
    
    def shape_repetitions(self, x = '', y = '', z = ''):
    	return self.SendCommand("shape_repetitions", [x, y, z])
    
    def shape_rotation(self, psi = '', theta = '', phi = ''):
    	return self.SendCommand("shape_rotation", [psi, theta, phi])
    
    def shape_set(self, name = '', dim_x = '', dim_y = '', dim_z = '', cpos_x = '', cpos_y = '', cpos_z = '', rot_psi = '', rot_theta = '', rot_phi = '', repeat_x = '', repeat_y = '', repeat_z = '', disp_x = '', disp_y = '', disp_z = '', method = ''):
    	return self.SendCommand("shape_set", [name, dim_x, dim_y, dim_z, cpos_x, cpos_y, cpos_z, rot_psi, rot_theta, rot_phi, repeat_x, repeat_y, repeat_z, disp_x, disp_y, disp_z, method])
    
    def shape_setangle(self, name = '', dim_x = '', dim_y = '', dim_z = '', cpos_x = '', cpos_y = '', cpos_z = '', rot_psi = '', rot_theta = '', rot_phi = '', repeat_x = '', repeat_y = '', repeat_z = '', disp_x = '', disp_y = '', disp_z = '', method = '', theta = '', polar = ''):
    	return self.SendCommand("shape_setangle", [name, dim_x, dim_y, dim_z, cpos_x, cpos_y, cpos_z, rot_psi, rot_theta, rot_phi, repeat_x, repeat_y, repeat_z, disp_x, disp_y, disp_z, method, theta, polar])
    
    def shape_setparam(self, paramname = '', name = '', dim_x = '', dim_y = '', dim_z = '', cpos_x = '', cpos_y = '', cpos_z = '', rot_psi = '', rot_theta = '', rot_phi = '', repeat_x = '', repeat_y = '', repeat_z = '', disp_x = '', disp_y = '', disp_z = '', method = '', scaling_value = ''):
    	return self.SendCommand("shape_setparam", [paramname, name, dim_x, dim_y, dim_z, cpos_x, cpos_y, cpos_z, rot_psi, rot_theta, rot_phi, repeat_x, repeat_y, repeat_z, disp_x, disp_y, disp_z, method, scaling_value])
    
    def shape_tetrahedron(self, len_x = '', len_y = '', len_z = '', cpos_x = '', cpos_y = '', cpos_z = ''):
    	return self.SendCommand("shape_tetrahedron", [len_x, len_y, len_z, cpos_x, cpos_y, cpos_z])
    
    def shape_torus(self, len_x = '', len_y = '', len_z = '', cpos_x = '', cpos_y = '', cpos_z = ''):
    	return self.SendCommand("shape_torus", [len_x, len_y, len_z, cpos_x, cpos_y, cpos_z])
    
    def shape_triangle(self, len_x = '', len_y = '', cpos_x = '', cpos_y = '', z_start = '', z_end = ''):
    	return self.SendCommand("shape_triangle", [len_x, len_y, cpos_x, cpos_y, z_start, z_end])
    
    def showa(self):
    	return self.SendCommand("showa")
    
    def showdata(self, dataname = '', meshname = '', rectangle = ''):
    	return self.SendCommand("showdata", [dataname, meshname, rectangle])
    
    def showk(self):
    	return self.SendCommand("showk")
    
    def showlengths(self):
    	return self.SendCommand("showlengths")
    
    def showmcells(self):
    	return self.SendCommand("showmcells")
    
    def showms(self):
    	return self.SendCommand("showms")
    
    def showtc(self):
    	return self.SendCommand("showtc")
    
    def skyposdmul(self, multiplier = '', meshname = ''):
    	return self.SendCommand("skyposdmul", [multiplier, meshname])
    
    def skyrmion(self, core = '', chirality = '', diameter = '', position = '', meshname = ''):
    	return self.SendCommand("skyrmion", [core, chirality, diameter, position, meshname])
    
    def skyrmionbloch(self, core = '', chirality = '', diameter = '', position = '', meshname = ''):
    	return self.SendCommand("skyrmionbloch", [core, chirality, diameter, position, meshname])
    
    def skyrmionpreparemovingmesh(self, meshname = ''):
    	return self.SendCommand("skyrmionpreparemovingmesh", [meshname])
    
    def ssolverconfig(self, s_convergence_error = '', s_iters_timeout = ''):
    	return self.SendCommand("ssolverconfig", [s_convergence_error, s_iters_timeout])
    
    def stages(self):
    	return self.SendCommand("stages")
    
    def startupscriptserver(self, status = ''):
    	return self.SendCommand("startupscriptserver", [status])
    
    def startupupdatecheck(self, status = ''):
    	return self.SendCommand("startupupdatecheck", [status])
    
    def statictransportsolver(self, status = ''):
    	return self.SendCommand("statictransportsolver", [status])
    
    def stochastic(self):
    	return self.SendCommand("stochastic")
    
    def stop(self):
    	return self.SendCommand("stop")
    
    def surfroughenjagged(self, depth = '', spacing = '', seed = '', sides = ''):
    	return self.SendCommand("surfroughenjagged", [depth, spacing, seed, sides])
    
    def tau(self, tau_11 = '', tau_22 = '', tau_12 = '', tau_21 = '', meshname = ''):
    	return self.SendCommand("tau", [tau_11, tau_22, tau_12, tau_21, meshname])
    
    def tcellsize(self, value = ''):
    	return self.SendCommand("tcellsize", [value])
    
    def temperature(self, value = '', meshname = ''):
    	return self.SendCommand("temperature", [value, meshname])
    
    def threads(self, number = ''):
    	return self.SendCommand("threads", [number])
    
    def tmodel(self, num_temperatures = '', meshname = ''):
    	return self.SendCommand("tmodel", [num_temperatures, meshname])
    
    def tsolverconfig(self, convergence_error = '', iters_timeout = ''):
    	return self.SendCommand("tsolverconfig", [convergence_error, iters_timeout])
    
    def updatemdb(self):
    	return self.SendCommand("updatemdb")
    
    def updatescreen(self):
    	return self.SendCommand("updatescreen")
    
    def vecrep(self, meshname = '', vecreptype = ''):
    	return self.SendCommand("vecrep", [meshname, vecreptype])
    
    def vortex(self, longitudinal = '', rotation = '', core = '', rectangle = '', meshname = ''):
    	return self.SendCommand("vortex", [longitudinal, rotation, core, rectangle, meshname])
    
    
    
    #################### Command Send / Data Receive #######################

    def SendCommand(self, command, values = None):

        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        
            sock.connect((self.serverip, self.serverport))
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
                    if self.verbose == True: sock.sendall(bytes(self.serverpwd + '>' + message, 'utf-8'))
                    else: sock.sendall(bytes(self.serverpwd + '*' + message, 'utf-8'))
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