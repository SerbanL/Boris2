#NetSocks Module Updated on : 12/07/2020
#Boris version : 2.8

import socket
import time
import matplotlib.pyplot as plt
import numpy as np

###############################################################################################################################

class NSClient:

    #################### DATA #######################
    
    # Create a TCP/IP socket
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

    maxLenMessage = 4096
    timeout_ms = 300000
    
    serverip = ''
    defaultPort = 1542
    
    connected = False
    
    pollTimer_ms = 0

    #################### CTOR / DTOR #######################

    def __init__(self, serverip = 'localhost', verbose = False):

        self.serverip = serverip        
        self.verbose = verbose

        print('connecting to %s:%s' % (serverip, self.defaultPort))

        try:
            self.sock.connect((serverip, self.defaultPort))
            self.sock.settimeout(self.timeout_ms / 1000)
            self.connected = True
        except:
            self.connected = False
            
    def __del__(self):
        #make sure any previous instance has the socket closed before starting a new one
        self.sock.close()

    #################### AUXILIARY #######################

    #use this on return parameters : returned parameters can either be numbers or a word (text without spaces)
    #the words by themselves cannot be converted to numbers
    #thus try to convert to a number first, if not must be a word
    def Convert_Returned_Parameter(self, string):
        try:    return float(string)
        except: return string

    #wait for data to be returned by server using a blocking socket (typically used with the run command to wait for simulation finished signal)
    def WaitForResponse(self):
        
        try: 
            self.sock.setblocking(True)
            data = str(self.sock.recv(self.maxLenMessage), 'utf-8')
            print('RX : %s' % data)
            self.sock.settimeout(self.timeout_ms / 1000)
        except:
            print("WaitForResponse: failed.")
            self.sock.close()
            
    #SendCommand also serves as auxiliary method
    
    #check if we can convert the row text to list of numbers; expecting tab-spaced data
    def can_convert(self, row_text):
        
        for entry in row_text.rstrip().split('\t'):
            try:    float(entry)
            except: return False

        return True
    
    #################### PLOTTING HELPERS #######################
    
    #load columns from tab-spaced data file, e.g. as outputted by a Boris simulation
    def Get_Data_Columns(self, fileName, column_indexes = ''):
        
        #Get data locally
        f = open(fileName, 'r')
        rows = [[float(number) for number in row.rstrip().split('\t')] for row in f.readlines() if self.can_convert(row)]
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
        
        plt.axes(xlabel = xlabel, ylabel = ylabel, title = title)
        plt.grid()
        plt.plot(x, y, label = label_)
        if len(label_): plt.legend()
        
        if len(imageFile): plt.savefig(imageFile + '.png')
        
        plt.show()
        
    def PlotPolar_Data(self, r, theta_deg, xlabel = '', ylabel = '', title = '', label_ = '', imageFile = ''):
        
        #plt.axes(xlabel = xlabel, ylabel = ylabel, title = title)
        plt.grid()
        plt.polar([np.radians(t) for t in theta_deg], r, label = label_)
        if len(label_): plt.legend()
        
        if len(imageFile): plt.savefig(imageFile + '.png')
        
        plt.show()

    #################### OVF2 HELPERS #######################

    #write an OVF2 file for a mesh with given rectangle (m), number of nodes and values in vec list ordered by x, then y, finally z.
    #rect_m must be a lsit with 6 elements : [xmin, ymin, zmin, xmax, ymax, zmax]
    #nodes must be a list with 3 integers : [xnodes, ynodes, znodes]; the cellsize is determined from rect_m and nodes
    #vec can be a scalar quantity (list of floats), or a vector quantity (list of 3-element lists)
    def Write_OVF2(self, fileName, vec, nodes, rect_m):
        
        lines = []
        
        vector = isinstance(vec[0], list) and len(vec[0]) == 3
        
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
        if vector: lines.append("# valuedim: 3")
        else: lines.append("# valuedim: 1")
        lines.append("#")
        lines.append("# End: Header")
        lines.append("#")
        lines.append("# Begin: data binary 8")
               
        with open(fileName, 'wb') as f:
            
            #write header
            for line in lines: f.write(bytes(line + '\n', 'utf-8'))
            
            #check value
            vec.insert(0, 123456789012345.0)
            #write vec_sca as 8-byte floats (including the check value)
            if vector:
                #write vector quantity
                for value in vec:
                    float_aray = np.array(value, 'float64')
                    float_aray.tofile(f)
            else:
                #write scalar quantity
                float_aray = np.array(vec, 'float64')
                float_aray.tofile(f)
            
            #write termination
            f.write(bytes("# End: data " + "binary 8" + '\n', 'utf-8'))
            f.write(bytes("# End: Segment" + '\n', 'utf-8'))

            return True
        
        return False

    #################### SPECIAL COMMANDS #######################

    #Send run command and wait for simulation to finish : blocking call
    def Run(self):

        self.SendCommand("run")
        self.WaitForResponse()
        
    #Save in given filename a new row containing parameters in dataList as tab-spaced characters
    #This uses the savecomment command to save in the local Boris data directory as currently configured in Boris
    def SaveDataToFile(self, fileName, dataList):

        command = 'savecomment ' + fileName + ' '

        for entry in dataList: command += str(entry) + '\t'
        
        #send command with parameters (remove tab if last character - not needed)
        self.SendCommand(command.rstrip('\t'))

    #################### CONSOLE COMMANDS <-> METHODS #######################

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
    
    def cuda(self, status = ''):
    	return self.SendCommand("cuda", [status])
    
    def curietemperature(self, curie_temperature = '', meshname = ''):
    	return self.SendCommand("curietemperature", [curie_temperature, meshname])
    
    def data(self):
    	return self.SendCommand("data")
    
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
    
    def display(self, name = '', meshname = ''):
    	return self.SendCommand("display", [name, meshname])
    
    def displaybackground(self, name = '', meshname = ''):
    	return self.SendCommand("displaybackground", [name, meshname])
    
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
    
    def dp_calcexchange(self):
    	return self.SendCommand("dp_calcexchange")
    
    def dp_calcsot(self, hm_mesh = '', fm_mesh = ''):
    	return self.SendCommand("dp_calcsot", [hm_mesh, fm_mesh])
    
    def dp_calctopochargedensity(self):
    	return self.SendCommand("dp_calctopochargedensity")
    
    def dp_cartesiantopolar(self, dp_in_x = '', dp_in_y = '', dp_out_r = '', dp_out_theta = ''):
    	return self.SendCommand("dp_cartesiantopolar", [dp_in_x, dp_in_y, dp_out_r, dp_out_theta])
    
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
    
    def dp_getpath(self, dp_index_in = '', dp_index_out = ''):
    	return self.SendCommand("dp_getpath", [dp_index_in, dp_index_out])
    
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
    
    def dp_mean(self, dp_index = ''):
    	return self.SendCommand("dp_mean", [dp_index])
    
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
    
    def dp_saveasrow(self, filename = '', dp_index = ''):
    	return self.SendCommand("dp_saveasrow", [filename, dp_index])
    
    def dp_saveappendasrow(self, filename = '', dp_index = ''):
    	return self.SendCommand("dp_saveappendasrow", [filename, dp_index])
    
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
    
    def evalspeedup(self, status = ''):
    	return self.SendCommand("evalspeedup", [status])
    
    def exchangecoupledmeshes(self, status = '', meshname = ''):
    	return self.SendCommand("exchangecoupledmeshes", [status, meshname])
    
    def excludemulticonvdemag(self, status = '', meshname = ''):
    	return self.SendCommand("excludemulticonvdemag", [status, meshname])
    
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
    
    def loadovf2disp(self, filename = ''):
    	return self.SendCommand("loadovf2disp", [filename])
    
    def loadovf2mag(self, renormalize_value = '', filename = ''):
    	return self.SendCommand("loadovf2mag", [renormalize_value, filename])
    
    def loadovf2mesh(self, renormalize_value = '', filename = ''):
    	return self.SendCommand("loadovf2mesh", [renormalize_value, filename])
    
    def loadovf2strain(self, filename_diag = '', filename_odiag = ''):
    	return self.SendCommand("loadovf2strain", [filename_diag, filename_odiag])
    
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
    
    def mcellsize(self, value = ''):
    	return self.SendCommand("mcellsize", [value])
    
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
    
    def ode(self):
    	return self.SendCommand("ode")
    
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
    
    def random(self, meshname = ''):
    	return self.SendCommand("random", [meshname])
    
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
    
    def roughenmesh(self, depth = '', axis = '', seed = ''):
    	return self.SendCommand("roughenmesh", [depth, axis, seed])
    
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
    
    def setdata(self, dataname = '', meshname = '', rectangle = ''):
    	return self.SendCommand("setdata", [dataname, meshname, rectangle])
    
    def setdefaultelectrodes(self):
    	return self.SendCommand("setdefaultelectrodes")
    
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
    
    def setmaterial(self, name = '', rectangle = ''):
    	return self.SendCommand("setmaterial", [name, rectangle])
    
    def setmesh(self, name = '', rectangle = ''):
    	return self.SendCommand("setmesh", [name, rectangle])
    
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
    

    #################### LEGACY METHODS : DON'T USE (directly) #######################

    #Legacy method / Auxiliary : send named command together with any parameters specified using a list
    #New version use methods named after command name directly : improves simulation script readability
    #You can also use this to send the run command as a non-blocking call
    def SendCommand(self, command, values = None):

        try: 
            self.sock.settimeout(self.timeout_ms / 1000)
        except: 
            self.sock.close()
            return

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
                if self.verbose == True: self.sock.sendall(bytes('>' + message, 'utf-8'))
                else: self.sock.sendall(bytes('*' + message, 'utf-8'))
                print('TX : %s' % message)
            except:
                print("SendCommand (send): timed out.")
                self.sock.close()
                return

            # Look for the response
            try:
                data = str(self.sock.recv(self.maxLenMessage), 'utf-8')
            except:
                print("SendCommand (receive): timed out.")
                self.sock.close()
                return
    
            #note, the returned data always starts with a tab
            fields = data.split('\t')

            if len(fields) >= 2 and fields[1] == 'stopped':
                    #if we received the 'stopped' message this means a simulation was running when we sent this command. 
                    #this caused the simulation to stop thus issuing the 'stopped' message since a client is connected.
                    #since this command is expecting another message to be returned, we need to receive it - issue recv call again
                    try:
                        data = str(self.sock.recv(self.maxLenMessage), 'utf-8')
                    except:
                        print("SendCommand (receive): timed out.")
                        self.sock.close()
                        return

                    #note, the returned data always starts with a tab
                    fields = data.split('\t')

            print('RX : %s' % data)

            #the received message should always have at least 2 fields since the message always starts with a tab
            if len(fields) >= 2:      

                #list of floats where there should be floats instead of strings (skip first entry always as this is empty)
                return_data = [self.Convert_Returned_Parameter(entry) for entry in fields[1:]]
                #if the returned data has multiple parameters then return it as a list, else as a single element
                if len(return_data) == 1: return return_data[0]
                else: return return_data

    #Legacy : send multiple commands without parameters
    def SendCommands(self, commands):

        for command in commands: self.SendCommand(command)

    #Legacy : use Run() instead.
    def IsSimulationRunning(self, simRunningPollInterval_ms = 500):

        if abs(time.clock()*1000 - self.pollTimer_ms) > simRunningPollInterval_ms:
            self.pollTimer_ms = time.clock()*1000
            return self.SendCommand('isrunning')
        else:
            #this method will typically be used in a while loop, so make sure user won't flood cpu with function calls
            time.sleep(simRunningPollInterval_ms/10000)

        return True
    
###############################################################################################################################
    
#Still keep this function for legacy compatibility : instead you should just index returned list directly if you're expecting a list
def Get(expected_list, index):

    try: return expected_list[index]
    except: return expected_list