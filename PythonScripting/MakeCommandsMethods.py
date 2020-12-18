"""
Process commands list dump ready for WinSocks.py methods list
"""

#commands list dump here
f_commandslist = open('commands.txt', 'r')

#output formatted methods ready to cut and paste in WinSocks.py here
f_methods = open('methods.txt', 'w')

lines = [line for line in f_commandslist.readlines()]

spec_start = "USAGE : "

#Structure of command methods (use this to generate them programatically after grabbing commands list with their USAGE from Boris):
#Thus if you add new commands in Boris you can just run a separate script to update this module automatically so you don't have to keep track of changes
#def name(self, param1 = '', param2 = '', ...):
#   return self.SendCommand("name", [param1, param2, ...]) 

for line in lines:
    
    if line[:len(spec_start)] == spec_start:
        
        command = line[len(spec_start):].rstrip('\n')
        
        if (len(command)):
            
            #fields are space separated and of the form command param1 param2 ...
            fields = command.split(' ')
            
            #the command name. there may not be any parameters
            command_name = fields[0]
            
            #exception for the run command : there is a special Run method already defined which is blocking
            #don't want to add a non-blocking run method as some users could inadvertently use it and expect it to be blocking
            if command_name == "run": continue
            
            params_list = []
            for param in fields[1:]:
                
                #some optional parameters may be specified as (param) - get rid of all outer brackets
                param = param.lstrip('(').rstrip(')')
                
                #some parameters may specify the directory name is optional as (directory\)name : get rid of it
                if param[:len("directory/)")] == "directory/)": param = param.lstrip("directory/)")
                
                #if there are multiple levels of optional parameters, e.g. (param1, (param2)), then a comma will be present now - get rid of it
                #if a list of parameters is possible this will be specified as param... - get rid of ... : in the Python method you can pass this as a list to the function
                param = param.rstrip(',').rstrip("...")
                
                #done with parameter : add it
                if len(param): params_list.append(param)
            
            #make method header with or without parameters
            method_header = "def " + command_name + "(self"
            
            for param in params_list:
                method_header += ", " + param + " = ''"
            
            method_header += "):\n"
            
            #write completed header
            f_methods.write(method_header)
            
            #make method code
            method_code = "\treturn self.SendCommand(" + '"' + command_name + '"'
            
            if len(params_list): method_code += ", ["
            
            for param in params_list:
                method_code += param + ', '
                
            if len(params_list):
                method_code = method_code.rstrip(', ') + "]"
                
            method_code += ")\n\n"
            
            f_methods.write(method_code)
        
f_methods.close()
f_commandslist.close()