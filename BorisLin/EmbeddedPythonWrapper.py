######################################################
# 
# This wrapper is needed so we can stop the embedded Python script early when needed.
# Also needed to redirect sys.stdout to Boris.
#
# NOTE : If you normally use an IDE to run Python, some required packages will likely not be available for embedded Python.
# In this case you need to install them using pip. e.g. from console run this (make sure Python is added to Path environment variable):
#
# python -m pip install --upgrade pip
# pip install numpy
# pip install matplotlib
#
# etc. 
# 
# Packages will then be available for embedded Python scripts to import in the usual way.
#
######################################################

import BPython, threading, ctypes, time, sys

#Redirect sys.stdout to Boris console (e.g. display Python print output in Boris console)
class StdoutCatcher:
    def write(self, text):
        BPython.sys_stdout_redirect(text)
sys.stdout = StdoutCatcher()

#Wrapper for user script
class thread_with_exception(threading.Thread):
    def __init__(self, name):
        threading.Thread.__init__(self)
        self.name = name
             
    def run(self):
        try:
            #user script goes here in its entirety - do not modify this line
            
        finally: pass
          
    def get_id(self):
        if hasattr(self, '_thread_id'): return self._thread_id
        for id, thread in threading._active.items():
            if thread is self: return id
  
    def raise_exception(self):
        thread_id = self.get_id()
        res = ctypes.pythonapi.PyThreadState_SetAsyncExc(thread_id, ctypes.py_object(SystemExit))
        if res > 1: ctypes.pythonapi.PyThreadState_SetAsyncExc(thread_id, 0)
      
#Run user script on a thread which catches exceptions
BPython_Thread = thread_with_exception('BPython Thread')
BPython_Thread.start()

#To stop the embedded Python script (using the stopscript command in Boris), then raise exception
while True:
    StopPythonScriptFlag = BPython.BPython_StopPythonScriptFlag()
    if StopPythonScriptFlag == 1 or not BPython_Thread.is_alive():
        BPython_Thread.raise_exception()
        BPython_Thread.join()
        break
    time.sleep(0.1)