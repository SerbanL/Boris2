#pragma once

#include "Boris.h"

#if PYTHON_EMBEDDING == 1

#define PY_SSIZE_T_CLEAN
#include <Python.h>

static PyObject* BPython_Message(PyObject *self, PyObject *args)
{
	char *pinput_string;
	if (!PyArg_ParseTuple(args, "s", &pinput_string)) return nullptr;

	std::string input_string(pinput_string);
	if (pSim) {

		std::string response = pSim->NewPythonMessage(input_string);
		return PyUnicode_FromString(response.c_str());
	}
	else return nullptr;
}

static PyObject* BPython_Message_Run(PyObject *self, PyObject *args)
{
	if (pSim) {

		pSim->NewPythonMessage_Run();
		return PyUnicode_FromString("");
	}
	else return nullptr;
}

static PyObject* BPython_Message_RunStage(PyObject *self, PyObject *args)
{
	int stage;
	if (!PyArg_ParseTuple(args, "i", &stage)) return nullptr;

	if (pSim) {

		pSim->NewPythonMessage_RunStage(stage);
		return PyUnicode_FromString("");
	}
	else return nullptr;
}

static PyObject* BPython_Message_RunWithCode(PyObject *self, PyObject *args)
{
	if (pSim) {

		int running = pSim->NewPythonMessage_RunWithCode();
		return PyLong_FromLong(running);
	}
	else return nullptr;
}

static PyObject* BPython_Message_PrepareRunWithCode(PyObject *self, PyObject *args)
{
	if (pSim) {

		pSim->NewPythonMessage_PrepareRunWithCode();
		return PyUnicode_FromString("");
	}
	else return nullptr;
}

static PyObject* BPython_StopPythonScriptFlag(PyObject *self, PyObject *args)
{
	if (pSim) return PyLong_FromLong(pSim->StopPythonScriptFlag());
	else return nullptr;
}

static PyObject* sys_stdout_redirect(PyObject *self, PyObject *args)
{
	const char *stringchars;
	if (!PyArg_ParseTuple(args, "s", &stringchars)) return nullptr;
	//after print statements, stdout issues a "\n" message, but DisplayConsoleMessage already adds a new line, so this simply results in an additional blank line after every print.
	//Stopping single "\n" messages does mean print("\n") will have no effect however.
	std::string text = std::string(stringchars);
	if (pSim && text != "\n") pSim->DisplayConsoleListing(text);
	
	Py_INCREF(Py_None);
	return Py_None;
}

static PyMethodDef PythonEmbeddedMethods[] = {

	{"BPython_Message", BPython_Message, METH_VARARGS,
	 "Receive and execute a message from embedded Python script, and return response."},
	 {"BPython_Message_Run", BPython_Message_Run, METH_VARARGS,
	 "Execute run command on simulation from embedded Python script."},
	 {"BPython_Message_RunStage", BPython_Message_RunStage, METH_VARARGS,
	 "Execute run command on simulation from embedded Python script for given stage number only (i.e. this is same as the runstage command)"},
	 {"BPython_Message_RunWithCode", BPython_Message_RunWithCode, METH_VARARGS,
	 "Execute run command on simulation from embedded Python script, in a way which allows further Python code to be executed after every iteration."},
	 {"BPython_Message_PrepareRunWithCode", BPython_Message_PrepareRunWithCode, METH_VARARGS,
	 "Execute run command on simulation from embedded Python script, in a way which allows further Python code to be executed after every iteration - first call this to prepare."},
	 {"BPython_StopPythonScriptFlag", BPython_StopPythonScriptFlag, METH_VARARGS,
	 "Embedded Python script will poll this flag to see if it needs to terminate execution."},
	 {"sys_stdout_redirect", sys_stdout_redirect, METH_VARARGS,
	 "sys stdout redirection so we can display messages in Boris console."},
	{nullptr, nullptr, 0, nullptr}
};

static PyModuleDef PythonEmbeddedModule = {

	PyModuleDef_HEAD_INIT, "BPython", nullptr, -1, PythonEmbeddedMethods,
	nullptr, nullptr, nullptr, nullptr
};

static PyObject* PyInit_PythonEmbeddedModule(void)
{
	return PyModule_Create(&PythonEmbeddedModule);
}

#endif


