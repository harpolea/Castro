/*
Compile me using

g++ -I/home/alice/anaconda3/include/python3.6m cpp_code.cpp -o cpp_code -L/home/alice/anaconda3/lib -lpython3.6m

or:

g++ `python3.6-config --cflags` cpp_code.cpp -o cpp_code `python3.6-config --ldflags`

*/

#include <AMReX.H>
#include "Python.h"

int main() {
    return 0;
}

int ca_initdata(const int& level, const double& time,
                const int* lo, const int* hi,
                const int& num_state,
                double* state, const int* slo, const int* shi,
                const double* dx, const double* xlo, const double* xhi)
{

    PyObject *pName, *pModule, *pDict, *pFunc;
    PyObject *pArgs, *pValue;
    //PyArrayObject *pyState, *pyLo, *pyHi, *pyState, *pyDx, *pyXlo, *pyXhi;

    const char* pymodule = "Prob_3d";
    const char* pyfunc = "ca_initdata";
    // initialise python interpreter
    Py_Initialize();
    // add current directory to python path
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\".\")");

    // build the name object
    pName = PyUnicode_FromString(pymodule);
    /* Error checking of pName left out */

    // build the module object
    pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule != NULL) {
        pDict = PyModule_GetDict(pModule);
        pFunc = PyDict_GetItemString(pDict, pyfunc);
        /* pFunc is a new reference */

        if (pFunc && PyCallable_Check(pFunc)) {
            PyArg_ParseTuple(pArgs, "ifiiffff", &level, &time, lo, hi, state, dx, xlo, xhi);

            PyObject_CallObject(pFunc, pArgs);
            Py_DECREF(pArgs);

        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\"\n", pyfunc);
        }
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
    }
    else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", pymodule);

    }

    Py_Finalize();
}

/*double *pyvector_to_Carrayptrs(PyArrayObject *arrayin)  {
  int n = arrayin->dimensions[0];
  return (double *) arrayin->data;  // pointer to arrayin data as double
}

int *intpyvector_to_Carrayptrs(PyArrayObject *arrayin)  {
  int n = arrayin->dimensions[0];
  return (int *) arrayin->data;  // pointer to arrayin data as double
}*/
