/*
Compile me using

g++ -I/home/alice/anaconda3/include/python3.6m cpp_code.cpp -o cpp_code -L/home/alice/anaconda3/lib -lpython3.6m

or:

g++ `python3.6-config --cflags` ca_Prob_3d.cpp -o cpp_code `python3.6-config --ldflags` -I/home/alice/Documents/amrex/Src/Base -I/home/alice/Documents/amrex/Src/AmrCore -I/home/alice/Documents/amrex/Src/Amr

*/

#include "Python.h"
#include <AMReX.H>
#include <AMReX_Utility.H>
//#include <AMReX_CONSTANTS.H>
#include "Castro_F.H"
#include "Castro_io.H"
#include <AMReX_ParmParse.H>
#include <Castro.H>
#include "AMReX_buildInfo.H"
//#include "AMReX_IntVect.H"


using namespace amrex;

//int main() {
//    return 0;
//}

void Castro::ca_initdata(const int& level, const double& time,
                const int* lo, const int* hi,
                const int& num_state,
                double* state, const int* slo, const int* shi,
                const double* dx, const double* xlo, const double* xhi)
{

    PyObject *pName, *pModule, *pDict, *pFunc;
    PyObject *pArgs, *pValue;
    PyObject *plo, *phi, *pslo, *pshi, *pdx, *pxlo, *pxhi, *pstate;

    const char* pymodule = "Prob_3d";
    const char* pyfunc = "ca_initdata";
    // initialise python interpreter
    Py_Initialize();
    // add current directory to python path
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\".\")");

    int dim = 3;

    plo = PyTuple_New(dim);
    phi = PyTuple_New(dim);
    pslo = PyTuple_New(dim);
    pshi = PyTuple_New(dim);
    pdx = PyTuple_New(dim);
    pxlo = PyTuple_New(dim);
    pxhi = PyTuple_New(dim);

    for (int i = 0; i < dim; i++) {
        pValue = PyLong_FromLong(lo[i]);
        PyTuple_SetItem(plo, i, pValue);

        pValue = PyLong_FromLong(hi[i]);
        PyTuple_SetItem(phi, i, pValue);

        pValue = PyLong_FromLong(slo[i]);
        PyTuple_SetItem(pslo, i, pValue);

        pValue = PyLong_FromLong(shi[i]);
        PyTuple_SetItem(pshi, i, pValue);

        pValue = PyFloat_FromDouble(dx[i]);
        PyTuple_SetItem(pdx, i, pValue);

        pValue = PyLong_FromLong(xlo[i]);
        PyTuple_SetItem(pxlo, i, pValue);

        pValue = PyLong_FromLong(xhi[i]);
        PyTuple_SetItem(pxhi, i, pValue);
    }

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

            pArgs = PyTuple_New(7);

            PyTuple_SetItem(pArgs, 0, plo);
            PyTuple_SetItem(pArgs, 1, phi);
            PyTuple_SetItem(pArgs, 2, pslo);
            PyTuple_SetItem(pArgs, 3, pshi);
            PyTuple_SetItem(pArgs, 4, pdx);
            PyTuple_SetItem(pArgs, 5, pxlo);
            PyTuple_SetItem(pArgs, 6, pxhi);

            pstate = PyObject_CallObject(pFunc, pArgs);
            Py_DECREF(pArgs);

            //pstate = PyObject_GetItem(pValue, 0);

            // TODO: need NQ, NVAR here
            for (int i = 0; i < (shi[0]-slo[0])*(shi[1]-slo[1])*(shi[2]-slo[2]); i++) {

                pValue = PyList_GetItem(pstate, i);
                state[i] = PyFloat_AsDouble(pValue);
            }
            Py_DECREF(pValue);

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
