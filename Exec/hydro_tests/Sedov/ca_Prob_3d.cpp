#include "Python.h"
#include <AMReX.H>
#include <AMReX_Utility.H>
#include "Castro_F.H"
#include "Castro_io.H"
#include <AMReX_ParmParse.H>
#include <Castro.H>
#include "AMReX_buildInfo.H"
#include <typeinfo>
#include <AMReX_PROB_AMR_F.H>

using namespace amrex;

void amrex_probinit (const int* init,
         const int* name,
         const int* namelen,
         const amrex_real* problo,
         const amrex_real* probhi) {

    Py_BEGIN_ALLOW_THREADS;

    // Make sure own the GIL
    PyGILState_STATE gil_state = PyGILState_Ensure();

    const char* pymodule = "Prob_3d";
    const char* pyfunc = "amrex_probinit";

    // add current directory to python path
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\".\")");

    int dim = 3;

    PyObject *pproblo = PyTuple_New(dim);
    PyObject *pprobhi = PyTuple_New(dim);

    for (int i = 0; i < dim; i++) {
        PyTuple_SetItem(pproblo, i, PyFloat_FromDouble(problo[i]));
        PyTuple_SetItem(pprobhi, i, PyFloat_FromDouble(probhi[i]));
    }

    char probin[*namelen];

    for (int i = 0; i < *namelen; i++) {
        probin[i] = char(name[i]);
    }

    PyObject *pprobin = PyUnicode_FromString(probin);

    // build the module object
    PyObject *pModule = PyImport_ImportModule(pymodule);

    if (pModule != NULL) {
        PyObject *pDict = PyModule_GetDict(pModule);
        PyObject *pFunc = PyDict_GetItemString(pDict, pyfunc);

        if (pFunc && PyCallable_Check(pFunc)) {

            PyObject *pArgs = PyTuple_New(3);

            PyTuple_SetItem(pArgs, 0, pproblo);
            PyTuple_SetItem(pArgs, 1, pprobhi);
            PyTuple_SetItem(pArgs, 2, pprobin);

            PyObject *pstate = PyObject_CallObject(pFunc, pArgs);
            Py_DECREF(pArgs);

        } else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\"\n", pyfunc);

            Py_DECREF(pproblo);
            Py_DECREF(pprobhi);
            Py_DECREF(pprobin);
        }

    } else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", pymodule);

    }

    Py_DECREF(pModule);

    // restore previous GIL state and return
    PyGILState_Release(gil_state);

    Py_END_ALLOW_THREADS;

}

void Castro::ca_initdata(int& level, amrex::Real& time,
                        const int* lo, const int* hi,
                        int& num_state,
                        double* state, const int* slo, const int* shi,
                        const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* xhi)
{

    PyObject *pValue;

    const char* pymodule = "Prob_3d";
    const char* pyfunc = "ca_initdata";

    // Make sure own the GIL
    PyGILState_STATE gil_state = PyGILState_Ensure();
    // add current directory to python path
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\".\")");

    int dim = 3;

    PyObject *plo = PyTuple_New(dim);
    PyObject *phi = PyTuple_New(dim);
    PyObject *pslo = PyTuple_New(dim);
    PyObject *pshi = PyTuple_New(dim);
    PyObject *pdx = PyTuple_New(dim);
    PyObject *pxlo = PyTuple_New(dim);
    PyObject *pxhi = PyTuple_New(dim);

    for (int i = 0; i < dim; i++) {
        PyTuple_SetItem(plo, i, PyLong_FromLong(lo[i]));

        PyTuple_SetItem(phi, i, PyLong_FromLong(hi[i]));

        PyTuple_SetItem(pslo, i, PyLong_FromLong(slo[i]));

        PyTuple_SetItem(pshi, i, PyLong_FromLong(shi[i]));

        PyTuple_SetItem(pdx, i, PyFloat_FromDouble(dx[i]));

        PyTuple_SetItem(pxlo, i, PyLong_FromLong(xlo[i]));

        PyTuple_SetItem(pxhi, i, PyLong_FromLong(xhi[i]));
    }

    // build the module object
    PyObject *pModule = PyImport_ImportModule(pymodule);

    amrex::ParmParse pp("amr");
    pp.query("probin_file", probin_file);
    variableReSetUp();

    if (pModule != NULL) {
        PyObject *pDict = PyModule_GetDict(pModule);
        PyObject *pFunc = PyDict_GetItemString(pDict, pyfunc);

        if (pFunc && PyCallable_Check(pFunc)) {

            PyObject *pArgs = PyTuple_New(7);

            PyTuple_SetItem(pArgs, 0, plo);
            PyTuple_SetItem(pArgs, 1, phi);
            PyTuple_SetItem(pArgs, 2, pslo);
            PyTuple_SetItem(pArgs, 3, pshi);
            PyTuple_SetItem(pArgs, 4, pdx);
            PyTuple_SetItem(pArgs, 5, pxlo);
            PyTuple_SetItem(pArgs, 6, pxhi);

            PyObject *pstate = PyObject_CallObject(pFunc, pArgs);
            Py_DECREF(pArgs);

            int NVAR;
            ca_get_nvar(&NVAR);

            for (int i = 0; i < (shi[0]+1)*(shi[1]+1)*(shi[2]+1)*NVAR; i++) {

                pValue = PyList_GetItem(pstate, i);

                if (pValue != NULL)
                    state[i] = PyFloat_AsDouble(pValue);
            }
            Py_DECREF(pstate);

        } else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\"\n", pyfunc);

            Py_DECREF(plo);
            Py_DECREF(phi);
            Py_DECREF(pslo);
            Py_DECREF(pshi);
            Py_DECREF(pdx);
            Py_DECREF(pxlo);
            Py_DECREF(pxhi);
        }

    }
    else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", pymodule);
    }

    Py_DECREF(pModule);

    PyGILState_Release(gil_state);

    StateData& statedata2 = (*this).get_state_data(0);
    const StateDescriptor& desc2 = *(statedata2.descriptor());
    std::cout << "in Castro::ca_initdata, end of function = " << desc2.nComp() << '\n';
}
