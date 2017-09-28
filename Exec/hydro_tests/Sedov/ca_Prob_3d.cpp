/*
Compile me using

g++ -I/home/alice/anaconda3/include/python3.6m cpp_code.cpp -o cpp_code -L/home/alice/anaconda3/lib -lpython3.6m

or:

g++ `python3.6-config --cflags` ca_Prob_3d.cpp -o cpp_code `python3.6-config --ldflags` -I/home/alice/Documents/amrex/Src/Base -I/home/alice/Documents/amrex/Src/AmrCore -I/home/alice/Documents/amrex/Src/Amr

*/

#include "Python.h"
#include <AMReX.H>
#include <AMReX_Utility.H>
#include "Castro_F.H"
#include "Castro_io.H"
#include <AMReX_ParmParse.H>
#include <Castro.H>
#include "AMReX_buildInfo.H"
#include <typeinfo>


using namespace amrex;

//int main() {
//    return 0;
//}

void Castro::ca_initdata(int& level, amrex::Real& time,
                        const int* lo, const int* hi,
                        int& num_state,
                        double* state, const int* slo, const int* shi,
                        const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* xhi)
{

    //std::cout << "Calling ca_initdata\n";

    // initialise python interpreter
    Py_Initialize();

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

    //while (pValue != NULL) Py_DECREF(pValue);

    // build the module object
    PyObject *pModule = PyImport_ImportModule(pymodule);

    // call this again to set up desc_lst as PyImport_ImportModule
    // manages to destroy all the global variables?

    amrex::ParmParse pp("amr");
    pp.query("probin_file", probin_file);
    variableReSetUp();

    /*std::cout << "After import module in Castro::ca_initdata, state variables are:\n";
    std::cout <<  "size = " << get_desc_lst().size() << '\n';
    for (int typ = 0; typ < get_desc_lst().size(); typ++)
    {
        const amrex::StateDescriptor& desc = get_desc_lst()[typ];

        for (int n = 0; n < desc.nComp(); n++)
        {
            std::cout << desc.name(n) << '\n';
        }
    }*/

    if (pModule != NULL) {
        PyObject *pDict = PyModule_GetDict(pModule);
        PyObject *pFunc = PyDict_GetItemString(pDict, pyfunc);
        /* pFunc is a new reference */

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

            //std::cout << "Called python object\n";

            int NVAR;
            ca_get_nvar(&NVAR);

            //std::cout << "size of state = " << sizeof(state)/sizeof(*state) << '\n';

            for (int i = 0; i < (shi[0]+1)*(shi[1]+1)*(shi[2]+1)*NVAR; i++) {

                pValue = PyList_GetItem(pstate, i); // doesn't INCREF pValue

                if (pValue != NULL)
                    state[i] = PyFloat_AsDouble(pValue);
                //std::cout << "Made double at " << i << " = " << PyFloat_AsDouble(pValue) <<'\n';
            }
            //Py_DECREF(pValue);; //borrowed ref from pstate
            Py_DECREF(pstate);

            //std::cout << "Made state\n";

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

        //Py_XDECREF(pFunc); //borrowed ref
        //Py_XDECREF(pDict); borrowed ref
    }
    else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", pymodule);

    }

    //Py_DECREF(plo); // pArgs will decref these
    //Py_DECREF(phi);
    //Py_DECREF(pslo);
    //Py_DECREF(pshi);
    //Py_DECREF(pdx);
    //Py_DECREF(pxlo);
    //Py_DECREF(pxhi);
    Py_DECREF(pModule);

    //Py_DECREF(pValue);
    // restore previous GIL state and return
    PyGILState_Release(gil_state);

    Py_Finalize();

    //std::cout <<  "exiting\n";
}
