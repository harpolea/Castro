import numpy as np
from eos import eos_type_module, eos_module
from riemann_util import riemann_util_module as riemann
from probdata import probdata_module as probdata
from prob_params import prob_params_module as prob_params
from meth_params import meth_params_module as meth_params

def amrex_probinit (probin):# bind(c)

    """use bl_constants_module
    use probdata_module
    use prob_params_module, only : center
    use eos_type_module, only: eos_t, eos_input, eos_input_rp
    use eos_module, only: eos
"""

    # set namelist defaults

    probdata.p_ambient = 1.e-5        # ambient pressure (in erg/cc)
    probdata.dens_ambient = 1.e0      # ambient density (in g/cc)
    probdata.exp_energy = 1.e0        # absolute energy of the explosion (in erg)
    probdata.r_init = 0.05e0          # initial radius of the explosion (in cm)
    probdata.nsub = 4
    probdata.temp_ambient = -1.e2     # Set original temp. to negative, which is overwritten in the probin file

    # set local variable defaults
    prob_params.center = np.zeros(3)
    prob_params.center[0] = 0.5*(problo[0] + probhi[0])
    prob_params.center[1] = 0.5*(problo[1] + probhi[1])
    prob_params.center[2] = 0.5*(problo[2] + probhi[2])

    # Read namelists
    params = read_probin(probin)

    extract_dict(params, 'p_ambient', probdata.p_ambient)
    extract_dict(params, 'dens_ambient', probdata.dens_ambient)
    extract_dict(params, 'temp_ambient', probdata.temp_ambient)
    extract_dict(params, 'exp_energy', probdata.exp_energy)
    extract_dict(params, 'r_init', probdata.r_init)
    extract_dict(params, 'nsub', probdata.nsub)

    probdata.xn_zone[:] = 0.0
    probdata.xn_zone[0] = 1.0

    eos_state = eos_type_module.Eos_T()

    # override the pressure iwth the temperature
    if (temp_ambient > 0.0):
        eos_state.rho = probdata.dens_ambient
        eos_state.xn[:] = probdata.xn_zone[:]
        eos_state.T = probdata.temp_ambient

        eos_module.eos(eos_type_module.eos_input_rt, eos_state)

        probdata.p_ambient = eos_state.p


    # Calculate ambient state data

    eos_state.rho = probdata.dens_ambient
    eos_state.p   = probdata.p_ambient
    eos_state.T   = 1.e5 # Initial guess for iterations
    eos_state.xn  = probdata.xn_zone

    eos_module.eos(eos_type_module.eos_input_rp, eos_state)

    probdata.e_ambient = eos_state.e


# ::: -----------------------------------------------------------
# ::: This routine is called at problem setup time and is used
# ::: to initialize data on each grid.
# :::
# ::: NOTE:  all arrays have one cell of ghost zones surrounding
# :::        the grid interior.  Values in these cells need not
# :::        be set here.
# :::
# ::: INPUTS/OUTPUTS:
# :::
# ::: level     => amr level of grid
# ::: time      => time at which to init data
# ::: lo,hi     => index limits of grid interior (cell centered)
# ::: nstate    => number of state components.  You should know
# :::		   this already#
# ::: state     <=  Scalar array
# ::: delta     => cell size
# ::: xlo,xhi   => physical locations of lower left and upper
# :::              right hand corner of grid.  (does not include
# :::		   ghost region).
# ::: -----------------------------------------------------------
def ca_initdata(lo, hi, slo, shi, delta, xlo, xhi):

    """use probdata_module
    use meth_params_module , only: NVAR, NQ, QRHO, QU, QV, QW, QREINT, QPRES, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
    use prob_params_module, only : center
    use network, only : nspec
    use eos_module, only : eos
    use eos_type_module, only : eos_t, eos_input_rp, eos_input_re
    use riemann_util_module, only: gr_cons_state
"""

    gamma_up = np.zeros(9)
    gamma_up[0] = 1.0
    gamma_up[4] = 1.0
    gamma_up[8] = 1.0

    # set explosion pressure -- we will convert the point-explosion energy into
    # a corresponding pressure distributed throughout the perturbed volume
    vctr  = 4.0/3.0 * np.pi * probdata.r_init**3

    e_zone = probdata.exp_energy / vctr / probdata.dens_ambient

    eos_state = eos_type_module.Eos_T()

    eos_state.e = e_zone
    eos_state.rho = probdata.dens_ambient
    eos_state.xn[:] = probdata.xn_zone[:]
    eos_state.T = 100.0  # initial guess

    eos_module.eos(eos_type_module.eos_input_re, eos_state)

    p_exp = eos_state.p

    NQ = meth_params.nq
    NVAR = meth_params.nvar

    state = np.zeros((shi[0]-slo[0]+1, shi[1]-slo[1]+1, shi[2]-slo[2]+1, NVAR))
    q = np.zeros((shi[0]-slo[0]+1, shi[1]-slo[1]+1, shi[2]-slo[2]+1, NQ))

    QRHO = meth_params.qrho
    QU = meth_params.qu
    QV = meth_params.qv
    QW = meth_params.qw
    QREINT = meth_params.qreint
    QPRES = meth_params.qpres
    UTEMP = meth_params.utemp
    UFS = meth_params.ufs
    URHO = meth_params.urho

    for k in range(lo[2], hi[2]):
        zmin = xlo[2] + delta[2]*(k-lo[2])

        for j in range(lo[1], hi[1]):
            ymin = xlo[1] + delta[1]*(j-lo[1])

            for i in range(lo[0], hi[0]):
                xmin = xlo[0] + delta[0]*(i-lo[0])

                npert = 0
                nambient = 0

                for kk in range(probdata.nsub):
                    zz = zmin + (delta[2]/probdata.nsub)*(kk + 0.5e0)

                    for jj in range(probdata.nsub):
                        yy = ymin + (delta[1]/probdata.nsub)*(jj + 0.5e0)

                        for ii in range(nsub):
                            xx = xmin + (delta[0]/probdata.nsub)*(ii + 0.5e0)

                            dist = (prob_params.center[0]-xx)**2 + (prob_params.center[1]-yy)**2 + (prob_params.center[2]-zz)**2

                            if(dist <= r_init**2):
                                npert = npert + 1
                            else:
                                nambient = nambient + 1


                p_zone = (npert*probdata.p_exp + nambient*probdata.p_ambient) / probdata.nsub**3

                eos_state.p = p_zone
                eos_state.rho = probdata.dens_ambient
                eos_state.xn[:] = probdata.xn_zone[:]

                eos_module.eos(eos_type_module.eos_input_rp, eos_state)

                eint = probdata.dens_ambient * eos_state.e

                q[i,j,k,QRHO] = probdata.dens_ambient
                q[i,j,k,QU] = 0.e0
                q[i,j,k,QV] = 0.e0
                q[i,j,k,QW] = 0.e0

                q[i,j,k,QREINT] = eint
                q[i,j,k,QPRES] = eos_state.p

                riemann.gr_cons_state(q[i,j,k,:], state[i,j,k,:], gamma_up)

                state[i,j,k,UTEMP] = eos_state.T

                state[i,j,k,UFS] = state[i,j,k,URHO]

    return np.ascontiguousarray(state)

def read_probin(filename):
    try:
        f = open(filename)
    except:
        sys.exit("error opening the input file")

    lines = f.readlines()
    f.close()
    params = {}

    # first check formatting and create dictionary
    for i, l in enumerate(lines):
        if l[0] == '&' or l[0] == '/' or re.match('\s+', l): #comment
            continue
        else:
            m = re.match('([\w\.]+)\s*\=\s*([\w\s\-\."]+)', l)
            if m is None:
                print('Line {} is not of the expected format: {}'.format(i+1, l))
            else:
                params[m.group(1)] = m.group(2).rstrip().split()
                values = params[m.group(1)].copy()

                try:
                    for j, v in enumerate(params[m.group(1)]):
                        values[j] = int(v)
                    params[m.group(1)] = values
                except ValueError: # try and see if they're floats
                    try:
                        for j, v in enumerate(params[m.group(1)]):
                            values[j] = float(v)
                        params[m.group(1)] = values
                    except ValueError:
                        # pass them all as strings instread
                        value = values
                if len(params[m.group(1)]) == 1:
                    params[m.group(1)] = params[m.group(1)][0]

    return params

def extract_dict(dictionary, param, variable):
    try:
        variable = dictionary[param]
    except KeyError:
        pass
