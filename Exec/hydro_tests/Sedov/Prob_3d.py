import numpy as np
from eos import eos_type_module, eos_module
from riemann_util import riemann_util_module as riemann
from probdata import probdata_module as probdata
from prob_params import prob_params_module as prob_params
from meth_params import meth_params_module as meth_params
import re

def amrex_probinit (problo, probhi, probin):# bind(c)

    """use bl_constants_module
    use probdata_module
    use prob_params_module, only : center
    use eos_type_module, only: eos_t, eos_input, eos_input_rp
    use eos_module, only: eos
"""

    print("\nCalling Prob_3d.py:amrex_probinit")

    # set namelist defaults

    probdata.p_ambient = 1.e-5        # ambient pressure (in erg/cc)
    probdata.dens_ambient = 1.e0      # ambient density (in g/cc)
    probdata.exp_energy = 1.e0        # absolute energy of the explosion (in erg)
    probdata.r_init = 0.05e0          # initial radius of the explosion (in cm)
    probdata.nsub = 4
    probdata.temp_ambient = -1.e2     # Set original temp. to negative, which is overwritten in the probin file

    # set local variable defaults
    #prob_params.center = np.zeros(3)
    prob_params.center[0] = 0.5*(problo[0] + probhi[0])
    prob_params.center[1] = 0.5*(problo[1] + probhi[1])
    prob_params.center[2] = 0.5*(problo[2] + probhi[2])

    # Read namelists

    # for some reason an extra garbage character is stuck on the end, so get rid of that
    params = py_read_probin(probin[:-1])

    if 'p_ambient' in params:
        probdata.p_ambient = params['p_ambient']
    if 'dens_ambient' in params:
        probdata.dens_ambient = params['dens_ambient']
    if 'temp_ambient' in params:
        probdata.temp_ambient = params['temp_ambient']
    if 'exp_energy' in params:
        probdata.exp_energy = params['exp_energy']
    if 'r_init' in params:
        probdata.r_init = params['r_init']
    if 'nsub' in params:
        probdata.nsub = params['nsub']

    probdata.xn_zone[:] = 0.0
    probdata.xn_zone[0] = 1.0

    eos_state = eos_type_module.Eos_T()

    if meth_params.qu == 0:
        meth_params.f_set_castro_method_params()

    eos_module.eos_init(1.e-20, 1.e-20)

    # override the pressure with the temperature
    if (probdata.temp_ambient > 0.0):
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

    print("Leaving Prob_3d.py:amrex_probinit\n")


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
def ca_initdata(lo, hi, slo, shi, delta, xlo, xhi, probin="probin.3d.sph"):

    """use probdata_module
    use meth_params_module , only: NVAR, NQ, QRHO, QU, QV, QW, QREINT, QPRES, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
    use prob_params_module, only : center
    use network, only : nspec
    use eos_module, only : eos
    use eos_type_module, only : eos_t, eos_input_rp, eos_input_re
    use riemann_util_module, only: gr_cons_state
"""

    print('Calling Prob_3d.py:ca_initdata')

    if meth_params.qu == 0:
        meth_params.f_set_castro_method_params()

    gamma_up = np.zeros(9)
    gamma_up[0] = 1.0
    gamma_up[4] = 1.0
    gamma_up[8] = 1.0

    xn_zone = np.zeros(2)
    xn_zone[:] = 0.0
    xn_zone[0] = 1.0

    params = py_read_probin(probin)

    p_ambient = params['p_ambient']
    dens_ambient = params['dens_ambient']
    if 'temp_ambient' in params:
        temp_ambient = params['temp_ambient']
    else:
        temp_ambient = 0.
    exp_energy = params['exp_energy']
    r_init = params['r_init']
    nsub = params['nsub']

    # set explosion pressure -- we will convert the point-explosion energy into
    # a corresponding pressure distributed throughout the perturbed volume
    vctr  = 4.0/3.0 * np.pi * r_init**3

    e_zone = exp_energy / vctr / dens_ambient

    eos_state = eos_type_module.Eos_T()
    if eos_module.get_initialized() != 0:
        eos_module.eos_init(1.e-20, 1.e-20)
    eos_state.e = e_zone
    eos_state.rho = dens_ambient
    eos_state.xn = xn_zone[:len(eos_state.xn)]
    eos_state.T = 100.0  # initial guess

    eos_module.eos(eos_type_module.eos_input_re, eos_state)

    p_exp = eos_state.p

    NQ = meth_params.nq
    NVAR = meth_params.nvar

    state = np.zeros((shi[0]+1, shi[1]+1, shi[2]+1, NVAR))
    q = np.zeros((shi[0]+1, shi[1]+1, shi[2]+1, NQ))

    #print(meth_params)

    # minus 1 as stupid fortran array numbering from 1
    QRHO = meth_params.qrho-1
    QU = meth_params.qu-1
    QV = meth_params.qv-1
    QW = meth_params.qw-1
    QREINT = meth_params.qreint-1
    QPRES = meth_params.qpres-1
    UTEMP = meth_params.utemp-1
    UFS = meth_params.ufs-1
    URHO = meth_params.urho-1

    q[:,:,:,QRHO] = dens_ambient

    centre = prob_params.center

    for k in range(lo[2], hi[2]+1):
        zmin = xlo[2] + delta[2]*(k-lo[2])

        for j in range(lo[1], hi[1]+1):
            ymin = xlo[1] + delta[1]*(j-lo[1])

            for i in range(lo[0], hi[0]+1):
                xmin = xlo[0] + delta[0]*(i-lo[0])

                zz = zmin + (delta[2]/nsub)*(np.array(range(nsub)) + 0.5)
                yy = ymin + (delta[1]/nsub)*(np.array(range(nsub)) + 0.5)
                xx = xmin + (delta[0]/nsub)*(np.array(range(nsub)) + 0.5)

                XX, YY, ZZ = np.meshgrid(xx, yy, zz)

                dist = (centre[0]-XX)**2 + (centre[1]-YY)**2 + (centre[2]-ZZ)**2

                npert = len(dist[dist <= r_init**2].flatten())
                nambient = len(dist.flatten()) - npert

                p_zone = (npert*p_exp + nambient*p_ambient) / nsub**3

                eos_state.p = p_zone
                eos_state.rho = dens_ambient
                eos_state.xn = xn_zone[:len(eos_state.xn)]

                eos_module.eos(eos_type_module.eos_input_rp, eos_state)

                eint = dens_ambient * eos_state.e

                q[i,j,k,QREINT] = eint
                q[i,j,k,QPRES] = eos_state.p

                riemann.gr_cons_state(np.asfortranarray(q[i,j,k,:]), np.asfortranarray(state[i,j,k,:]), np.asfortranarray(gamma_up))

                state[i,j,k,UTEMP] = eos_state.T

    q[:,:,:,QU] = 0.0
    q[:,:,:,QV] = 0.0
    q[:,:,:,QW] = 0.0

    if UFS < len(state[lo[0],lo[1],lo[2],:]):
        state[:,:,:,UFS] = state[:,:,:,URHO]

    meth_params.f_finalize_meth_params()

    print("Leaving Prob_3d.py:ca_initdata\n")

    return list(np.ascontiguousarray(np.ndarray.flatten(state.T)))

def py_read_probin(filename):

    params = {}

    try:
        with open(filename, "r") as f:

            lines = f.readlines()

            # first check formatting and create dictionary
            for i, l in enumerate(lines):
                if l[0] == '&' or l[0] == '/' or re.fullmatch('\s+', l): #comment
                    continue
                else:
                    m = re.match('\s*([\w\.]+)\s*\=\s*([\w\s\-\."]+)', l)
                    if m is None:
                        print(f'Line {i+1} is not of the expected format: {l}')
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
                                    # replace d by e
                                    if 'd' in v:
                                        indx = v.find('d')
                                        v = v[:indx] + 'e' + v[indx+1:]
                                    values[j] = float(v)
                                params[m.group(1)] = values
                            except ValueError:
                                # pass them all as strings instread
                                value = values
                        if len(params[m.group(1)]) == 1:
                            params[m.group(1)] = params[m.group(1)][0]

    except:
        sys.exit("error opening the input file")

    return params
