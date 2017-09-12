import numpy as np
from eos import eos_type_module, eos_module
from riemann_util import riemann_util_module as riemann

def amrex_probinit (init,probin,namlen,problo,probhi):# bind(c)

    """use bl_constants_module
    use probdata_module
    use prob_params_module, only : center
    use bl_error_module
    use eos_type_module, only: eos_t, eos_input, eos_input_rp
    use eos_module, only: eos
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: init, namlen
    integer, intent(in) :: name(namlen)
    real(rt), intent(in) :: problo[2], probhi[2]

    integer :: untin, i

    type(eos_t) :: eos_state

    namelist /fortin/ p_ambient, dens_ambient, exp_energy, &
       r_init, nsub, temp_ambient

    # Build "probin" filename -- the name of file containing fortin namelist.
    integer, parameter :: maxlen = 256
    character :: probin*(maxlen)

    if (namlen .gt. maxlen) then
     call bl_error('probin file name too long')
    end if"""

    # set namelist defaults

    p_ambient = 1.e-5        # ambient pressure (in erg/cc)
    dens_ambient = 1.e0      # ambient density (in g/cc)
    exp_energy = 1.e0        # absolute energy of the explosion (in erg)
    r_init = 0.05e0          # initial radius of the explosion (in cm)
    nsub = 4
    temp_ambient = -1.e2     # Set original temp. to negative, which is overwritten in the probin file

    # set local variable defaults
    center = np.zeros(3)
    center[0] = 0.5*(problo[0] + probhi[0])
    center[1] = 0.5*(problo[1] + probhi[1])
    center[2] = 0.5*(problo[2] + probhi[2])

    # Read namelists
    params = read_probin(probin)

    extract_dict(params, 'p_ambient', p_ambient)
    extract_dict(params, 'dens_ambient', dens_ambient)
    extract_dict(params, 'temp_ambient', temp_ambient)
    extract_dict(params, 'exp_energy', exp_energy)
    extract_dict(params, 'r_init', r_init)
    extract_dict(params, 'nsub', nsub)

    xn_zone[:] = 0.0
    xn_zone[0] = 1.0

    eos_state = eos_type_module.Eos_T()

    # override the pressure iwth the temperature
    if (temp_ambient > 0.0):
        eos_state.rho = dens_ambient
        eos_state.xn[:] = xn_zone[:]
        eos_state.T = temp_ambient

        eos(eos_input_rt, eos_state)

        p_ambient = eos_state.p


    # Calculate ambient state data

    eos_state.rho = dens_ambient
    eos_state.p   = p_ambient
    eos_state.T   = 1.e5 # Initial guess for iterations
    eos_state.xn  = xn_zone

    eos_module.eos(eos_type_module.eos_input_rp, eos_state)

    e_ambient = eos_state.e


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
def ca_initdata(level, time, lo, hi,
               state_dims, q_dims, delta, xlo, xhi):

    """use probdata_module
    use bl_constants_module, only: M_PI, 4.0/3.0, 0.0, 1.0
    use meth_params_module , only: NVAR, NQ, QRHO, QU, QV, QW, QREINT, QPRES, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
    use prob_params_module, only : center
    use amrex_fort_module, only : rt => amrex_real
    use network, only : nspec
    use eos_module, only : eos
    use eos_type_module, only : eos_t, eos_input_rp, eos_input_re
    use riemann_util_module, only: gr_cons_state

    implicit none

    integer :: level, nscal
    integer :: lo[2], hi[2]
    integer :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
    real(rt) :: xlo[2], xhi[2], time, delta[2]
    real(rt) :: state[state_l1:state_h1, &
                    state_l2:state_h2, &
                    state_l3:state_h3,NVAR)
    real(rt)         :: q[state_l1:state_h1, &
                    state_l2:state_h2, &
                    state_l3:state_h3,NQ)

    real(rt) :: xmin,ymin,zmin
    real(rt) :: xx, yy, zz
    real(rt) :: dist
    real(rt) :: eint, p_zone
    real(rt) :: vctr, p_exp

    integer :: i,j,k, ii, jj, kk
    integer :: npert, nambient
    real(rt) :: e_zone, gamma_up[8], gamma
    type(eos_t) :: eos_state"""

    gamma_up = np.zeros(9)
    gamma_up[0] = 1.0
    gamma_up[4] = 1.0
    gamma_up[8] = 1.0

    # set explosion pressure -- we will convert the point-explosion energy into
    # a corresponding pressure distributed throughout the perturbed volume
    vctr  = 4.0/3.0 * np.pi * r_init**3

    e_zone = exp_energy / vctr / dens_ambient

    eos_state = eos_type_module.Eos_T()

    eos_state.e = e_zone
    eos_state.rho = dens_ambient
    eos_state.xn[:] = xn_zone[:]
    eos_state.T = 100.0  # initial guess

    eos_module.eos(eos_type_module.eos_input_re, eos_state)

    p_exp = eos_state.p

    state = np.zeros((state_dims[0], state_dims[1], state_dims[2], state_dims[3]))
    q = np.zeros((q_dims[0], q_dims[1], q_dims[2], q_dims[3]))

    for k in range(lo[2], hi[2]):
        zmin = xlo[2] + delta[2]*(k-lo[2])

        for j in range(lo[1], hi[1]):
            ymin = xlo[1] + delta[1]*(j-lo[1])

            for i in range(lo[0], hi[0]):
                xmin = xlo[0] + delta[0]*(i-lo[0])

                npert = 0
                nambient = 0

                for kk in range(nsub):
                    zz = zmin + (delta[2]/nsub)*(kk + 0.5e0)

                    for jj in range(nsub):
                        yy = ymin + (delta[1]/nsub)*(jj + 0.5e0)

                        for ii in range(nsub):
                            xx = xmin + (delta[0]/nsub)*(ii + 0.5e0)

                            dist = (center[0]-xx)**2 + (center[1]-yy)**2 + (center[2]-zz)**2

                            if(dist <= r_init**2):
                                npert = npert + 1
                            else:
                                nambient = nambient + 1


                p_zone = (npert*p_exp + nambient*p_ambient) / nsub**3

                eos_state.p = p_zone
                eos_state.rho = dens_ambient
                eos_state.xn[:] = xn_zone[:]

                eos_module.eos(eos_type_module.eos_input_rp, eos_state)

                eint = dens_ambient * eos_state.e

                q[i,j,k,QRHO] = dens_ambient
                q[i,j,k,QU] = 0.e0
                q[i,j,k,QV] = 0.e0
                q[i,j,k,QW] = 0.e0

                q[i,j,k,QREINT] = eint
                q[i,j,k,QPRES] = eos_state.p

                riemann.gr_cons_state(q[i,j,k,:], state[i,j,k,:], gamma_up)

                state[i,j,k,UTEMP] = eos_state.T

                state[i,j,k,UFS] = state[i,j,k,URHO]

    return np.ascontiguousarray(q), np.ascontiguousarray(state)

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
