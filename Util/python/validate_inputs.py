import sys
import re
import numpy as np

def check_parameter_file(filename, dim=3):
    try:
        f = open(filename)
    except:
        sys.exit("error opening the input file")

    lines = f.readlines()
    params = {}

    # first check formatting and create dictionary
    for i, l in enumerate(lines):
        if l[0] == '#' or re.match('\s+', l): #comment
            continue
        else:
            m = re.match('([\w\.]+)\s*\=\s*([\w\s\-\."]+)', l)
            if m is None:
                # check to see if they contain some input data
                m1 = re.match('(\d+\.?\d*E?\-?\d*)', l)
                m2 = re.match('(\d+\.?\d*E?\-?\d*)\s(\d+\.?\d*E?\-?\d*)', l)
                if not m1 and not m2:
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

    # check that all the essential parameters are in there

    essential_params = [['max_step', 'stop_time'], 'geometry.is_periodic', 'geometry.prob_lo', 'geometry.prob_hi',
                        'amr.n_cell', 'amr.max_level', 'castro.lo_bc', 'castro.hi_bc']

    for p in essential_params:
        if type(p) == list:
            found_option = False
            for i in p:
                if i in params:
                    found_option = True
                    break
            if not found_option:
                sys.exit("Input file did not define one of the essential parameters: {}".format(p))

        else:
            if p not in params:
                sys.exit("Input file did not define one of the essential parameters: {}".format(p))

    validate_params(params, dim)

    return params

def validate_params(params, dim):
    # validate parameters

    def check_real_array(param, min_length=dim):
        if param in params:
            if min_length == 1 and np.isscalar(params[param]):
                ps = [params[param]] # make list so can loop over
            else:
                ps = params[param]

            assert len(ps) >= min_length
            assert np.isreal(ps).all()
            assert np.isfinite(ps).all()

    def check_int_array(param, minimum=None, maximum=None, min_length=dim):
        if param in params:
            if min_length == 1 and np.isscalar(params[param]):
                ps = [params[param]] # make list so can loop over
            else:
                ps = params[param]

            assert len(ps) >= min_length
            for i in range(min_length):
                assert type(ps[i]) == int
                if minimum:
                    assert ps[i] >= minimum
                if maximum:
                    assert ps[i] <= maximum

    def check_bool_array(param, min_length=dim):
        check_int_array(param, 0, 1, min_length)

    def check_string_array(param, min_length=dim):
        if param in params:
            assert len(params[param]) >= min_length

            for p in params[param]:
                assert type(p) == str

    def check_int(param, minimum=None, maximum=None):
        if param in params:
            assert np.isscalar(params[param])

            assert type(params[param]) == int
            if minimum:
                assert params[param] >= minimum
            if maximum:
                assert params[param] <= maximum

    def check_real(param, minimum=None, maximum=None):
        if param in params:
            assert np.isscalar(params[param])
            assert np.isreal(params[param])
            if minimum:
                assert params[param] >= minimum
            if maximum:
                assert params[param] <= maximum
            assert np.isfinite(params[param])

    def check_bool(param):
        check_int(param, 0, 1)

    def check_string(param):
        if param in params:
            assert np.isscalar(params[param])
            assert type(params[param]) == str

    check_real('stop_time', 0., 1.e15)
    check_int('max_step', 0, 1e10)
    check_bool('castro.use_stopping_criterion')
    check_bool('castro.use_retry')

    check_bool_array('geometry.is_periodic')
    check_real_array('geometry.prob_lo')
    check_real_array('geometry.prob_hi')
    check_int('geometry.coord_sys', 0, 0)

    check_real_array('castro.center')

    check_int_array('amr.n_cell', 1, 1e10)
    check_int('amr.max_level', 0, 20)

    check_int_array('castro.lo_bc', 0, 5)
    check_int_array('castro.hi_bc', 0, 5)

    check_bool('castro.do_hydro')
    check_int('castro.ppm_type', 0, 3)
    check_bool('castro.allow_negative_energy')
    check_real('castro.cfl', 0., 1.)
    check_real('castro.init_shrink', 0., 1.)
    check_real('castro.change_max', 0.)
    check_real('castro.dt_cutoff', 0.)

    check_real('castro.small_dens', 0.)
    check_real('castro.small_temp', 0.)

    check_int('castro.sum_interval')
    check_bool('castro.v')
    check_bool('amr.v')
    check_string('amr.grid_log')
    check_string_array('amr.data_log', 1)

    if 'amr.max_level' in params and params['amr.max_level'] > 0:
        refinement_levels = params['amr.max_level']

        check_int_array('amr.ref_ratio', 1, 10, min_length=refinement_levels)
        check_int_array('amr.regrid_int', -1, 1e5, min_length=1)
        check_int('amr.block_factor', 1)
        check_int_array('amr.n_error_buf', 0, min_length=refinement_levels)
        check_bool('castro.use_post_step_regrid')
        check_int_array('amr.max_grid_size', 1, 10000, min_length=1)

    check_bool('castro.track_grid_losses')
    check_bool('amr.plotfile_on_restart')
    check_bool('amr.checkpoint_on_restart')

    check_bool('amr.checkpoint_files_output')
    check_string('amr.check_file')
    check_real('amr.check_per', 0)
    check_int('amr.check_int')

    check_bool('amr.plot_files_output')
    check_string('amr.plot_file')
    check_real('amr.plot_per', 0)
    check_int('amr.plot_int')

    check_string_array('amr.plot_vars', min_length=1)
    check_string_array('amr.derive_plot_vars', min_length=1)

    check_string('amr.probin_file')

    # check that dx=dy=dz
    if dim > 1:
        dxs = (np.array(params['geometry.prob_hi'])[:dim] - np.array(params['geometry.prob_lo'])[:dim]) / \
                np.array(params['amr.n_cell'])[:dim]
        assert np.isclose(dxs.min(), dxs.max(), rtol=1.e-9)
    
