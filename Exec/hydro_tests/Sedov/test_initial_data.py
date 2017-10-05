from unittest.mock import patch, Mock
import unittest
from probdata import probdata_module as probdata
from prob_params import prob_params_module as prob_params
from meth_params import meth_params_module as meth_params
import numpy as np
import sys

# we need to mock riemann_util here as it relies on stuff from meth_params and the fortran wrapper does not set the variables in the actual fortran object properly

class mocked_riemann:
    class riemann_util_module:
        def gr_cons_state(q, state, c):
            state[:] = q[:len(state)]

sys.modules['riemann_util'] = mocked_riemann

import Prob_3d

class TestCase(unittest.TestCase):

    def test_amrex_probinit(self):

        probin = "probin.3d.testa"

        problo = [0,0,0]
        probhi = [1,1,1]

        Prob_3d.amrex_probinit(problo, probhi, probin)

        print(f'probdata = {probdata}')

        centre = 0.5 * (np.array(problo) + np.array(probhi))

        np.testing.assert_array_equal(prob_params.center, centre)
        np.testing.assert_equal(probdata.r_init, 1.0e-3)
        np.testing.assert_equal(probdata.p_ambient, 1.e-5)
        np.testing.assert_equal(probdata.exp_energy, 1.e-1)
        np.testing.assert_equal(probdata.dens_ambient, 1.e-1)
        np.testing.assert_equal(probdata.nsub, 10)
        np.testing.assert_equal(probdata.xn_zone[0], 1.0)
        np.testing.assert_equal(probdata.e_ambient, 1.5e-4)


    def test_ca_initdata(self):

        lo = [0,0,0]
        hi = [5,5,5]
        delta = 0.1 * np.ones(3)
        xlo = [0,0,0]
        xhi = [1,1,1]

        probin = "probin.3d.testa"

        Prob_3d.amrex_probinit(xlo, xhi, probin)

        state = Prob_3d.ca_initdata(lo, hi, lo, hi, delta, xlo, xhi, probin=probin[:-1])
        state = np.reshape(state, (hi[0]+1, hi[1]+1, hi[2]+1, 7), order='F')

        np.testing.assert_array_equal(state[:,:,:,0], np.ones_like(state[:,:,:,0]) * 0.1)

        np.testing.assert_array_equal(state[:,:,:,1:4], np.zeros_like(state[:,:,:,1:4]))

        np.testing.assert_array_equal(state[:,:,:,4], np.zeros_like(state[:,:,:,4]))

        np.testing.assert_array_equal(state[:,:,:,5], np.ones_like(state[:,:,:,5]) * 1.e-5)

        np.testing.assert_array_equal(state[:,:,:,6], np.ones_like(state[:,:,:,6]) * 1.e2)


if __name__ == "__main__":
    unittest.main()
