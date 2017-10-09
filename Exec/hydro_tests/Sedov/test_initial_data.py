from unittest.mock import patch
import unittest
from probdata import probdata_module as probdata
from prob_params import prob_params_module as prob_params
from meth_params import meth_params_module as meth_params
import numpy as np
import sys, os

# we need to mock riemann_util here as it relies on stuff from meth_params and the fortran wrapper does not set the variables in the actual fortran object properly

class TestCase(unittest.TestCase):

    def setUp(self, root_dir=os.getcwd()):
        self.probin = root_dir + "/" + "probin.3d.testa"

    def test_amrex_probinit(self):

        from Prob_3d import amrex_probinit

        problo = [0,0,0]
        probhi = [1,1,1]

        amrex_probinit(problo, probhi, self.probin)

        centre = 0.5 * (np.array(problo) + np.array(probhi))

        np.testing.assert_array_equal(prob_params.center, centre)
        np.testing.assert_equal(probdata.r_init, 1.0e-3)
        np.testing.assert_equal(probdata.p_ambient, 1.e-5)
        np.testing.assert_equal(probdata.exp_energy, 1.e-1)
        np.testing.assert_equal(probdata.dens_ambient, 1.e-1)
        np.testing.assert_equal(probdata.nsub, 10)
        np.testing.assert_equal(probdata.xn_zone[0], 1.0)
        np.testing.assert_equal(probdata.e_ambient, 1.5e-4)

    def mocked_gr_cons_state(q, state, c):
        state[:min(len(state), len(q))] = q[:min(len(state), len(q))]

    @patch('Prob_3d.riemann.gr_cons_state', new=mocked_gr_cons_state)
    def test_ca_initdata(self):

        from Prob_3d import amrex_probinit, ca_initdata

        lo = [0,0,0]
        hi = [5,5,5]
        delta = 0.1 * np.ones(3)
        xlo = [0,0,0]
        xhi = [1,1,1]

        amrex_probinit(xlo, xhi, self.probin)

        state = ca_initdata(lo, hi, lo, hi, delta, xlo, xhi, probin=self.probin[:-1])
        state = np.reshape(state, (hi[0]+1, hi[1]+1, hi[2]+1, 8), order='F')

        np.testing.assert_array_equal(state[:,:,:,0], np.ones_like(state[:,:,:,0]) * 0.1)

        np.testing.assert_array_equal(state[:,:,:,1:4], np.zeros_like(state[:,:,:,1:4]))

        np.testing.assert_array_equal(state[:,:,:,4], np.zeros_like(state[:,:,:,4]))

        np.testing.assert_array_equal(state[:,:,:,5], np.ones_like(state[:,:,:,5]) * 1.e-5)

        np.testing.assert_array_equal(state[:,:,:,6], np.ones_like(state[:,:,:,6]) * 1.e2)


if __name__ == "__main__":
    unittest.main()
