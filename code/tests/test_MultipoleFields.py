# test some functionality of MultipoleFields.py
# coding: utf-8

import unittest
import numpy as np
from numpy.testing import assert_array_almost_equal
from MultipoleFields import HarmonicField

class MultipoleFields_coefficients(unittest.TestCase):

    def test_get_m_blocks(self):

        lmax = 2
        coeffs = np.array([[i] for i in range(2*lmax*(lmax+2))])
        m0_coeffs = np.array([[1], [5], [9], [13]])
        mM1_coeffs= np.array([[0], [4], [8], [12]])
        mP1_coeffs= np.array([[2], [6], [10], [14]])
        mM2_coeffs = np.array([[3], [11]])
        mP2_coeffs = np.array([[7], [15]])
        m_list = [mM2_coeffs, mM1_coeffs, m0_coeffs,
                  mP1_coeffs, mP2_coeffs]

        H = HarmonicField(coeffs, lmax)
        for m in range(-lmax, lmax+1):
            mBlock = (H.coeffs[H.get_coef_at(m)]).real
            assert_array_almost_equal(mBlock, m_list[m+lmax])

    def test_lmax_reduction(self):

        lmax = 5
        n0 = lmax*(lmax + 2)
        coeffs = np.arange(2*n0)[:, np.newaxis]
        hf = HarmonicField(coeffs, lmax)
        
        l1_coeffs = np.append(coeffs[0:3,:], coeffs[n0:n0+3,:], axis=0)
        hf.reduce_lmax(1)
        assert_array_almost_equal(l1_coeffs, hf.coeffs)

    def test_rotation_translation(self):
        """ test if a multipole field translated by d in positive
            z-direction is the same as a field rotated by pi around the
            y-axis and then translated a distance -d in the z-direction
        """
        lmax = 15
        lmax_red = 5
        kd = 1 
        sign_z = 1
        regreg = 1
        coeffs = np.arange(lmax*(lmax + 2)*2)[:,np.newaxis]
        coeffs1 = coeffs.copy()
        coeffs2 = coeffs.copy()
        hf1 = HarmonicField(coeffs1, lmax)
        hf1.z_translate_matrix(kd, sign_z, regreg)
        hf1.reduce_lmax(lmax_red)
        hf2 = HarmonicField(coeffs2, lmax)
        hf2.rotate_matrix(np.pi, 0)
        hf2.z_translate_matrix(kd, -sign_z, regreg)
        hf2.rotate_matrix(-np.pi, 0)
        hf2.reduce_lmax(lmax_red)
        assert_array_almost_equal(hf1.coeffs, hf2.coeffs, decimal=4)
if __name__ == "__main__":
    unittest.main()

