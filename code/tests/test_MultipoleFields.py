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
        coeffs = np.array([[i] for i in range(2*lmax*(lmax+2))])
        hf = HarmonicField(coeffs, lmax)
        
        l1_coeffs = np.append(coeffs[0:3,:], coeffs[n0:n0+3,:], axis=0)
        hf.reduce_lmax(1)
        assert_array_almost_equal(l1_coeffs, hf.coeffs)

if __name__ == "__main__":
    unittest.main()

