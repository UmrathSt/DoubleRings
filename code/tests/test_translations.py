# coding: utf-8
# test the translation matrices for vector multipole
# fields
import numpy as np
from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_almost_equal
import unittest
from MultipoleTransformations import translation_matrix, translate_l1l2
from MultipoleTransformations import translation_matrix_debug


class Translation_functionality(unittest.TestCase):
    """ Run functionality tests
    """
    
    def setUp(self):
        self.l1max = 20
        self.l2max = 20
        self.regreg = 1
    
    def test_identity(self):
        """ A Translation of a distance +d in z-direction
            followed by a translation of -d in z-direction
            should lead to a unit-matrix
        """

        d = 0.42 # distance of z-shift times wavenumber
        sign_z = 1 # translation in pos. z-direction
        regreg = 1 # regular -> regular translations
        l1max = self.l1max
        l2max = self.l2max
        m = 1 
        dim = 2*(l1max - max(abs(m), 1) + 1)
        T = translation_matrix(l1max, l2max, m, d, sign_z, regreg)
        Ti= translation_matrix(l2max, l1max, m, d, -sign_z, regreg)
        assert_array_almost_equal(np.dot(T, Ti), np.eye(dim))

    def test_double_translation(self):
        """ Test if a translation of distance d, twice, is the same
            as a single translation of a distance 2*d
        """
        l1max = self.l1max
        l2max = self.l2max
        m = 2
        d = 0.51
        regreg = self.regreg
        sign_z = 1
        T1d = translation_matrix(l1max, l2max, m, d, sign_z, regreg)
        T2d = translation_matrix(l2max, l1max, m, d, sign_z, regreg)
        Tdd = translation_matrix(l1max, l1max, m, 2*d, sign_z, regreg)
        assert_array_almost_equal(np.dot(T1d, T2d), Tdd)
        

if __name__ == "__main__":
    unittest.main()
