# coding: utf-8
# test the translation matrices for vector multipole
# fields
import numpy as np
from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_almost_equal
import unittest
from MultipoleTransformations import translation_matrix, translate_l1l2
from MultipoleTransformations import translation_matrix_debug
from MultipoleTransformations import full_translation_matrix

class Translation_functionality(unittest.TestCase):
    """ Run functionality tests
    """
    
    def setUp(self):
        self.l1max = 30 
        self.l2max = self.l1max 
        self.regreg = 1
    
    def test_elementwise(self):
        """ test versus the per-element evaluation
        """
        l1max, l2max = 3, 3
        kd = 2
        sign_z = 1
        regreg = 1
        for m in range(0, 5):
            lmin = max(abs(m), 1)
            Di, Dj = l1max - lmin, l2max - lmin
            dim1 = l1max - lmin + 1
            dim2 = l2max - lmin + 1
            Telementwise = np.zeros((2*dim1, 2*dim2), dtype=np.complex128)
            for l1 in range(lmin, l1max+1):
                for l2 in range(lmin, l2max+1):
                    i, j = l1-lmin, l2-lmin
                    PP, PQ = translate_l1l2(l1, l2, m, kd, sign_z, regreg)
                    Telementwise[i, j] = PP
                    Telementwise[i+dim1, j+dim2] = PP
                    Telementwise[i, j+dim2] = PQ
                    Telementwise[i+dim1, j] = PQ

            Tmatrix = translation_matrix(l1max, l2max, m, kd, sign_z, regreg)
            assert_array_almost_equal(Tmatrix, Telementwise)
            

    def test_identity(self):
        """ A Translation of a distance +d in z-direction
            followed by a translation of -d in z-direction
            should lead to a unit-matrix
        """

        d = 0.42 # distance of z-shift times wavenumber
        sign_z = 1 # translation in pos. z-direction
        regreg = 1 # regular -> regular translations
        l1max = self.l1max
        l2max = self.l2max + 5
        m = 1 
        dim = 2*(l1max - max(abs(m), 1) + 1)
        T = full_translation_matrix(l2max, l2max, d, sign_z, regreg)
        Ti= full_translation_matrix(l2max, l2max, d, -sign_z, regreg)
        dot = np.dot(T, Ti)
        assert_array_almost_equal(dot, np.eye(dot.shape[0]), decimal=3)

    def test_double_pos_translation(self):
        """ Test if a translation of distance d, twice, is the same
            as a single translation of a distance 2*d
        """
        l1max = self.l1max
        l2max = self.l1max
        d = 0.51
        regreg = self.regreg
        sign_z = 1
        T1d = full_translation_matrix(l1max, l2max, d, sign_z, regreg)
        T2d = full_translation_matrix(l2max, l1max, d, sign_z, regreg)
        Tdd = full_translation_matrix(l1max, l1max, 2*d, sign_z, regreg)
        Tdot = np.dot(T1d, T2d)
        Di = l1max 
        Dj = l2max
        assert_array_almost_equal(Tdot[0:Di-2,0:Dj-2], Tdd[0:Di-2,0:Dj-2], decimal = 2)
        assert_array_almost_equal(Tdot[Di:-2,Dj:-2], Tdd[Di:-2,Dj:-2], decimal = 2)
        assert_array_almost_equal(Tdot[0:Di-2,Dj:-2], Tdd[0:Di-2,Dj:-2], decimal = 2)
        assert_array_almost_equal(Tdot[Di:-2,0:Dj-2], Tdd[Di:-2,0:Dj-2], decimal = 2)

    def test_double_neg_translation(self):
        """ Test if a translation of distance -d, twice, is the same
            as a single translation of a distance -2*d
        """
        l1max = self.l1max
        l2max = self.l2max
        m = 2
        d = 0.51
        regreg = self.regreg
        sign_z = -1
        T1d = translation_matrix(l1max, l2max, m, d, sign_z, regreg)
        T2d = translation_matrix(l2max, l1max, m, d, sign_z, regreg)
        Tdd = translation_matrix(l1max, l1max, m, 2*d, sign_z, regreg)
        Tdot = np.dot(T1d, T2d)
        Di = l1max 
        Dj = l1max
        assert_array_almost_equal(Tdot[0:Di-2,0:Dj-2], Tdd[0:Di-2,0:Dj-2], decimal = 2)
        assert_array_almost_equal(Tdot[Di:-2,Dj:-2], Tdd[Di:-2,Dj:-2], decimal = 2)
        assert_array_almost_equal(Tdot[0:Di-2,Dj:-2], Tdd[0:Di-2,Dj:-2], decimal = 2)
        assert_array_almost_equal(Tdot[Di:-2,0:Dj-2], Tdd[Di:-2,0:Dj-2], decimal = 2)


if __name__ == "__main__":
    unittest.main()
