# coding: utf-8
# calculate VSH in rotated and/or translated  frames of reference

import numpy as np
from math import sqrt
from transformation_coefficients import WignerD_matrices as WignerDs
from transformation_coefficients import translation_matrix as translations
from transformation_coefficients import translation_matrix_debug as translationsD


class HarmonicField:
    def __init__(self, coeffs, lmax):
        """ declare a harmonic field as linear combination of 
            Hansens Mlm and Nlm vector fields represented as a 
            vector with index ordering according to:
            coeffs = [M(l=1, m=-1), M(l=1, m=0), ...
                     , M(l=lmax, m=-lmax),..., M(l=lmax, m=lmax),
                     N(l=1, m=-1), N(l=1, m=0),...
                     , N(l=lmax, m=lmax)]
        """
        self.no_coeffs = (lmax + 2) * lmax * 2
        self.lmax = lmax
        assert type(coeffs) == np.ndarray
        assert np.shape(coeffs) == (self.no_coeffs, 1)
        self.coeffs = coeffs
        self.coeffs.dtype = np.complex128
        # determine the offset to jump from magnetic to electric
        self.jump = lmax*(lmax + 2)
    
    def get_coef_at(self, m):
        """ get a mask of magnetic and electric coefficients at 
            a specific value of m
        """
        valid_l = np.arange(max(abs(m), 1), self.lmax + 1)
        pos_M = np.array([l*(l+1) + m -1 for l in valid_l])
        pos_M = np.append(pos_M, pos_M+self.jump, axis=0)
        if m >=0:
            return pos_M
        return pos_M - 2*m

    def rotate(self, theta, phi):
        """ calculate the multipole coefficients after a rotation of
            theta around the y-axis followed by a rotation of 
            phi around the z-axis
        """
        Ds = WignerDs(self.lmax, theta, phi)
        # get the index, where a block of fixed l starts for
        # magnetic coefficients
        start_idx = lambda L: L**2 - 1
        lwidth = lambda L: 2*L + 1
        jump = self.jump 
        for l in range(1, self.lmax + 1):
            i0 = l**2 - 1
            i1 = i0 + 2*l + 1
            self.coeffs[i0:i1] = np.dot(Ds[l-1], self.coeffs[i0:i1])
            self.coeffs[i0+jump:i1+jump] = np.dot(Ds[l-1], 
                                self.coeffs[i0+jump:i1+jump])

    def z_translate(self, kd, sign_z, regreg, debug = False):
        """ calculate the VSH-coefficients in a frame of 
            reference which is translated kd/k along the 
            z-axis, k beeing the modulus of the wavevector
            Formulas of the translation coefficients according to:
            R. C. Wittmann, IEEE Trans. Antennas Propag. 34 (1988)
        """
        trans = lambda l1, l2, m, kd, sign_z, regreg: translations(l1,l2,m,kd,sign_z,regreg)
        if debug:
            trans = lambda l1, l2, m, kd, sign_z, regreg: translationsD(l1,l2,m,kd,sign_z,regreg)

        for m in range(1, self.lmax + 1):
            # existing values of l
            # translation matrix
            Tmatrix = trans(self.lmax, self.lmax, 
                    m, kd, sign_z, regreg)
            # update the coefficients belonging to a fixed value of m
            # and use the fact that translations are sign(m) invariant
            pos_M = self.get_coef_at(m) 
            neg_M = self.get_coef_at(-m)
            self.coeffs[pos_M] = np.dot(Tmatrix, self.coeffs[pos_M])
            self.coeffs[neg_M] = np.dot(Tmatrix, self.coeffs[neg_M])
        # now also translate m = 0
        zero_M = self.get_coef_at(0)
        Tmatrix = translations(self.lmax, self.lmax, 0, kd, sign_z, regreg)
        self.coeffs[zero_M] = np.dot(Tmatrix, self.coeffs[zero_M])

if __name__ == "__main__":
    np.set_printoptions(precision=2)
    lmax = 5 
    ncoeffs = 2*lmax*(lmax+2)
    m = np.zeros(ncoeffs, dtype = np.complex128).reshape(ncoeffs, 1)
    l1_coeffs =  np.array([-1, 0.01, 1], dtype = np.complex128)[:, np.newaxis]
    m[:l1_coeffs.size,:] = l1_coeffs
    coeffs = m
    coeffs2= m.copy()
    field = HarmonicField(coeffs, lmax)
    field2= HarmonicField(coeffs.copy(), lmax)
    kd = 0.15
    field.z_translate(kd, sign_z = 1, regreg = 1, debug = False)

    print("applying inverse translation \n")
    field.z_translate(kd, sign_z = 1, regreg = 1, debug = False)
    field2.z_translate(2*kd, sign_z = 1, regreg = 1, debug = False)
    print("Coeffs almost equal", np.allclose(field.coeffs, field2.coeffs, atol=1e-4))
    print("difference \n", field.coeffs- field2.coeffs)
