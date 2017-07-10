# coding: utf-8
# calculate VSH in rotated and/or translated  frames of reference

import numpy as np
from math import sqrt
from MultipoleTransformations import WignerD_matrices as WignerDs
from MultipoleTransformations import translation_matrix as translations
from MultipoleTransformations  import translation_matrix_debug as translationsD


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
        if not self.coeffs.dtype == np.complex128:
            self.coeffs = self.coeffs + 0j 
        # determine the offset to jump from magnetic to electric
        # self.jump will be the line-index of first electric 
        # coefficient QE(l=1, m=-1)
        self.jump = lmax*(lmax + 2)
    
    def reduce_lmax(self, new_lmax):
        assert self.lmax >= new_lmax
        if self.lmax == new_lmax:
            pass
        else:
            new_max_idx = new_lmax*(new_lmax + 2)
            E0_idx = self.lmax*(self.lmax + 2)
            Qmagnetic = self.coeffs[0:new_max_idx, :]
            new_max_idx += E0_idx
            Qelectric = self.coeffs[E0_idx:new_max_idx,:]
            self.coeffs = np.append(Qmagnetic, Qelectric, axis=0)
    
    def get_coef_at(self, m):
        """ get a mask of magnetic and electric coefficients at 
            a specific value of m
        """
        valid_l = np.arange(max(abs(m), 1), self.lmax + 1)
        pos_M = np.array([l*(l+1) + m -1 for l in valid_l])
        pos_M = np.append(pos_M, pos_M + self.jump, axis=0)
        return pos_M

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
            pos_m_idx = self.get_coef_at(m) 
            neg_m_idx = self.get_coef_at(-m)
            self.coeffs[pos_m_idx, 0] = np.dot(Tmatrix, self.coeffs[pos_m_idx])[0,0]
            self.coeffs[neg_m_idx, 0] = np.dot(Tmatrix, self.coeffs[neg_m_idx])[0,0]
        # now also translate m = 0
        m0_idx = self.get_coef_at(0)
        Tmatrix = translations(self.lmax, self.lmax, 0, kd, sign_z, regreg)
        self.coeffs[m0_idx, 0] = np.dot(Tmatrix, self.coeffs[m0_idx])[0,0]

if __name__ == "__main__":
    np.set_printoptions(precision=3)
    lmax = 2 
    ncoeffs = 2*lmax*(lmax+2)
    m = np.zeros(ncoeffs, dtype = np.complex128).reshape(ncoeffs, 1)
    l1_coeffs =  np.array([-1, 0.01, 1, -20, -10, 0.1, 10, 20], dtype = np.complex128)
    m[:l1_coeffs.size,0] = l1_coeffs
    coeffs = m.copy()
    field = HarmonicField(coeffs, lmax)
    print("Koeffizienten")
    print(field.coeffs)
    kd = 5
    field.z_translate(kd, sign_z = 1, regreg = 1)
    field.z_translate(kd, sign_z =-1, regreg = 1)
    field.reduce_lmax(1)
    #print("Coeffs almost equal", np.allclose(field.coeffs, m, atol=1e-4))
    print("Dipole-coeffs \n", field.coeffs)
