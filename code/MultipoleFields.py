# coding: utf-8
# calculate VSH in rotated and/or translated  frames of reference

import numpy as np
from math import sqrt
from transformation_coefficients import WignerD_matrices as WignerDs
from transformation_coefficients import translation_matrix as translations



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

    def z_translate(self, kd, sign_z, regreg):
        """ calculate the VSH-coefficients in a frame of 
            reference which is translated kd/k along the 
            z-axis, k beeing the modulus of the wavevector
            Formulas of the translation coefficients according to:
            R. C. Wittmann, IEEE Trans. Antennas Propag. 34 (1988)
        """
        for m in range(1, self.lmax + 1):
            # existing values of l
            valid_l = np.arange(max(abs(m), 1), self.lmax + 1)
            # translation matrix
            Tmatrix = translations(self.lmax, self.lmax, 
                    m, kd, sign_z, regreg)
            # update the coefficients belonging to a fixed value of m
            # and use the fact that translations sign(m) invariant
            pos_M = np.array([l*(l+1) + m -1 for l in valid_l])
            neg_M = pos_M - 2*m
            pos_m_coeffs = np.append(pos_M, pos_M + self.jump, axis = 0)
            neg_m_coeffs = np.append(neg_M, neg_M + self.jump, axis = 0)
            print("m>0 block for m = %i: \n" %m, self.coeffs[pos_m_coeffs])
            print("m<0 block for |m| = %i: \n" %m, self.coeffs[neg_m_coeffs])
            self.coeffs[pos_m_coeffs] = np.dot(Tmatrix, 
                                        self.coeffs[pos_m_coeffs])
            
            self.coeffs[neg_m_coeffs] = np.dot(Tmatrix,
                                        self.coeffs[neg_m_coeffs])
        # now also translate m = 0
        zero_M = np.array([l*(l+1) - 1 for l in range(1, self.lmax+1)])
        zero_m_coeffs = np.append(zero_M, zero_M + self.jump, axis = 0)
        Tmatrix = translations(self.lmax, self.lmax, 0, kd, sign_z, regreg)
        self.coeffs[zero_m_coeffs] = np.dot(Tmatrix, 
                                    self.coeffs[zero_m_coeffs])

if __name__ == "__main__":
    np.set_printoptions(precision=2)
    m = np.array([-1, 0.01, 1, -20, -10, 0.1, 10, 20], dtype = np.complex128)[:, np.newaxis]
    n = np.array([0, 0, 0, 0, 0, 0, 0, 0], dtype = np.complex128)[:, np.newaxis]
    coeffs = np.append(m, n, axis=0)
    field = HarmonicField(coeffs, 2)
    kd = 0.05
    print("applying translation to Coeffs = \n", coeffs)
    field.z_translate(kd, sign_z = 1, regreg = 1)
    print("applying inverse translation \n")
    field.z_translate(kd, sign_z = -1, regreg = 1)
    print("Coeffs")
    print(field.coeffs)
