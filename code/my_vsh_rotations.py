# coding: utf-8
# calculate VSH in rotated and/or translated  frames of reference

import numpy as np
from math import sqrt
from transformation_coefficients import WignerD_matrices as WignerDs, translation_matrix as translations



class HarmonicField:
    def __init__(self, mcoeffs, ncoeffs, lmax):
        """ declare a harmonic field as linear combination of 
            Hansens Mlm and Nlm vector fields represented as a 
            vector with index ordering according to:
            mcoeffs = [M(l=1, m=-1), M(l=1, m=0), ...
                     , M(l=lmax, m=-lmax),..., M(l=lmax, m=lmax)]
        """
        assert type(mcoeffs) == np.ndarray
        assert np.shape(mcoeffs) == np.shape(ncoeffs)
        self.no_coeffs = (lmax + 2) * lmax
        self.lmax = lmax
        assert mcoeffs.shape == (self.no_coeffs, 1)
        self.mcoeffs = mcoeffs
        self.ncoeffs = ncoeffs

    def rotate(self, theta, phi):
        """ calculate the m- and ncoeffs after a rotation of
            theta around the y-axis followed by a rotation of 
            phi around the z-axis
        """
        Ds = WignerDs(self.lmax, theta, phi)
        mcoeffs = np.zeros(self.mcoeffs.shape, dtype=np.complex128)
        ncoeffs = np.zeros(self.ncoeffs.shape, dtype=np.complex128)
        start_idx = lambda L, M: (M)*(2*L + 1) + L**2 -1 
        for l in range(1, self.lmax + 1):
            i0 = int(4*l**3 - l - 3)
            i1 = i0 + (2*l + 1)**2
            mcoeffs[i0:i1] = np.dot(Ds[l-1], self.mcoeffs[i0:i1])
            ncoeffs[i0:i1] = np.dot(Ds[l-1], self.ncoeffs[i0:i1])
        self.mcoeffs = mcoeffs
        self.ncoeffs = ncoeffs


    def z_translate(self, kd, sign_z):
        """ calculate the VSH-coefficients in a frame of 
            reference which is translated kd/k along the 
            z-axis, k beeing the wavevector
            Translation coefficients according to:
            R. C. Wittmann, IEEE Trans. Antennas Propag. 34 (1988)
        """
        l1 = np.arange(1, self.lmax+1).reshape(self.lmax, 1)
        l2 = l1.transpose()
        mcoeffs = self.mcoeffs.copy()
        ncoeffs = self.ncoeffs.copy()
        for l in range(1, self.lmax+1):
            for m in range(-self.lmax, 1):
                for l2 in range(1, self.lmax+1):
                    PP, PQ = self.Vl1l2_m(l, l2, m, kd, sign_z)
                    idx = l2*(l2+1) - 1 + m
                    mcoeffs[l*(l+1)-1+m] += (self.mcoeffs[idx] * PP +
                                             self.ncoeffs[idx] * PQ)
                    ncoeffs[l*(l+1)-1+m] += (self.ncoeffs[idx] * PP +
                                             self.mcoeffs[idx] * PQ)
                    mcoeffs[l*(l+1)-1-m] += (self.mcoeffs[idx-2*m] * PP +
                                             self.ncoeffs[idx-2*m] * PQ)
                    ncoeffs[l*(l+1)-1-m] += (self.ncoeffs[idx-2*m] * PP +
                                             self.mcoeffs[idx-2*m] * PQ)
        self.mcoeffs = mcoeffs
        self.ncoeffs = ncoeffs

                                        

        

if __name__ == "__main__":
    m = np.array([1j, 0, 1j])[:, np.newaxis]
    n = np.array([1, 0, -1])[:, np.newaxis]
    field = HarmonicField(m, n, 1)
    field.rotate(np.pi/2, 0)
    print("applying rotation to mcoeffs = ", m, "ncoeffs = ", n)
    print("M-Coeffs")
    print(field.mcoeffs)
    print("N-Coeffs")
    print(field.ncoeffs)
    print("rotating back")
    field.rotate(-np.pi/2,0)
    print("M-Coeffs")
    print(field.mcoeffs)
    print("N-Coeffs")
    print(field.ncoeffs)
