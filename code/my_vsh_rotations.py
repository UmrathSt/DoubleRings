# coding: utf-8
# calculate VSH in rotated and/or translated  frames of reference

import numpy as np
from math import sqrt
import quaternion
import spherical_functions as sf
from scipy.special import jv, hankel1

jn = lambda l, x: jv(l + 0.5, x)*np.sqrt(np.pi/(2*x))
h1 = lambda l, x. hankel1(l + 0.5, x)*np.sqrt(np.pi/(2*x))


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
        assert len(mcoeffs) == self.no_coeffs
        self.mcoeffs = mcoeffs
        self.ncoeffs = ncoeffs

    def rotate(self, theta, phi):
        """ calculate the m- and ncoeffs after a rotation of
            theta around the y-axis followed by a rotation of 
            phi around the z-axis
        """
        rotor = quaternion.from_spherical_coords(theta, phi)
        D = sf.WignerD.Wigner_D_matrices(rotor, 1, self.lmax)
        mcoeffs = np.zeros(self.mcoeffs.shape)
        ncoeffs = np.zeros(self.ncoeffs.shape)
        start_idx = lambda L, M: (L + M)*(2*L + 1) 
        for l in range(self.lmax):
            for m in range(2*l + 1):
                istart = start_idx(l, m)
                istop = istart + 2*l+1
                self.mcoeffs[l**2 + l - 1 + m] = (
                        self.mcoeffs[istart:istop] 
                         * D[istart:istop]).sum()
                self.ncoeffs[l**2 + l - 1 + m] = (
                        self.ncoeffs[istart:istop] 
                         * D[istart:istop]).sum()
    def z_translate(self, kd):
        """ calculate the VSH-coefficients in a frame of 
            reference which is translated kd/k along the 
            z-axis, k beeing the wavevector
            Translation coefficients according to:
            R. C. Wittmann, IEEE Trans. Antennas Propag. 34 (1988)
        """
        l = np.arange(1, self.lmax).reshape(self.lmax, 1)
        nu = l.transpose()
        common = lambda nu, l, m: ((-1)**m*4*np.pi*1j**(l-nu)/
                                    sqrt(l*(l+1)*nu*(nu+1))
                                  )
        for 
        

if __name__ == "__main__":
    field = HarmonicField(np.array([1,1,1], dtype = np.complex128), 
            np.array([0,0,0], dtype = np.complex128), 1)
    field.rotate(np.pi, 0)
    print(field.mcoeffs)
