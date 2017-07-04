# coding: utf-8
# calculate VSH in rotated and/or translated  frames of reference

import numpy as np
from math import sqrt
import quaternion
import spherical_functions as sf
from scipy.special import jv, hankel1
from cython_gaunt import gaunt

jn = lambda l, x: jv(l + 0.5, x)*np.sqrt(np.pi/(2*x))
h1 = lambda l, x: hankel1(l + 0.5, x)*np.sqrt(np.pi/(2*x))


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
        mcoeffs = np.zeros(self.mcoeffs.shape, dtype=np.complex128)
        ncoeffs = np.zeros(self.ncoeffs.shape, dtype=np.complex128)
        start_idx = lambda L, M: (M)*(2*L + 1) + L**2 -1 
        for l in range(1, self.lmax + 1):
            for j in range(0, 2*l + 1):
                m = j - l 
                m0 = l**2 - 1 
                m1 = m0 + 2*l + 1
                D0 = (2*l + 1)*(m + l) + l**2 - 1
                D1 = D0 + 2*l +1
                Dlm = D[D0:D1]
                mcoeffs[l**2  - 1 + j] = (
                        self.mcoeffs[m0:m1] * Dlm).sum()
                ncoeffs[l**2  - 1 + j] = (
                        self.ncoeffs[m0:m1] * Dlm).sum()
        self.mcoeffs = mcoeffs
        self.ncoeffs = ncoeffs
    def Vl1l2_m(self, l1, l2, m, kd, sign_z):
        """return the (l1, m) contribution of a vector spherical 
           harmonic M_l1,m (N_l1,m) when it is translated along (sign_z = +1)
           or against (sign_z = -1) the z-direction an adimensional distance
           kd. 
           Returns the polarization conserving contribution M_l1,m (N_l1,m) 
           as the first and the polarization mixing contribution N_l1,m (M_l1,m)
           as the second value
        """
        
        common =  ((-1)**m*4*np.pi * 1j**(l2-l1)/
                   sqrt(l2*(l2+1) * l1*(l1+1)))
        alpha = np.arange(abs(l1-l2), l1+l2+1, 2)
        gaunts = gaunt(l1, l2, m)
        h1_sph= h1(alpha+0.5, kd)
        fak_PP = ((l1*(l1+1) + l2*(l2+1) - alpha*(alpha+1))*0.5 * h1_sph *
                    (sign_z*1j)**alpha*gaunts*np.sqrt((2*alpha+1)/(4*np.pi))).sum()
        fak_PQ = (-m*kd*(1j)**(alpha+1) * gaunts * np.sqrt((2*alpha+1)/(4*np.pi))
                  *h1_sph * sign_z**(alpha+1)).sum()
        PP = fak_PP * common
        PQ = fak_PQ * common
        return PP, PQ

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
    m = np.array([1j, 0, 1j])
    n = np.array([1, 0, -1])
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
