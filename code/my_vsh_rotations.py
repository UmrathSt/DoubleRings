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
        assert np.shape(coeffs) == np.shape(self.no_coeffs, 1)
        self.coeffs = coeffs

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
        # determine the offset to jump from magnetic to electric
        jump = lmax*(lmax + 2)
        for l in range(1, self.lmax + 1):
            i0 = l**2 - 1
            i1 = i0 + 2*l + 1
            self.coeffs[i0:i1] = np.dot(Ds[l-1], self.coeffs[i0:i1])
            self.coeffs[i0+jump:i1+jump] = np.dot(Ds[l-1], 
                                self.coeffs[i0+jump:i1+jump])


    def z_translate(self, kd, sign_z, regreg):
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
        [PPm, PQm] = translations(self.lmax, self.lmax, kd, sign_z, regreg)
        for m in range(-self.lmax, 1):
            PP, PQ = PPm[self.lmax+m], PQm[self.lmax+m] # square matrices for 
            mc_neg = np.array([[self.mcoeffs[l2*(l2+1)-1+m,0]] 
                      for l2 in range(max(abs(m),1), self.lmax+1)])
            nc_neg = np.array([[self.ncoeffs[l2*(l2+1)-1+m,0]] 
                      for l2 in range(max(abs(m),1), self.lmax+1)])
            mc_pos = np.array([[self.mcoeffs[l2*(l2+1)-1-m,0]] 
                      for l2 in range(max(abs(m),1), self.lmax+1)])
            nc_pos = np.array([[self.ncoeffs[l2*(l2+1)-1-m,0]] 
                      for l2 in range(max(abs(m),1), self.lmax+1)])
            print("m = %i" %m)
            l_min = max(abs(m), 1)
            for l in range(l_min, self.lmax+1):  
                i = l - l_min
                j = l*(l+1)+m-1
                print("summing line for (l,m)=(%i,%i)" %(l,m))
                print("mc_neg", mc_neg)
                mcoeffs[j,0] = (PP[i,:] * mc_neg + PQ[i,:] * nc_neg).sum()
                ncoeffs[j,0] = (PP[i,:] * nc_neg + PQ[i,:] * mc_neg).sum()
                if not m == 0:
                    mcoeffs[j+abs(2*m),0] = (PP[i,:] * mc_pos + 
                                             PQ[i,:] * nc_pos).sum()
                    ncoeffs[j+abs(2*m),0] = (PP[i,:] * nc_pos + 
                                             PQ[i,:] * mc_pos).sum()
        self.mcoeffs = mcoeffs
        self.ncoeffs = ncoeffs
    
        

if __name__ == "__main__":
    np.set_printoptions(precision=2)
    m = np.array([-1, 0.01, 1, -20, -10, 0.1, 10, 20], dtype = np.complex128)[:, np.newaxis]
    n = np.array([0, 0, 0, 0, 0, 0, 0, 0], dtype = np.complex128)[:, np.newaxis]
    field = HarmonicField(m, n, 2)
    kd = 0.05
    print("applying translation to mcoeffs = \n", m)
    field.z_translate(kd, 1, 1)
    print("applying inverse translation \n")
    field.z_translate(kd, -1, 1)
    print("M-Coeffs")
    print(field.mcoeffs)
