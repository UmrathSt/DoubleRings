# coding: utf-8
# calculate VSH in rotated and/or translated  frames of reference

import numpy as np
from math import pi
from math import sqrt
from MultipoleTransformations import WignerD_matrices as WignerDs
from MultipoleTransformations import translation_matrix as translations
from MultipoleTransformations  import translation_matrix_debug as translationsD
from MultipoleTransformations import full_translation_matrix as Tfull
from MultipoleTransformations import full_rotation_matrix as Rfull
from plot_MultipoleField import plot_multipole_field as pmf


class HarmonicField:
    def __init__(self, coeffs, freq, lmax):
        """ declare a harmonic field as linear combination of 
            Hansens Mlm and Nlm vector fields represented as a 
            vector with index ordering according to:
            coeffs = [M(l=1, m=-1), M(l=1, m=0), ...
                     , M(l=lmax, m=-lmax),..., M(l=lmax, m=lmax),
                     N(l=1, m=-1), N(l=1, m=0),...
                     , N(l=lmax, m=lmax)]
        """
        if not freq >= 0:
            raise ValueError("frequency must be a positive number")
        self.freq = freq 
        self.k = 2*np.pi * freq / 299792458 # speed of light in vacuum
        self.no_coeffs = (lmax + 2) * lmax * 2
        if not type(lmax) == int and lmax >= 1:
            raise TypeError("maximum angular momentum is expected to be a natural number")
        self.lmax = lmax
        assert type(coeffs) == np.ndarray
        assert np.shape(coeffs) == (self.no_coeffs,)
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
    def rotate_matrix(self, theta, phi):
        """ applies a rotation by multiplication of the 
            multipole-coefficients with the full rotation 
            matrix
        """
        R = Rfull(self.lmax, theta, phi)
        self.coeffs = np.dot(R, self.coeffs)
 
    def z_translation_matrix(self, d, regreg, debug = False):
        """ do a translation and use matrixmultiplication by 
            the full translation matrix and decide the sign whether
            the translation is in + or - z-direction upon the sign
            of d
        """
        kd = self.k * abs(d)
        T = Tfull(self.lmax, self.lmax, kd, np.sign(d), regreg)
        self.coeffs = np.dot(T, self.coeffs)
           
    
    def z_translate_matrix(self, kd, sign_z, regreg, debug = False):
        """ do a translation and use matrixmultiplication by 
            the full translation matrix
        """
        T = Tfull(self.lmax, self.lmax, kd, sign_z, regreg)
        self.coeffs = np.dot(T, self.coeffs)
    
    def z_translation(self, d, regreg, debug = False):
        """ calculate the VSH-coefficients in a frame of 
            reference which is translated kd/k along the 
            z-axis, k beeing the modulus of the wavevector
            Formulas of the translation coefficients according to:
            R. C. Wittmann, IEEE Trans. Antennas Propag. 34 (1988)
        """
        trans = lambda l1, l2, m, kd, sign_z, regreg: translation(l1,l2,m,kd,regreg)
        if debug:
            trans = lambda l1, l2, m, kd, sign_z, regreg: translationD(l1,l2,m,kd,regreg)
        kd = self.k * d
        for m in range(1, self.lmax + 1):
            # existing values of l
            # translation matrix
            Tmatrix = trans(self.lmax, self.lmax, 
                    m, kd, regreg)
            # update the coefficients belonging to a fixed value of m
            # and use the fact that translations are sign(m) invariant
            pos_m_idx = self.get_coef_at(m) 
            neg_m_idx = self.get_coef_at(-m)
            self.coeffs[pos_m_idx, 0] = np.dot(Tmatrix, self.coeffs[pos_m_idx])[0,0]
            self.coeffs[neg_m_idx, 0] = np.dot(Tmatrix, self.coeffs[neg_m_idx])[0,0]
        # now also translate m = 0
        m0_idx = self.get_coef_at(0)
        Tmatrix = translations(self.lmax, self.lmax, 0, kd, regreg)
        self.coeffs[m0_idx, 0] = np.dot(Tmatrix, self.coeffs[m0_idx])[0,0]


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
    
    def get_field_at(self,x, y, z, reg):
        """ Get the X, Y and Z components of the given multipole 
            field at the coordinates specified as numpy arrays x, 
            y and z for regular (reg = 1) or outgoing fields (reg = 0)
        """
        field = 0
        idx_offset = self.lmax*(self.lmax+2)
        for l in range(1, self.lmax+1):
            for m in range(-l, l+1):
                Midx = l*(l+1) - 1 + m 
                Nidx = Midx + idx_offset
                Mcoeff = self.coeffs[Midx]
                Ncoeff = self.coeffs[Nidx]
                F = pmf(l, m, self.k, x, y, z, reg)
                if abs(Mcoeff) >= 1e-14:
                    field += Mcoeff * np.array(F.get_Mlm())
                if abs(Ncoeff) >= 1e-14:
                    field += Ncoeff * np.array(F.get_Nlm())
        return field


def get_z_multipole_coeffs(polarization, lmax):
    """ Get the multiple coefficients of a z-propagating
        plane wave with maximum considered angular momentum 
        "lmax" and electric field vector in "x" 
        (polarization=0) or "y"-direction (polarization=1)
    """
    if polarization not in [0, 1]:
        raise ValueError("polarization must bei either 0 or 1, given: ", polarization)
    if polarization == 0:
        Mfact = lambda m: 1
        Nfact = lambda m: m
    if polarization == 1:
        Mfact = lambda m: -1j*m
        Nfact = lambda m: 1j*m
    coeffs = np.zeros(2*lmax*(lmax+2), dtype=np.complex128)
    offset = lmax*(lmax + 2)
    for l in range(1, lmax+1):
        prefactor = 1j**(l+1) * sqrt((2*l+1)*pi)
        for m in [-1, 1]:
            idx = l*(l+1) -1 + m
            coeffs[idx] = prefactor * Mfact(m)
            coeffs[idx + offset] = prefactor * Nfact(m)
    return coeffs

def PlaneMultipoleWave(k_theta, k_phi, polarization, freq, lmax):
    """ return a HarmonicField object which represents
        a plane wave propagating in the direction of (k_theta, k_phi)
        polarized in x/y direction (polarization = 0/1) with frequency
        freq and maximum considered angular momentum lmax
    """
    F = HarmonicField(get_z_multipole_coeffs(polarization, lmax), freq, lmax)
    F.rotate(k_theta, k_phi)
    return F

if __name__ == "__main__":
    np.set_printoptions(precision=3, suppress=True)
    from matplotlib import pyplot as plt
    PW = PlaneMultipoleWave(k_theta = np.pi/32, k_phi = 0, 
                            polarization = 0, freq = 0.5e9, lmax = 15)
    x = np.linspace(-1, 1, 100)
    z = np.linspace(-2, 2, 100)
    y = np.array([0])
    field = PW.get_field_at(x, y, z, 1)
    X, Z = np.meshgrid(x, z)
    for f in [field, X, Z]:
        print("shape is ", f.shape)
    plt.pcolor(X, Z, np.sqrt(np.real(field[0,:,0,:])**2 +
                             np.real(field[1,:,0,:])**2 +
                             np.real(field[2,:,0,:])**2).transpose())
    
    plt.colorbar()
    plt.quiver(X[:,:], Z[:,:], np.real(field[0,:,0,:]).transpose(), np.real(field[2,:,0,:]).transpose())
    plt.show()
