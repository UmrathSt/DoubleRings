import numpy as np
from matplotlib import pyplot as plt
from scipy.special import sph_harm 
from scipy.special import hankel1, jv
from math import sqrt

def Ylm(l, m, theta, phi):
    if abs(m) > l:
        return 0
    return sph_harm(m, l, phi, theta)

def zl(l, kr, reg):
    """ return regular (reg=1) or outgoing
        spherical bessel function (reg=0)
        of integer order l with argument kr
    """
    factor = np.sqrt(np.pi/(2*kr))
    if reg:
        return jv(l+0.5,kr)*factor
    if reg == 0:
        return hankel1(l+0.5, kr)*factor
    else:
        raise ValueError("Must be called with 0 or 1, but got", reg)

def zlD(l, kr, reg):
    return -l*zl(l, kr, reg) + kr * zl(l-1, kr, reg)


class plot_multipole_field:
    def __init__(self, l, m, k, x, y, z, reg):
            assert type(l) == type(m) and type(l) == int
            assert type(x) == type(y) and type(x) == np.ndarray
            assert type(z) == type(x) 
            self.l = l
            self.m = m
            self.k = k
            self.norm = np.sqrt(l*(l+1))
            self.XX, self.YY, self.ZZ = np.meshgrid(x, y, z, indexing = "ij")
            self.r = np.sqrt(self.XX**2 + self.YY**2 + self.ZZ**2)
            self.theta = np.arccos(self.ZZ/self.r)
            self.phi = np.arctan2(self.YY, self.XX)
            self.zl = zl(self.l, self.k*self.r, reg)
            self.zlD = zlD(self.l, self.k*self.r, reg)

    def get_Mlm(self):
        """ Get the values of the first-order vector
            spherical harmonic of angular momentum l
            and order m at the values given by the 
            two numpy arrays x and y at a fixed value
            of z.
        """
        theta = self.theta
        phi = self.phi
        r = self.r
        l = self.l
        m = self.m
        kr = self.k*self.r
        zl = self.zl    
        Etheta = 1j*m*Ylm(l,m, theta, phi)/np.sin(theta) * zl 
        Ephi = -(m / np.tan(theta) * Ylm(l, m, theta, phi) + 
                sqrt((l - m)*(l + m + 1)) * np.exp(-1j*phi) * 
                Ylm(l, m+1, theta, phi)) * zl
        Ex = Etheta * np.cos(theta) * np.cos(phi) - Ephi * np.sin(phi)
        Ey = Etheta * np.cos(theta) * np.sin(phi) + Ephi * np.cos(phi)
        Ez = -Etheta * np.sin(theta) 
        return [Ex/self.norm, Ey/self.norm, Ez/self.norm]
    
    def get_Nlm(self):
        """ Get the values of the second-order vector
            spherical harmonic of angular momentum l
            and order m at the values given by the 
            numpy arrays x, y, z
        """
        theta = self.theta
        phi = self.phi
        r = self.r
        l = self.l
        m = self.m
        kr = self.k * self.r
        zl = self.zl
        zlD = self.zlD
        Ephi = 1j*m*Ylm(l,m, theta, phi)/np.sin(theta) * zlD/kr
        Etheta = (m / np.tan(theta) * Ylm(l, m, theta, phi) + 
                sqrt((l - m)*(l + m + 1)) * np.exp(-1j*phi) * 
                Ylm(l, m+1, theta, phi)) * zlD/kr
        Er = Ylm(l, m, theta, phi)*sqrt(l*(l+1))*zl/kr
        Ex = (Er * np.sin(theta) * np.cos(phi) +
              Etheta * np.cos(theta) * np.cos(phi) +
              Ephi * (-1) * np.sin(phi)
                )
        Ey = (Er * np.sin(theta) * np.sin(phi) + 
                Etheta * np.cos(theta) * np.sin(phi) + Ephi * np.cos(phi)
                )
        Ez = Er * np.cos(theta) - Etheta * np.sin(theta)
        return [Ex/self.norm, Ey/self.norm, Ez/self.norm]



if __name__ == "__main__":
    N = 100
    val1 = 10 
    minval = 0 
    val2 =20 
    wavevector = 2
    x, z = np.linspace(-val1, val1, N), np.linspace(0, val2, N)
    y = np.array([0])
    lmax = 50 
    from matplotlib import pyplot as plt
    
    field = [0]*3
    for l in range(1, lmax+1):
        flm = plot_multipole_field(l=l, m=1, k=wavevector, x=x, y=y, z=z, reg=1)
        M, N = flm.get_Mlm(), flm.get_Nlm()
        for i in range(3):
            field[i] += (1j**(l+1)*np.sqrt((2*l+1)*np.pi)* (
                    2j*np.imag(M[i]) +  2*np.real(N[i])))
    plt.pcolor(flm.XX[:,0,:], flm.ZZ[:,0,:], np.sqrt(np.real(field[0][:,0,:])**2+
        np.real(field[1][:,0,:])**2+np.real(field[2][:,0,:])**2))
    plt.colorbar()
    skip_x, skip_z = 5,2
    plt.quiver(flm.XX[::skip_x,0,::skip_z], flm.ZZ[::skip_x,0,::skip_z], 
            np.real(field[0][::skip_x,0,::skip_z]), np.real(field[2][::skip_x,0,::skip_z]), color="Teal", headlength=6)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
