import numpy as np
from matplotlib import pyplot as plt
from scipy.special import sph_harm 
from scipy.special import hankel1, jv
from math import sqrt

def Ylm(l, m, theta, phi):
    if abs(m) > l:
        return 0
    return sph_harm(m, l, theta, phi)

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
            assert type(z) == float
            self.l = l
            self.m = m
            self.k = k
            XX, YY = np.meshgrid(x, y)
            self.r = np.sqrt(XX**2 + YY**2 + z**2)
            self.theta = np.arccos(z/self.r)
            self.phi = np.arctan2(YY, XX)
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
        return [Ex, Ey, Ez]
    
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
                Etheta * np.cos(theta) * np.cos(phi) - Ephi * np.sin(phi)
                )
        Ey = (Er * np.sin(theta) * np.sin(phi) + 
                Etheta * np.cos(theta) * np.sin(phi) + Ephi * np.cos(phi)
                )
        Ez = Er * np.cos(theta) - Etheta * np.sin(theta)
        return [Ex, Ey, Ez]



if __name__ == "__main__":
    N = 100 
    val = 10 
    x, y = np.linspace(-val, val, N), np.linspace(-val, val, N)
    z = 0.001
   
    from matplotlib import pyplot as plt

    f = plot_multipole_field(l=1, m=1, k=2, x=x, y=y, z=z, reg=1)
    XX, YY = np.meshgrid(x, y)
    for l in range(1, 15):
        for m in [-1, 1]:
            flm = plot_multipole_field(l=l, m=m, k=2, x=x, y=y, z=z, reg=1)
            M, N = flm.get_Mlm(), flm.get_Nlm()
            field = []
            for i in range(3):
                field.append(1j**(l+1)*np.sqrt((2*l+1)*np.pi)* (
                    M[i] +m * N[i]))

    plt.pcolor(XX, YY, np.real(field[0]))
    plt.colorbar()
    plt.show()
