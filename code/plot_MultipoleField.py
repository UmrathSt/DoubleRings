import numpy as np
from matplotlib import pyplot as plt
from scipy.special import sph_harm 
from scipy.special import hankel1, jv
from math import sqrt

def Ylm(l, m, theta, phi):
    return sph_harm(m, l, theta, phi)

def get_Mlm(l, m, kr, x, y, z, reg):
    """ Get the values of the first-order vector
        spherical harmonic of angular momentum l
        and order m at the values given by the 
        numpy arrays x, y, z
    """
    r = np.sqrt(x**2+y**2+z**2)
    theta = np.arccos(z/r)
    phi = np.arctan2(y, x)
    size = x.size
    assert size == y.size and size == z.size
    r = r.reshape(size, 1, 1)
    theta = theta.reshape(1, size, 1)
    phi = phi.reshape(1, 1, size)
    sph_factor = np.sqrt(np.pi/(2*kr))
    if reg == 1:
        zl = jv(l+0.5, kr) * sph_factor
    if reg == 0:
        zl = hankel1(l+0.5, kr) * sph_factor 

    Etheta = 1j*m*Ylm(l,m, theta, phi)/np.sin(theta) * zl
    Ephi = -(m / np.tan(theta) * Ylm(l, m, theta, phi) + 
            sqrt((l - m)*(l + m + 1)) * np.exp(-1j*phi)*Ylm(l, m+1, theta, phi)) * zl
    Ex = Etheta * np.cos(theta) * np.cos(phi) - Ephi * np.sin(phi)
    Ey = Etheta * np.cos(theta) * np.sin(phi) + Ephi * np.cos(phi)
    Ez = -Etheta * np.sin(theta)
    return [Ex, Ey, Ez]

def get_Nlm(l, m, kr, x, y, z, reg):
    """ Get the values of the second-order vector
        spherical harmonic of angular momentum l
        and order m at the values given by the 
        numpy arrays x, y, z
    """
    r = np.sqrt(x**2+y**2+z**2)
    theta = np.arccos(z/r)
    phi = np.arctan2(y, x)
    size = x.size
    assert size == y.size and size == z.size
    r = r.reshape(size, 1, 1)
    theta = theta.reshape(1, size, 1)
    phi = phi.reshape(1, 1, size)
    sph_factor = np.sqrt(np.pi/(2*kr))
    if reg == 1:
        zl = jv(l+0.5, kr) * sph_factor
        zlD= -l * zl + kr*jv(l-0.5, kr)*sph_factor
    if reg == 0:
        zl = hankel1(l+0.5, kr) * sph_factor
        zlD= -l * zl + kr*hankel1(l-0.5, kr)*sph_factor

    Ephi = 1j*m*Ylm(l,m, theta, phi)/np.sin(theta) * zlD
    Etheta = (m / np.tan(theta) * Ylm(l, m, theta, phi) + 
            sqrt((l - m)*(l + m + 1)) * np.exp(-1j*phi)*Ylm(l, m+1, theta, phi)) * zlD
    Er = Ylm(l, m, theta, phi)
    Ex = (Er * np.sin(theta) * np.cos(phi) +
            Etheta * np.cos(theta) * np.cos(phi) - Ephi * np.sin(phi)
            )
    Ey = (Er * np.sin(theta) * np.sin(phi) + 
            Etheta * np.cos(theta) * np.sin(phi) + Ephi * np.cos(phi)
            )
    Ez = Er * np.cos(theta) - Etheta * np.sin(theta)
    return [Ex, Ey, Ez]



if __name__ == "__main__":
    print(get_Nlm(10, 3, 2, np.array([0.1, 0.2, 0.3]), np.array([0.1, 0.2, 0.3]), np.array([0.001, 0.001, 0.001]), reg = 1))

