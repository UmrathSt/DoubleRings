# coding: utf-8
# implementation of the factor occurring in an
# xy-periodic array of scatterers in a spherical
# wave basis
import sympy as sy
from sympy.functions.special.spherical_harmonics import Ynm as Ylm
from sympy.functions.special.bessel import jn
theta, phi, r, l = sy.symbols(["theta", "phi", "r", "l"])


def periodic_xy(Lxy, omega, lmax):
    coef = lambda l_, m_: (4*sy.pi*sy.I*jn(l_, 2*sy.pi*r/Lxy)*
                    Ylm(l_, m_, theta, phi)*sy.conjugate(
                        Ylm(l_, m_, sy.pi/2, sy.pi/4)))
    result = 0
    for ll in range(0, lmax+1):
        for mm in range(0, ll+1):
            result += coef(ll, mm)
    return result





if __name__ == "__main__":

    print("result:\n", periodic_xy(1, 1, 1).expand(func=True).simplify())
