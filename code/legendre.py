# coding: utf-8
# Associated Legendre Polynomials at argument zero
# and there derivative at argument zero
from scipy.special import gamma
from scipy.misc import factorial
from numpy import sin, cos, sqrt, pi
import numpy as np
from scipy.special import jv, hankel1
from scipy.special import gammaln as lngamma

jl = lambda l, x: sqrt(pi/(2*x)) * jv(l+0.5,x)
hl = lambda l, x: sqrt(pi/(2*x)) * hankel1(l+0.5,x)
jdl = lambda l, x: -l*jl(l,x) + x * jl(l-1,x)
hdl = lambda l, x: -l*hl(l,x) + x * hl(l-1,x)
Nlm = lambda l, m : sqrt((2*l+1)*pi/(l*(l+1))*np.exp(lngamma(l-abs(m)+1)-lngamma(l+abs(m)+1)))


def diffPlm(l, m):
    """ Return the value of an associated Legendre polynomial
        Plm at argument zero
        Parameters
        ----------
        l, m, (integer) : degree l, order m

        Returns
        -------
        Plm(0) (float) 
    """
    if not (l+m)%2:
        return 0
    else:
     return -sin((l+m)*pi/2) * np.exp(np.log(2)*(m+1) + 
                    lngamma((l+m)/2+1)-lngamma((l-m+1)/2)) /sqrt(pi)

def Plm(l, m):
    """ Return the value of the deriviative of the associated Legendre polynomial
        Plm at argument zero
        Parameters
        ----------
        l, m, (integer) : degree l, order m

        Returns
        -------
        derivative(Plm(x), x=0) (float) 
    """
    if (l+m)%2:
        return 0
    else:
        return cos((l+m)*pi/2) * np.exp(np.log(2)*m + lngamma((l+m+1)/2)- lngamma((l-m)/2+1))/sqrt(pi)

def ringAl(l, m, kr):
    """ Return the diagonal TM T-Matrix element of
        a ring scatterer of radius r at wavevector k
    """
    factor = float(sqrt(2*pi/(l*(l+1)))*Plm(l,m) )
    factor1= float((2*pi/(l*(l+1))) * Plm(l, m))
    return factor1 * Nlm(l, m)**2

def ringBl(l, m, kr):
    """ Return the diagonal TE T-Matrix element of
        a ring scatterer of radius r at wavevector k
    """
    factor = float(sqrt(2*pi/(l*(l+1)))*Plm(l,m) )
    factor1= float(sqrt(2*pi/(l*(l+1))) * Plm(l, m))
    return factor1 * Nlm(l, m)**2

def testfactor(l, m):
    """ Trial to guess the amplitude ratio of TE and TM scattering
        coefficients of ring/sphere
    """
    fac1 = (2*l+1)/(2*l*(l+1))
    fac2 = factorial(l-abs(m))/factorial(l+abs(m))
    fac3 = 2**(2*m+2)/pi
    fac4 = (gamma((l+m)/2+1)/gamma((l-m+1)/2))**2
    return fac1*fac2*fac3*fac4

if __name__ == "__main__":
    kr = 0.1 
    for l in range(1, 11):
        print("l=%i" %l)
        m = 1 
        print("Testfaktor: ", 0.526*3*testfactor(l,m)/(2*l+1))
        print("factor(l=%i, m=%i) = %.2e" %(l,m, Plm(l,m)/(diffPlm(l,m)+Plm(l,m))))
