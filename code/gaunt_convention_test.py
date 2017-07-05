# coding: utf8
# test the sign conventions used for the
# gaunt coefficients
from sympy.physics.wigner import wigner_3j
from math import sqrt
from cython_gaunt import gaunt
from math import pi
import numpy as np

def wigner_gaunt(l1, l2, m):
    """ integral of three spherical harmonics
        int Y_lm Y_pn Y_qo 
    """
    pref = sqrt((2*l1 + 1)*(2*l2 + 1)/(4*pi))
    return np.array([pref*sqrt(2*lpp + 1)*float(wigner_3j(l1,l2,lpp,m,-m,0)*wigner_3j(l1,l2,lpp,0,0,0))
                    for lpp in range(abs(l1-l2), l1+l2+1, 2)])





if __name__ == "__main__":
    t1 = [(2,5,1), (2,5, 2),
          (2,4,1), (2,4, 2),
          (1,4,2), (1,4, 1)]
    for t in t1:
        print(wigner_gaunt(t[0], t[1], t[2]))
        print(gaunt(t[0], t[1], t[2]))
