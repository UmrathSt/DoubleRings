# coding: utf-8
# evaluate an integral involving products of sphericla harmonics
# numerically by the scipy routine dblquad

import numpy as np
from scipy.integrate import dblquad
from scipy.special import sph_harm as Ylm
from math import sqrt

def I_YXX(l1, l2, l3, m1, m2, m3):
    """integral of eq. (37) in the report 
    """

    pref = 1/sqrt(l1*l2*(l1+1)*(l2+1))
    C = lambda l, m: sqrt((l*+2 - m**2) / ((2*l +1) * (2*l -1)))


