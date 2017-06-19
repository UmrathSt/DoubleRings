# inductance of a circular wire-loop in front of a perfectly conducting plane
from scipy.special import ellipk, ellipe
from math import pi
from math import sqrt

def inductance(R1, R2, L, muR):
    mu0 = 4*pi*1e-7
    pref = mu0*muR*sqrt(L**2 + (R1 + R2)**2)
    prefK = (1 - 2*R1*R2/(L**2 + (R1 + R2)**2))
    arg = sqrt(4*R1*R2/(L**2 + (R1 + R2)**2))
    return pref*(prefK*ellipk(arg) - ellipe(arg))


if __name__ == "__main__":
    R1 = 9.8e-3
    R2 = 9.8e-3
    L = 2e-3
    print("Inductance for R1= %.2f mm, R2 = %.2f mm, L = %.2f mm = %.3e\n" %(R1*1e3, 
        R2*1e3, L*1e3, inductance(R1, R2, L, 1)))
