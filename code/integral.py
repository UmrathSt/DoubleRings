from legendre import jl, hl, Nlm, jdl, hdl, diffPlm, Plm
import numpy as np

def b_correction(kR, lmax):
    denom =   (Nlm(1,0)*diffPlm(1,0))**2*hl(1,kR)*jl(1,kR)
    factor = 1
    l = np.arange(2, lmax+1)
    l_e = l[::2]
    l_o = l[1::2]
    factor += np.array([(Nlm(l,0)*diffPlm(l,0))**2*hl(l,kR)*jl(l,kR) for l in l_o]).sum()
    factor -= np.array([(Nlm(l,0)*Plm(l,0)/kR)**2*hdl(l,kR)*jdl(l,kR) for l in l_e]).sum()

    return 1/factor



if __name__ == "__main__":
    lmax = 100
    result = np.zeros(lmax)
    for lmax in np.arange(1, lmax+1):
        result[lmax-1] = b_correction(0.1, lmax)
    print(result)
