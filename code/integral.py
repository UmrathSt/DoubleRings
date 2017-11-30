from legendre import jl, hl, Nlm, jdl, hdl, diffPlm, Plm
import numpy as np

def b_correction(kR, l1, m1, lmax):
    denom =   (Nlm(l1, m1)*diffPlm(l1, m1))**2*hl(l1,kR)*jl(l1,kR)
    factor = 1
    ls = np.arange(max(1,abs(m1)), lmax+1)
    l_e = ls[1::2]
    l_o = ls[::2]
    factor += np.array([(Nlm(l, m1)*diffPlm(l, m1))**2*hl(l, kR)*jl(l,kR) for l in l_o]).sum()
    factor += np.array([(Nlm(l, m1)*Plm(l, m1)/kR*m1)**2*hdl(l, kR)*jdl(l,kR) for l in l_e]).sum()
    factor /= denom
    return 1/factor



if __name__ == "__main__":
    LMAX = 200
    kr = 1 
    m1 = 0 
    l1 = 1 
    from matplotlib import pyplot as plt
    l = np.arange(1, LMAX+1)
    print(l)
    result = np.zeros(LMAX, dtype=np.complex128)
    for lmax in l:
        print("filling idx ", lmax-1)
        result[lmax-1] = b_correction(kr, l1, m1, lmax)
    plt.plot(l, np.abs(result), "ko")
    print(result)
    plt.show()
