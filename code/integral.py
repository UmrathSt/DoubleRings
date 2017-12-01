from legendre import jl, hl, Nlm, jdl, hdl, diffPlm, Plm
import numpy as np

def b_correction(kR, l1, l2, m1, lmax):
    if not (l1+m1)%2:
        return 0
    denom =   (Nlm(l1, m1)*diffPlm(l1, m1)*Nlm(l2,m1)*diffPlm(l2,m1))*hl(l1,kR)*jl(l2,kR)
    factor = 1
    ls = np.arange(max(1,abs(m1)), lmax+1)
    factor += np.array([(Nlm(l, m1)*diffPlm(l, m1))**2*hl(l, kR)*jl(l,kR) for l in ls]).sum()
    factor += np.array([(Nlm(l, m1)*Plm(l, m1)/kR*m1)**2*hdl(l, kR)*jdl(l,kR) for l in ls]).sum()
    factor /= denom
    return 1/factor

def a_correction(kR, l1, l2, m1, lmax):
    if (l1+m1)%2 or m1 == 0:
        return 0
    denom = (Nlm(l1, m1)*Plm(l1, m1)*Nlm(l2,m1)*Plm(l2, m1)/kR**2*m1**2)*hdl(l1,kR)*jdl(l2,kR)
    factor = 1
    ls = np.arange(max(1,abs(m1)), lmax+1)
    factor += np.array([(Nlm(l, m1)*diffPlm(l, m1))**2*hl(l, kR)*jl(l,kR) for l in ls]).sum()
    factor += np.array([(Nlm(l, m1)*Plm(l, m1)/kR*m1)**2*hdl(l, kR)*jdl(l,kR) for l in ls]).sum()
    factor /= denom
    return 1/factor


if __name__ == "__main__":
    from matplotlib import pyplot as plt
    plt.rc("text", usetex=True)
    LMAX = 100
    lmax = 21
    kr = 0.1 
    mTE, mTM = 0, 1 
    l = np.arange(1, lmax+1)
    result_a = np.zeros(lmax, dtype=np.complex128)
    result_b = np.zeros(lmax, dtype=np.complex128)

    for L in l:
        result_a[L-1] = a_correction(kr, L, L, mTM, LMAX) 
        result_b[L-1] = b_correction(kr, L, L, mTE, LMAX) 
    plt.plot(l, np.abs(result_a), linestyle="", markeredgewidth="1.5", markersize=8,
            label="$a_\ell$", markeredgecolor="r", marker="o", markerfacecolor="None")
    plt.plot(l, np.abs(result_b), linestyle="", markeredgewidth="1.5", markersize=8,
            label="$b_\ell$", markeredgecolor="b", marker="x", markerfacecolor="None")
    plt.xlabel("$\ell$", fontsize=14)
    plt.ylabel("$a^\mathrm{ring}_\ell/a^\mathrm{sphere}_\ell, \
                 b^\mathrm{ring}_\ell/b^\mathrm{sphere}_\ell$", fontsize=14)
    plt.legend(loc="best", fontsize=14).draw_frame(False)
    plt.show()
