from matplotlib import pyplot as plt
import argparse
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument("--eps", dest="eps", type=float)
parser.add_argument("--kappa", dest="kappa", type=float)
parser.add_argument("--lz", dest="lz", type=float)
parser.add_argument("--fmin", dest="fmin", type=float, default=0)
parser.add_argument("--fmax", dest="fmax", type=float, default=1e18)
parser.add_argument("--filename", dest="filename", type=str)
parser.add_argument("--onlyZm", dest="onlyZm", type=bool, default=False)


def get_Zsubs(eps,d,f):
    """ Return the impedance of a dielectric slab
        with permittivity eps, thickness d at frequency
        f
    """
    w = 2*np.pi*f
    return -1j*Z0/np.sqrt(eps)*np.tan(w*np.sqrt(eps)/C0)

def get_Zfss(Z0, Zd, Gamma0):
    """ Return the impedance of a frequency selective surface
        which is modelled in parallel to an impedance Zsubs
    """
    r = (1-Gamma0)/(1+Gamma0)
    return Z0*Zd/(r*Zd-Z0) 


if __name__ == "__main__":
    args = parser.parse_args()
    kappa = args.kappa
    d = np.loadtxt(args.filename, delimiter=",")
    f, Gamma0 = d[:,0], d[:,1]+1j*d[:,2]
    fmin, fmax = args.fmin, args.fmax
    i_min = (np.abs(f-fmin)).argmin() 
    i_max = (np.abs(f-fmax)).argmin() 
    f = f[i_min:i_max]
    Gamma0 = Gamma0[i_min:i_max]
    Z0 = 376.73
    C0 = 299792458
    EPS0 = 8.85e-12
    eps = args.eps
    lz = args.lz
    eps = eps+kappa*1j/(2*np.pi*f*EPS0)
    Zd = get_Zsubs(eps, args.lz, f)
    Zm = Z0*(1+Gamma0)/(1-Gamma0)
    Zfss = Zm*Zd/(Zd-Zm)
    print(Zm)
    if not args.onlyZm:

        plt.plot(f/1e9, np.real(Zfss),"r-",linewidth=2,
              label="$\mathcal{R}(Z_\mathrm{fss})$")
        plt.plot(f/1e9, np.imag(Zfss),"b-",linewidth=2,
              label="$\mathcal{I}(Z_\mathrm{fss})$")
        plt.plot(f/1e9, np.real(Zd),"k-",linewidth=2,
            label="$\mathcal{R}(Z_\mathrm{d})$")
        plt.plot(f/1e9, np.imag(Zd),"k--",linewidth=2,
            label="$\mathcal{I}(Z_\mathrm{d})$")
    plt.plot(f/1e9, np.real(Zm),"m-",linewidth=2,
            label="$\mathcal{R}(Z_\mathrm{m})$")
    plt.plot(f/1e9, np.imag(Zm),"c-",linewidth=2,
            label="$\mathcal{I}(Z_\mathrm{m})$")
    plt.xlabel("f [GHz]")
    plt.ylabel("$\mathcal{R}(Z_\mathrm{fss}),\mathcal{I}(Z_\mathrm{fss}),\mathcal{R}(Z_\mathrm{m}),\mathcal{I}(Z_\mathrm{m}$")
    plt.xlim(f[[0,-1]]/1e9)
    plt.legend(loc="best").draw_frame(False)
    plt.show()
