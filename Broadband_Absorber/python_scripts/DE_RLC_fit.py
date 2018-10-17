import numpy as np
import sys
hdir = "/home/stefan_dlr/Arbeit/openEMS/metamaterials/python_scripts"
sys.path.append(hdir)
from tl_fit import Rresiduals, Zresiduals, fit_func
from scipy.optimize import differential_evolution
from matplotlib import pyplot as plt
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--filename", dest="fname", type=str)
parser.add_argument("--fmin", dest="fmin", type=float)
parser.add_argument("--fmax", dest="fmax", type=float)
parser.add_argument("--plot", dest="plot", type=bool, default=False)
args = parser.parse_args()

data = np.loadtxt(args.fname, delimiter=",")
f, Rnum = data[:,0], data[:,1]+1j*data[:,2]
eps, tand, D = 4.6, 0.02, 0.15e-3
I = (np.abs(f-args.fmin)).argmin()
J = (np.abs(f-args.fmax)).argmin()


sumRresiduals = lambda x: Rresiduals(x, Rnum[I:J], f[I:J], D, eps, tand).sum()
solution = differential_evolution(sumRresiduals, bounds=[(0,1e1),(1e-13,1e-10),(1e-14, 1e-10)])
print(solution)
Znum = 376*(1+Rnum)/(1-Rnum)
Zfit = fit_func(solution.x, f, D, eps, tand)
Rfit = (Zfit-376)/(Zfit+376)
num = Rnum
fit = Rfit
plt.plot(f/1e9, np.abs(num), "b-", label="$\mathcal{R}(Z_\mathrm{num})$")
plt.plot(f/1e9, np.abs(fit), "m--", label="$\mathcal{R}(Z_\mathrm{fit})$")
#plt.plot(f/1e9, num.real, "r-", label="$\mathcal{R}(Z_\mathrm{num})$")
#plt.plot(f/1e9, num.imag, "b-", label="$\mathcal{R}(Z_\mathrm{num})$")
#plt.plot(f/1e9, fit.real, "m--", label="$\mathcal{R}(Z_\mathrm{fit})$")
#plt.plot(f/1e9, fit.imag, "c--", label="$\mathcal{R}(Z_\mathrm{fit})$")
plt.legend(loc="best").draw_frame(False)
plt.xlim([f[0]/1e9, f[-1]/1e9])
plt.xlabel("f [GHz]")
plt.ylabel("$\mathcal{R}(Z),\, \mathcal{I}(Z_\mathrm{fit})$")
plt.show()
