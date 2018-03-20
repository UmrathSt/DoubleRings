import numpy as np
from matplotlib import pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--measdata", dest="measdata", type=str)
parser.add_argument("--simdata", dest="simdata", type=str)
parser.add_argument("--resultname", dest="resultname", type=str)
parser.add_argument("--xmin", dest="xmin", type=float)
parser.add_argument("--xmax", dest="xmax", type=float)

args = parser.parse_args()
xdim = 5
fig = plt.figure(figsize=(xdim,xdim/1.6))

meas = np.loadtxt(args.measdata, delimiter=",")
#sim  = np.loadtxt(args.sim, delimiter=",")

plt.plot(meas[:,0]/1e9, 1-np.abs(meas[:,1]+1j*meas[:,2])**2, "k-", linewidth=2)

plt.xlim([args.xmin, args.xmax])
plt.ylim([0,1])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel(r"$\mathbf{Frequency}\; \mathbf{[GHz]}$", fontsize = 14)
plt.ylabel(r"$\mathbf{Absorbed\; Power}\; \mathbf{[\%]}$", fontsize = 14)
fig.tight_layout()
fig.savefig(args.resultname, format="pdf")
