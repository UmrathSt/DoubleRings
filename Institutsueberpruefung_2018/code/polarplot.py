from matplotlib import pyplot as plt
import numpy as np

d = np.loadtxt("rcs_PVC_cube_L_5cm_alpha_22.5.dat", delimiter=",")
alpha = d[:,0]
rcs = d[:,1]
xdim = 6
fig = plt.figure(figsize=(xdim,xdim/1.6))
ax = plt.subplot(111, polar=True)
ax.plot(alpha/180*np.pi, 10*np.log10(rcs),"k-", linewidth=2)
ax.set_xlim([0,np.pi])
ax.set_ylim([-50, 24])
ax.set_xlabel(r"$\mathbf{Scattering\; Angle\; [^\circ]}$", fontsize=14)
ax.set_ylabel(r"$\mathbf{RCS\; [dBsm]}$", fontsize=14)
ax.yaxis.set_label_coords(-0.15,0.5)
ax.set_rlabel_position(292.5)
ax.set_thetagrids(np.arange(0,360,45), frac=1.15)
ax.set_theta_zero_location("N")
ax.set_yticks([-2,12])
plt.xticks(fontsize=12)
plt.tight_layout()
plt.show()
