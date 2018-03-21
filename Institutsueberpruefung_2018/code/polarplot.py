from matplotlib import pyplot as plt
import numpy as np

d = np.loadtxt("rcs_PVC_cube_L_5cm_alpha_22.5.dat", delimiter=",")
m = np.loadtxt("RCS_PVC_Bistatisch.dat")
alpha = d[:,0]
rcs = d[:,1]
xdim = 6
fig = plt.figure(figsize=(xdim,xdim/1.6))
ax = plt.subplot(111, polar=True)
ax.plot(alpha/180*np.pi, 10*np.log10(rcs),"k-", linewidth=2,
        label="sim.")
sangle = (m[:,0])/180*np.pi

# find the angles between 90 and 180 degree
ww3 = sangle < np.pi/2
ax.plot(-sangle[ww3]-np.pi, m[ww3,1], "r-", linewidth=1)
# find the angles between 0 and 90 degree
ww2 = sangle >= np.pi/2
ax.plot(-sangle[ww2]+np.pi, m[ww2,1], "r-", linewidth=1,
        label="meas.")

ax.set_xlim([0,np.pi])
ax.set_ylim([-50, 24])
ax.set_xlabel(r"$\mathbf{Scattering\; Angle\; [^\circ]}$", fontsize=14)
ax.set_ylabel(r"$\mathbf{RCS\; [dBsm]}$", fontsize=14)
ax.yaxis.set_label_coords(-0.15,0.5)
ax.set_rlabel_position(292.5)
ax.set_thetagrids(np.arange(0,360,45), frac=1.15)
ax.set_theta_zero_location("N")
ax.set_yticks([-2,12])
plt.legend(loc=9,bbox_to_anchor=(1.3,1.2))
plt.xticks(fontsize=12)
plt.tight_layout()
fig.savefig("scattering_plot.pdf", format="pdf")
plt.show()
