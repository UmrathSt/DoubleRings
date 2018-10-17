import numpy as np
from matplotlib import pyplot as plt

alpha = 22.5/180*np.pi
Adiag = 25e-4*np.abs((np.cos(alpha)+np.sin(alpha)))
lamb = 3e8/92.5e9
fRCS = 4*np.pi*Adiag**2/lamb**2

d1 = np.loadtxt("bistrcs_dielectriccube_22.5Grad_eps_2.74.dat", delimiter=",")
d2 = np.loadtxt("Messdaten_bistatischer_RCS_22.5Grad.txt")
plt.plot([0,360], [10*np.log10(fRCS)]*2, "k--", label=r"$4\pi\frac{A^2}{\lambda^2}$")
print("RCS: ", 10*np.log10(fRCS))
theta = np.linspace(-2,178,len(d2))
plt.plot(d2[:,0]-4,d2[::-1,1], "b-", label="Messung")
plt.legend(loc="best", fontsize=12).draw_frame(False)
theta = np.linspace(0,360,len(d1[1:]))
plt.plot(theta, 10*np.log10(d1[1:]), "ro-", linewidth=0.5, markersize=1, label="Simulation")
plt.xlabel(r"scattering angle [°]", fontsize=14)
plt.ylabel(r"RCS [dBsm]",fontsize=14)
plt.title(r"RCS eines PVC-Würfels, L=5 cm, $\varepsilon^\mathrm{sim}=2.74+0.026\mathrm{i}$")
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#plt.plot(theta,d2[::-1,1], "b-")
plt.savefig("bistatic_RCS_Sim_Mess.pdf", format="pdf")
plt.show()
