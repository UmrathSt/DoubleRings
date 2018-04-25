import numpy as np
from matplotlib import pyplot as plt

d1 = np.loadtxt("bistrcs_dielectriccube_1_eps_2.74.dat", delimiter=",")
d2 = np.loadtxt("Messdaten_bistatischer_RCS_22.5Grad.txt")

plt.plot(d1[1:])
plt.show()
