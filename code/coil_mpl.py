# coding: utf-8

import matplotlib.pyplot as plt
import matplotlib.patches as patches

def rect(xc, yc, lx, ly, color="black"):
    patch = patches.Rectangle((xc-lx/2, yc-ly/2), lx, ly, facecolor=color)
    return patch

fig1 = plt.figure()
ax1 = fig1.add_subplot(111, aspect='equal')


L0 = 0.3
w = 0.025
dL = 0.05

Uxc, Uyc = 0.5, 0.5
xc, yc = -0.4, 0
g0 = 0.2 # initial gap between the bar and 
         # the edge of the unit-cell


for i in range(4):
    s = (-1)**i
    Lx = w*(0.5*(1-s)) + L0*(0.5*(1+s))
    Ly = L0*(0.5*(1-s)) + w*(0.5*(1+s))
    yc, xc = s*0.5*(s+1)*L0, s*0.5*(1-s)*L0 
    ax1.add_patch(rect(xc + Uxc, yc + Uyc, s*Lx, s*Ly))
    L0 -= i//3*dL
    

plt.xlim([-2,2])
plt.ylim([-2,2])
plt.show()
