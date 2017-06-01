import matplotlib.pyplot as plt
from math import sin, cos, pi
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def endpoint(R, angle):
    x, y = cos(angle)*R, sin(angle)*R
    return x, y


R1, w1 = 9.8, 1.5
R2, w2 = 5.1, 0.5
SubstLz = 2
CopperLz= 0.5
UCDim = 10
tiltangle1 = 30 # degree
tiltangle2 = 5 # degree
fz = 17
hw = 0.4
hl = 0.9
circle1 = plt.Circle((0, 0), R1, color='k', alpha=0.4)
circle2 = plt.Circle((0, 0), R1-w1, color='w')
circle3 = plt.Circle((0, 0), R2, color='k', alpha=0.4)
circle4 = plt.Circle((0, 0), R2-w2, color='w')


fig = plt.figure(figsize=(8,8))
ax1 = fig.add_subplot(111)
ax1.set_xlim([-UCDim-0.5, UCDim+0.5])
ax1.set_ylim([-UCDim-0.5, UCDim+0.5])
ax1.add_artist(circle1)
ax1.add_artist(circle2)
ax1.add_artist(circle3)
ax1.add_artist(circle4)
X1, Y1 = endpoint(R1-1, tiltangle1*pi/180)
X11, Y11 = endpoint(R1-w1-1, tiltangle2*pi/180)
X2, Y2 = endpoint(R2-1, (180-tiltangle1)*pi/180)
X22, Y22 = endpoint(R2-w2-1, (180-tiltangle2)*pi/180)

ax1.text(X1*2/3, Y1*2/3+1, '$R_1$', fontsize=fz, rotation=tiltangle1)
ax1.text(X2*2/3, Y2*2/3+0., '$R_2$', fontsize=fz, rotation=-tiltangle1)
ax1.text(X11*2/3-1, Y11*2/3-1, '$R_1-w_1$', fontsize=fz, rotation=tiltangle2)
ax1.text(X22*2/3-1, Y22*2/3-1, '$R_2-w_2$', fontsize=fz, rotation=-tiltangle2)

ax1.arrow(0, 0, X1, Y1, head_width=0.5, head_length=hl, fc="k", ec="k")
ax1.arrow(0, 0, X2, Y2, head_width=0.5, head_length=hl, fc="k", ec="k")
ax1.arrow(0, 0, X11, Y11, head_width=0.5, head_length=hl, fc="k", ec="k")
ax1.arrow(0, 0, X22, Y22, head_width=0.5, head_length=hl, fc="k", ec="k")

ax1.set_xticks([])
ax1.set_yticks([])
ax1.plot([-UCDim, -UCDim, UCDim, UCDim], [-UCDim, UCDim, UCDim, -UCDim], "k-")
ax1.arrow(-UCDim, -UCDim, 2*UCDim-1, 0, head_width=hw, head_length=hl, fc="k", ec="k")
ax1.arrow(UCDim, -UCDim, -2*UCDim+1, 0, head_width=hw, head_length=hl, fc="k", ec="k")
ax1.text(0, -UCDim-1.5, "$L^\mathrm{UC}_\mathrm{x}$", fontsize=fz)
ax1.arrow(-UCDim, -UCDim, 0, 2*UCDim-1, head_width=hw, head_length=hl, fc="k", ec="k")
ax1.arrow(-UCDim, UCDim, 0, -2*UCDim+1, head_width=hw, head_length=hl, fc="k", ec="k")
ax1.text(-UCDim-1.85, 0, "$L^\mathrm{UC}_\mathrm{y}$", fontsize=fz)
ax1.set_aspect("equal")
ax1.axis("off")

fig.savefig("double_ring_sketch.pdf", format="pdf")
plt.show()
