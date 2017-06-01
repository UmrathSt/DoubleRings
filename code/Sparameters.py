import numpy as np
from matplotlib import pyplot as plt
from filewalk import fileList
import sys
from matplotlib import rc
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
majorLocator = MultipleLocator(5)
minorLocator = MultipleLocator(1)

rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


argv = sys.argv
print("lenargv = ", len(argv))
if not len(argv) == 4:
    raise ValueError("Falsche Anzahl an Parametern!\n")
sweep_type = argv[3]
if sweep_type == "w1":
    w1 = "changing"
else:
    w1 = "0.5"
if sweep_type == "w2":
    w2 = "changing"
else:
    w2 = "1.5"

source = argv[1]
beginswith = argv[2]
filenames = fileList(source, beginswith) 
# current and voltage as functions of time
fig = plt.figure()
ax1 = fig.add_subplot(111)

cm = plt.get_cmap("rainbow")
col = [cm(val) for val in np.linspace(0,1, len(filenames))]
print(filenames)
data = [np.loadtxt(file_, delimiter=",") for file_ in filenames]
S_parameters = [d[:,1]+1j*d[:,2] for d in data]
minS11 = -10
wdct = {"w1": r"w_1", "w2": r"w_2"}
for i, S in enumerate(S_parameters):
    s11db = 20*np.log10(abs(S))
    minS11 = min(minS11, np.min(s11db))
    label = filenames[i][45:-3]
    if sweep_type == "w2":
        Rlong = label.split(r"R2_", 1)[1]
        R = Rlong[0:3]
        w = Rlong.split(r"%s_" %sweep_type, 2)[1][0:3]
        R2str = r"R_2=%s, " %R
    else:
        w = label.split(r"%s_" %sweep_type, 2)[1][0:3]
        R2str = ""
    ax1.plot(data[i][:,0]/1e9, s11db, color=col[i],linestyle="-", linewidth=2, label = r"$%s%s=%s$ mm" %(R2str,wdct[sweep_type], w))

if sweep_type == "w1":
    ncol = 3
else:
    ncol = 2
ax1.legend(ncol = ncol,loc = "best", fontsize=12).draw_frame(False)
ax1.set_xlabel(r"$f\, \mathrm{[GHz]}$", fontsize=14)
ax1.set_ylabel(r"$20\,\log(|S_{11}|)$", fontsize=14)
ax1.set_title("parameter study for concentric ring absorber", fontsize=14)
ax1.set_ylim([1.1*minS11,0])
annotate_string = r"parameters:"
annotate_string += "\n"
annotate_string += r"square unit-cell, $L_\mathrm{UC}=20$ mm"
annotate_string += "\n"
annotate_string += r"ring 1: $R_1 = 9.8$ mm"
if not sweep_type == "w1":
    annotate_string += r", $w_1=1.5$ mm"

if not sweep_type == "w2":
    annotate_string += "\n"
    annotate_string += r"ring 2: $R_2 = 5.1$ mm"
if not sweep_type == "w2":
    annotate_string += r", $w_2=0.5$ mm"
annotate_string += "\n"
annotate_string += r"FR4 substrate $l_z=2$ mm"
ax1.annotate(annotate_string, fontsize = 14, xy=(6,-20))
ax1.xaxis.set_major_locator(majorLocator)
ax1.xaxis.set_minor_locator(minorLocator)
ax1.tick_params(axis="both", which="major", labelsize=14)
plt.tight_layout()
plt.savefig("dual-wifi_absorber_%s.pdf" %sweep_type, format="pdf")
plt.show()
