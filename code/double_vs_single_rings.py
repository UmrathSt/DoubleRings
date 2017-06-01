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


# current and voltage as functions of time
fig = plt.figure()
ax1 = fig.add_subplot(111)
data_folder = "../data"
filenames = [data_folder+"/S_f_UCDim_20_lz_2_R1_0_w1_0_R2_5.1_w2_0.5_sweep.txt", data_folder+"/S_f_UCDim_20_lz_2_R1_9.8_w1_1.5_R2_0_w2_0_sweep.txt", data_folder+"/R1_9.8_R2_5.1_width_sweep/w1_sweep/S_f_UCDim_20_lz_2_R1_9.8_w1_1.5_R2_5.1_w2_0.5_sweep.txt"] 

cm = plt.get_cmap("rainbow")
col = ["b", "r", "k"]
lstyle = ["-", "-", "--"]
lw = [2, 2, 1]
lab = [r"only ring 2", r"only ring 1", r"both rings"]
data = [np.loadtxt(file_, delimiter=",") for file_ in filenames]
S_parameters = [d[:,1]+1j*d[:,2] for d in data]
minS11 = -10
for i, S in enumerate(S_parameters):
    s11db = 20*np.log10(abs(S))
    minS11 = min(minS11, np.min(s11db))
    label = filenames[i][45:-3]
    ax1.plot(data[i][:,0]/1e9, s11db, color=col[i],linestyle=lstyle[i], linewidth=lw[i], label = lab[i] )

ax1.legend(loc = "best", fontsize=14).draw_frame(False)
ax1.set_xlabel(r"$f\, \mathrm{[GHz]}$", fontsize=14)
ax1.set_ylabel(r"$20\,\log(|S_{11}|)$", fontsize=14)
ax1.set_title("parameter study for concentric ring absorber", fontsize=14)
ax1.set_ylim([1.1*minS11,0])
annotate_string = r"parameters:"
annotate_string += "\n"
annotate_string += r"square unit-cell, $L_\mathrm{UC}=20$ mm"
annotate_string += "\n"
annotate_string += r"ring 1: $R_1 = 9.8$ mm"
annotate_string += r", $w_1=1.5$ mm"
annotate_string += "\n"
annotate_string += r"ring 2: $R_2 = 5.1$ mm"
annotate_string += r", $w_2=0.5$ mm"
annotate_string += "\n"
annotate_string += r"FR4 substrate $l_z=2$ mm"
ax1.annotate(annotate_string, fontsize = 14, xy=(6,-20))
ax1.xaxis.set_major_locator(majorLocator)
ax1.xaxis.set_minor_locator(minorLocator)
ax1.tick_params(axis="both", which="major", labelsize=14)
plt.tight_layout()
plt.savefig("wifi_absorber_single_double_rings.pdf", format="pdf")
plt.show()
