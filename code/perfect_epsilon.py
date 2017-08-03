import numpy as np
from  min_helpers import get_min_mask, get_plot_data
from filewalk import fileList
from matplotlib import pyplot as plt

filenames = fileList("/home/stefan/Arbeit/latex/DoubleRings/code/double_ring_eps_sweep/tand_0.015/", "S11", ".txt")
datasets = get_plot_data(filenames, "eps_")
measured = np.loadtxt("/home/stefan/Arbeit/latex/DoubleRings/code/"+
        "Messung_S11_Re_Im_1-18_GHz.txt")[:,0:2]
measured[:,1] = 10**(measured[:,1]/20)
#print(datasets)

#plt.plot(measured[:,0], measured[:,1])
Mmask = get_min_mask(measured[:,1], 0.5)
#plt.plot(measured[mask,0], measured[mask,1], "ko")
#plt.show()
fig = plt.figure()
ax = fig.add_subplot(111)
measured_f = measured[Mmask,0]
colors = ["r","b","g","m","c"]
counter = 0 # count datasets in order to only get one lable indicating
            # the resonance frequency
for dataset in datasets:
    eps = dataset[0]
    mask = get_min_mask(dataset[1][:,1], 0.5)
    f, S = dataset[1][mask,0], dataset[1][mask,1]
    for idx, f in enumerate(f):
        label = r"$f_%i=%.2f$ GHz" %(idx+1, f/1e9)
        if counter > 0:
            label = ""
        if not eps == 4.45:
            ax.plot(eps, f/1e9/measured_f[idx], 
               "ko", color=colors[idx], label=label)
    counter += 1
#ax.set_ylim([0, 7])
#no_peaks = len(measured[Mmask,0])
#for i in range(no_peaks):
#    ax.plot([4, 4.6], measured[Mmask,0][i]*np.ones(2), "k--", label="measured")
ax.set_xlim([4, 4.6])
ax.set_ylim([0.9,1.1])
ax.set_xlabel(r"$\epsilon_\mathrm{r}^\mathrm{FR4}$", fontsize=16)
ax.set_ylabel(r"$f_i^\mathrm{sim.} / f_i^\mathrm{meas.} $", fontsize=16)
ax.tick_params(axis="both", labelsize=16)
ax.legend(loc="best").draw_frame(False)
#plt.show()
fig.savefig(r"fi_sim_meas.pdf", format="pdf")
