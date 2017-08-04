# coding: utf-8
# get the dependence of the depth and position of the
# absorption frequencies as a function of the real and imaginary part
# of the substrate permittivity for a double-ring geometry

import numpy as np
from matplotlib import pyplot as plt
from filewalk import fileList
from scipy.signal import argrelmin, general_gaussian, fftconvolve
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)




def get_minima_positions(dataset, above, N=21, P=0.5, Sig=25):
    """ get the position of the local minima
        in the 1d numpy array >dataset< after
        smoothing it with a gaussian
        Returns: indices of the local minima,
        where the value of dataset[indices] is 
        smaller than >above<
    """
    window = general_gaussian(N, p=P, sig=Sig)
    filtered = fftconvolve(window, dataset)
    filtered = (np.average(dataset) / np.average(filtered)) * filtered
    filtered = np.roll(filtered, int((N-1)/2))
    min_pos = argrelmin(filtered)[0][1:-1]
    mask = min_pos[filtered[min_pos] < above]
    return mask

def get_plot_data(file_list, val):
    """ get a list of tuples corresponding to the 
        value of > val < and a 2d numpy array holding
        the values for the abscissa and the ordinate
    """
    if val == "eps_":
        string = r"$\epsilon_\mathrm{r}=%.2f$"
        endmark = "_"
    if val == "tand_":
        string = r"$\tan \delta = %.4f$"
        endmark = ".txt"
    if val == "UCDim_":
        string = r"$L^\mathrm{UC} = %.2f$"
        endmark = "_lz" 
    result = []
    for filE in file_list:
        daten = np.loadtxt(filE, delimiter=",")
        f, S11 = daten[:,0], np.abs(daten[:,1] + 1j*daten[:,2])
        dataset = np.append(f[:,np.newaxis], S11[:,np.newaxis], axis=1)
        filename = filE[filE.find("S11_f_UCDim_"):]
        idx0 = filename.find(val) + len(val)
        s = filename[idx0:]
        print("working on file: ", s)
        idx1 = s.find(endmark)
        s = float(s[:idx1])
        print("%s=" %val, s)
        # make it a sorted list of tuples
        # with first entry equal to the float
        # value of "val" i. e. eps, tand ...
        # and the second entry a numpy array
        # containing frequencies vs. |S11|
        # in the first/second column
        result.append((s, dataset))
    result.sort(key = lambda x: x[0]) 
    return result


files = fileList("/home/stefan/Arbeit/latex/DoubleRings/code/double_ring_eps_sweep/tand_sweep_eps_4.1/", "S11_f_UCDim_20_lz_2", ".txt")
plot_data = get_plot_data(files, "tand_")
col = ["r", "b", "g", "m", "c"]
fig = plt.figure()
ax = fig.add_subplot(111)

for index in range(2): # increase index to study more resonances
    counter = 0
    for data in plot_data:
        dset = data[1]
        val = data[0] # the parameter which is swept
        mask = get_minima_positions(dset[:,1]**2, 0.5, N=1)
        S11 = 20*np.log10(dset[mask,1])[index]
        f = dset[mask,0][index]
        if counter == 0:
            norm = S11 
            ax.plot(val, S11/norm, label="$f=%.2f$ GHz" %(f/1e9),
                    marker="o", color=col[index])
        else:
            ax.plot(val, S11/norm, marker="o", color=col[index])
        counter += 1

ax.set_title(r"Doppelringabsorber, Einfluss von $\tan\delta^\mathrm{FR4}$ auf $|S_{11}|$")
ax.set_ylabel(r"$\log S_{11}(\tan\delta)/\log S_{11}(\tan\delta=0.005)$", fontsize=16)
#plt.xlim([3.99, 4.651])
#ax.plot([1, 2], [1, 1], "k--")
ax.set_xlabel(r"$\tan\delta^\mathrm{FR4}$", fontsize=16)
ax.tick_params(axis="both", labelsize=14)
#ax.set_ylim([0.99, 1.2])
#ax.set_xlim([1, 2.01])
ax.legend(loc="best").draw_frame(False)
fig.savefig("Einfluss_tand_absS11.pdf", format="pdf")
#plt.show()


if __name__ == "__main__":
    pass
  
