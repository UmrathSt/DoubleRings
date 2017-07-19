# coding: utf-8
# get the dependence of the depth and position of the
# absorption frequencies as a function of the real and imaginary part
# of the substrate permittivity for a double-ring geometry

import numpy as np
from matplotlib import pyplot as plt
from filewalk import fileList
from scipy.signal import argrelmin, general_gaussian, fftconvolve

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
    mask = min_pos[filtered[min_pos] < 0.9]
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
        result.append((s, dataset))
    result.sort(key = lambda x: x[0])
    return result


files = fileList("/home/stefan/Arbeit/latex/DoubleRings/code/double_ring_eps_sweep/tand_sweep/", "S11_f_UCDim", ".txt")
plot_data = get_plot_data(files, "tand_")
col = ["r", "b", "g", "m", "c"]
for index in range(3):
    counter = 0
    for data in plot_data:
        dset = data[1]
        eps = data[0]
        mask = get_minima_positions(dset[:,1], 0.9)
        f = dset[mask,0][index]
        if counter == 0:
            normalization = f
            plt.plot(eps, f/normalization, label="$f=%.2f$ GHz" %(f/1e9), marker="o", color=col[index])
        else:
            plt.plot(eps, f/normalization, marker="o", color=col[index])
        counter += 1
plt.legend(loc="best").draw_frame(False)
plt.show()


if __name__ == "__main__":
    pass
  
