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
from smooth import smooth



def get_min_mask(dataset, below, N=21, P=0.5, Sig=25):
    """ get the position of the local minima
        in the 1d numpy array >dataset< after
        smoothing it with a gaussian
        Returns: indices of the local minima,
        where the value of dataset[indices] is 
        smaller than >above<
    """
    filtered = smooth(dataset,31)[15:-15]
    min_pos = argrelmin(filtered)[0][1:-1]
    mask = min_pos[filtered[min_pos] < below]
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
        result.append((s, dataset))
    result.sort(key = lambda x: x[0])
    return result
