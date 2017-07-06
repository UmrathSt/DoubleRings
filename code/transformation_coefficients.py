import numpy as np
import spherical_functions as sf
import quaternion
from cython_gaunt import gaunt
from scipy.special import hankel1, jv
from math import sqrt, pi


def WignerD_matrices(l, theta_y, phi_z):
    rotor = quaternion.from_spherical_coords(theta_y, phi_z)
    D = sf.WignerD.Wigner_D_matrices(rotor, 1, l)
    D1 = D[0:9]
    result = [np.array([D1[3*i:(i+1)*3] for i in range(3)], dtype = np.complex128)]
    for alpha in range(2, l + 1):
        start = int((4*alpha**3 - alpha - 3)/3)
        line = 2*alpha + 1
        Dl = D[start:start+line**2]
        for j in range(1, line):
            Matrix = np.array([Dl[j*line:(j+1)*line] for j in range(line)], 
                                dtype = np.complex128)
        result.append(Matrix)
    return result

def translate_l1l2(l1, l2, m, kd, sign_z, regreg):
    """return the (l1, m) contribution of a vector spherical 
       harmonic M_l1,m (N_l1,m) when it is translated along (sign_z = +1)
       or against (sign_z = -1) the z-direction an adimensional distance
       kd. 
       Returns the polarization conserving contribution M_l1,m (N_l1,m) 
       as the first and the polarization mixing contribution N_l1,m (M_l1,m)
       as the second value
    """
    # common has extra "i" since M/N are defined 
    # differently compared to Wittann
    if abs(m) > max(l1, l2):
        raise ValueError("translate_l1l2 should not be called \
        with |m| > max(l1, l2), but got l1=%i, l2=%i, m=%i" %(l1,l2,m))
    common =  ((-1)**m * (sign_z*1j)**(l1-l2)*2*   
               sqrt(pi/(l1*(l1+1)*l2*(l2+1))))
    alpha = np.arange(abs(l1-l2), l1+l2+1, 2)
    gaunts = gaunt(l1, l2, m)
    spherical_factor = sqrt(pi/(2*kd))
    if regreg:
        bessel = jv(alpha+0.5, kd)
    else: 
        bessel = hankel1(alpha+0.5, kd)
    common_alpha_factor = (np.sqrt((2*alpha+1)*pi/(2*kd)) * 
                    1j**alpha * bessel * gaunts)
    fak_PP = ((l1*(l1+1) + l2*(l2+1) - alpha*(alpha+1)) * 0.5 
               * common_alpha_factor).sum()
    fak_PQ = (-m*1j*kd * common_alpha_factor).sum()
    PP = fak_PP * common
    PQ = fak_PQ * common
    return [PP, PQ]

def translation_matrix(l1_max, l2_max, m, kd, sign_z, regreg): 
    """calculate the full translation matrix with dimensions:
       [2 * (l1_max - max(|m|, 1)) x 2 * (l2_max - max(|m|, 1))]
       wavenumber*distance = kd and translation in positive
       (sign_z = +1) or negative (sign_z = -1) z-direction
       regreg should be either +1 for regular-wave to regular-wave
       or -1 for regular -> outgoing waves
       Returns a rectangular Matrix of block-structure:
       result = [[TETE, TETM],
                 [TMTE, TMTM]]
       since the result does not depend upon the sign of m 
       it should be not recalculated if needed for m = -lmax..+lmax
    """
    transpose = False
    if l1_max > l2_max:
        transpose = True
        l1_max, l2_max = l2_max, l1_max
    smaller_l = min(l1_max, l2_max)
    lmin = max(abs(m), 1)
    dim1 = l1_max - lmin + 1
    dim2 = l2_max - lmin + 1
    # polarization conserving block
    result = np.zeros((2*dim1, 2*dim2), dtype = np.complex128)
    # polarization mixing block
    for l1 in range(lmin, l1_max+1):
        for l2 in range(l1, l2_max+1):
            PP, PQ = translate_l1l2(l1, l2, m, kd, sign_z, regreg)
            # polarization offset in lines (columns) i (j)
            Di = l1_max - lmin + 1
            Dj = l2_max - lmin + 1
            result[l1-lmin, l2-lmin] = PP
            result[l1-lmin+Di, l2-lmin+Dj] = PP
            result[l1-lmin, l2-lmin+Dj] = PQ
            result[l1-lmin+Di, l2-lmin] = PQ
            l1l2_switch = -1j*sign_z
            if ([l2-lmin+1+Dj, l1-lmin+1+Di] <= [result.shape[0]/2,
                        result.shape[1]/2] and not l1 == l2):
                result[l2-lmin, l1-lmin] = l1l2_switch*PP
                result[l2-lmin+Dj, l1-lmin+Di] = l1l2_switch*PP
                result[l2-lmin+Dj, l1-lmin] = l1l2_switch*PQ 
                result[l2-lmin, l1-lmin+Di] = l1l2_switch*PQ
    if transpose:
        result.transpose()
    return result
            


if __name__ == "__main__":
    #print(WignerD_matrices(2, 0, 1.5))
    
    T = translation_matrix(5, 3, 0.1, 1, 1)
    print(T[0][0])
    print(T[1][0])
    print(T[0][-2])
    print(T[1][-2])

