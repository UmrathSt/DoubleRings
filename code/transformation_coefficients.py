import numpy as np
import spherical_functions as sf
import quaternion
from cython_gaunt import gaunt
from scipy.special import hankel1, jv


def WignerD_matrices(l, theta_y, phi_z):
    rotor = quaternion.from_spherical_coords(theta_y, phi_z)
    D = sf.WignerD.Wigner_D_matrices(rotor, 1, l)
    D1 = D[0:9]
    print(D1)
    result = [np.array([D1[3*i:(i+1)*3] for i in range(3)], dtype = np.complex128)]
    for alpha in range(2, l + 1):
        start = int((4*alpha**3 - alpha - 3)/3)
        line = 2*alpha + 1
        Dl = D[start:start+line**2]
        print("shape D_l=%i: " %alpha, np.shape(Dl))
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
    common =  ((-1)**m*1j**(l1-l2)/   
               np.sqrt(l2*(l2+1) * l1*(l1+1)/(np.pi*4)))
    alpha = np.arange(abs(l1-l2), l1+l2+1, 2)
    gaunts = gaunt(l1, l2, m)
    if regreg:
        sphB = jv(alpha+0.5, kd)*np.sqrt(np.pi/(2*kd))
    else: 
        sphB = hankel1(alpha+0.5, kd)*np.sqrt(np.pi/(2*kd))
    fak_PP = ((l1*(l1+1) + l2*(l2+1) - alpha*(alpha+1))*0.5 * sphB *
                 (sign_z*1j)**alpha*gaunts*np.sqrt(2*alpha+1)).sum()
    fak_PQ = (-m*kd*(1j)**(alpha+1) * gaunts * np.sqrt(2*alpha+1)
              *sphB * sign_z**(alpha+1)).sum()
    PP = fak_PP * common
    PQ = fak_PQ * common
    return [PP, PQ]

def translation_matrix(l1_max, l2_max, kd, sign_z, regreg): 
    """calculate the full translation matrix assuming that lmax is sufficient
       in both coordinate systems and return a list:
       [[PP(m=-lmax), PP(m=-lmax+1..., PP(m=0)], [PQ(m=-lmax),...,PQ(m=0)]] 
       square matrices PP and PQ
       for m = -l, -l+1, ..., 0 which is sufficient since translations in z-direction
       do not depend upon the sign of m and don't mix in m
    """
    PP_result = []
    PQ_result = []
    smaller_l = min(l1_max, l2_max)
    for m in range(-smaller_l, 1):
        lmin = max(abs(m), 1)
        dim1 = l1_max - lmin + 1
        dim2 = l2_max - lmin + 1
        PPm = np.zeros((dim1, dim2), dtype = np.complex128)
        PQm = np.zeros(PPm.shape, dtype = np.complex128)
        for l1 in range(lmin, l1_max+1):
            for l2 in range(l1, l2_max+1):
                PP, PQ = translate_l1l2(l1, l2, m, kd, 
                                        sign_z, regreg)
                PPm[l1-lmin, l2-lmin] = PP
                PQm[l1-lmin, l2-lmin] = PQ * (-1j) 
                if [l2-lmin+1, l1-lmin+1] <= list(PPm.shape):
                    PPm[l2-lmin, l1-lmin] = -1j*PP
                    PQm[l2-lmin, l1-lmin] = -1j*PQ * (-1j)
        PP_result.append(PPm)
        PQ_result.append(PQm)
    return [PP_result, PQ_result]
            


if __name__ == "__main__":
    #print(WignerD_matrices(2, 0, 1.5))
    
    T = translation_matrix(5, 3, 0.1, 1, 1)
    print(T[0][0])
    print(T[1][0])
    print(T[0][-2])
    print(T[1][-2])

