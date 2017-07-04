import numpy as np
import spherical_functions as sf
import quaternion
from cython_gaunt import gaunt
from scipy.special import hankel1

h1 = lambda l, x: hankel1(l+0.5, x)*np.sqrt(np.pi/(2*x))

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

def translate_l1l2(l1, l2, m, kd, sign_z):
    """return the (l1, m) contribution of a vector spherical 
       harmonic M_l1,m (N_l1,m) when it is translated along (sign_z = +1)
       or against (sign_z = -1) the z-direction an adimensional distance
       kd. 
       Returns the polarization conserving contribution M_l1,m (N_l1,m) 
       as the first and the polarization mixing contribution N_l1,m (M_l1,m)
       as the second value
    """
    common =  ((-1)**m*4*np.pi * 1j**(l2-l1)/
               np.sqrt(l2*(l2+1) * l1*(l1+1)))
    alpha = np.arange(abs(l1-l2), l1+l2+1, 2)
    gaunts = gaunt(l1, l2, m)
    h1_sph= h1(alpha+0.5, kd)
    fak_PP = ((l1*(l1+1) + l2*(l2+1) - alpha*(alpha+1))*0.5 * h1_sph *
                 (sign_z*1j)**alpha*gaunts*np.sqrt((2*alpha+1)/(4*np.pi))).sum()
    fak_PQ = (-m*kd*(1j)**(alpha+1) * gaunts * np.sqrt((2*alpha+1)/(4*np.pi))
              *h1_sph * sign_z**(alpha+1)).sum()
    PP = fak_PP * common
    PQ = fak_PQ * common
    return [PP, PQ]

def translation_matrix(lmax, kd, sign_z):
    """calculate the full translation matrix assuming that lmax is sufficient
       in both coordinate systems and return a list:
       [[PP(m=-lmax), PP(m=-lmax+1..., PP(m=0)], [PQ(m=-lmax),...,PQ(m=0)]] 
       square matrices PP and PQ
       for m = -l, -l+1, ..., 0 which is sufficient since translations in z-direction
       do not depend upon the sign of m and don't mix in m
    """
    PP_result = []
    PQ_result = []
    for m in range(-lmax, 1):
        dim = max(lmax - abs(m) + 1, 1)
        if abs(m) <= 1:
            dim = lmax
        PPm = np.zeros((dim, dim), dtype = np.complex128)
        PQm = np.zeros(PPm.shape, dtype = np.complex128)
        lmin = max(abs(m), 1)
        for l1 in range(lmin, lmax+1):
            for l2 in range(l1, lmax+1):
                PP, PQ = translate_l1l2(l1, l2, m, kd, sign_z)
                PPm[l1-lmin, l2-lmin] = PP
                PPm[l2-lmin, l1-lmin] = -1j*PP
                PQm[l1-lmin, l2-lmin] = PQ
                PQm[l2-lmin, l1-lmin] = -1j*PQ
        PP_result.append(PPm)
        PQ_result.append(PQm)
    return [PP_result, PQ_result]
            


if __name__ == "__main__":
    #print(WignerD_matrices(2, 0, 1.5))
    print(translation_matrix(2, 1, 1))

