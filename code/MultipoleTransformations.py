import numpy as np
import spherical_functions as sf
import quaternion
from cython_gaunt import gaunt
from scipy.special import hankel1, jv
from math import sqrt, pi
from scipy.linalg import block_diag


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

def full_rotation_matrix(lmax, theta_y, phi_z):
    """ Get the full Rotation Matrix
    """
    rotor = quaternion.from_spherical_coords(theta_y, phi_z)
    blocklist = WignerD_matrices(lmax, theta_y, phi_z)
    return block_diag(*blocklist)

def translate_l1l2(l1, l2, m, kd, sign_z, regreg):
    """return the (l1, m) contribution of a vector spherical 
       harmonic M_l1,m (N_l1,m) when it is translated along (sign_z = +1)
       or against (sign_z = -1) the z-direction an adimensional distance
       kd. 
       Returns the polarization conserving contribution M_l1,m (N_l1,m) 
       as the first and the polarization mixing contribution N_l1,m (M_l1,m)
       as the second value
    """
    if abs(m) > max(l1, l2):
        raise ValueError("translate_l1l2 should not be called \
        with |m| > max(l1, l2), but got l1=%i, l2=%i, m=%i" %(l1,l2,m))
    common =  ((-1)**m * (sign_z*1j)**(l1-l2)*2*   
               sqrt(pi/(l1*(l1+1)*l2*(l2+1))))
    alpha = np.arange(abs(l1-l2), l1+l2+1, 2)
    gaunts = gaunt(l1, l2, m)
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
    for l1 in range(lmin, l1_max+1):
        for l2 in range(l1, l2_max+1):
            # fill elements in the "upper" triangle
            i, j = l1-lmin, l2-lmin
            PPij, PQij = translate_l1l2(l1, l2, m, kd, sign_z, regreg)
            result[i, j] = PPij
            result[i, j+dim2] = PQij
        for l2 in range(lmin+1, l1_max+1):
            transpose_factor = (-1)**(l1-l2)
            j = l2-lmin
            # use V(l1,l2) = (-1)^(l1-l2)* V(l2,l1)
            result[j, i] = result[i, j] * transpose_factor
            result[j, i+dim2] = result[i, j+dim2] * transpose_factor
            # now the first line in the matrix [[PP, PQ],
            # is filled and ready to copy       [PQ, PP]]
    # copy the PP and PQ of the first line to the second
    # line in the correct order
    result[dim1:, dim2:] = result[0:dim1, 0:dim2] # PP to lower right corner
    result[dim1:, 0:dim2] = result[0:dim1, dim2:] # PQ to lower left corner
    if transpose:
        result = result.transpose()
    return result
 
def full_translation_matrix(l1max, l2max, kd, sign_z, regreg):
    """ get the full translation matrix
    """
    valid_ms = np.arange(0, min(l1max, l2max) + 1)
    Tmatrices = [translation_matrix(l1max, l2max, m, kd, sign_z, regreg) for
                 m in valid_ms]
    Tneg = Tmatrices[1:].copy()
    Tneg.reverse()
    Tneg.extend(Tmatrices)
    return block_diag(*Tneg)

def basisChange_l_to_m(lmax):
    """ return the basis-change matrix from l-ordered
        to m-ordered
    """
    valid_ms = np.arange(-lmax, lmax+1)
    dim = 2*lmax*(lmax+2)
    dh = lmax*(lmax+2)
    Bml = np.zeros((dim, dim))
    i = 0
    for m in valid_ms:
        lmin = max(abs(m), 1)
        for l in range(lmin, lmax+1):
            zero = np.zeros(dh)
            zero[l*(l+1) + m -1] =1
            Bml[i, 0:dh] = zero
            Bml[i+dh, dh:] = zero
            i += 1
    return Bml



def translation_matrix_debug(l1_max, l2_max, m, kd, sign_z, regreg): 
    """Calculate the matrix without using symmetries
       calculate the full translation matrix with dimensions:
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
    lmin = max(abs(m), 1)
    dim1 = l1_max - lmin + 1
    dim2 = l2_max - lmin + 1
    # polarization conserving block
    result = np.zeros((2*dim1, 2*dim2), dtype = np.complex128)
    # polarization mixing block
    # polarization offset in lines (columns) i (j)
    Di = l1_max - lmin + 1
    Dj = l2_max - lmin + 1
    for l1 in range(lmin, l1_max+1):
        for l2 in range(lmin, l2_max+1):
            # fill elements in the "upper" triangle
            i, j = l1-lmin, l2-lmin
            PPij, PQij = translate_l1l2(l1, l2, m, kd, sign_z, regreg)
            result[i, j] = PPij
            result[i, j+Dj] = PQij
    result[Di:, Dj:] = result[0:Di, 0:Dj] # PP to lower right corner
    result[Di:, 0:Dj] = result[0:Di, Dj:] # PQ to lower left corner
    return result
 

if __name__ == "__main__":
    #print(WignerD_matrices(2, 0, 1.5))
    l1_max, l2_max = 4, 10 
    m, kd = 1, 0.5
    sign_z, regreg = 1, 1
    T1 = translation_matrix_debug(l1_max, l2_max, m, kd, +sign_z, regreg)
    T2 = translation_matrix_debug(l2_max, l1_max, m, kd, -sign_z, regreg)
    print(np.dot(T1, T2))
