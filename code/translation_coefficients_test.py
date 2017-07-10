import numpy as np
from MultipoleTransformations import translate_l1l2 as trans
from MultipoleTransformations import translation_matrix
from math import isclose

def test_symmetry(lmax, kd):
    passed = True
    for sign_z in [-1, 1]:
        for regreg in [1, 0]:
            for m in range(-lmax, lmax+1):
                for l1 in range(max(abs(m), 1), lmax+1):
                    for l2 in range(l1, lmax+1):
                        PPl1l2, PQl1l2 = trans(l1, l2, m, kd, sign_z, regreg)
                        PPl2l1, PQl2l1 = trans(l2, l1, m, kd, sign_z, regreg)
                        PPl2l1 *= (-1)**(l1-l2)
                        PQl2l1 *= (-1)**(l1-l2)
                        if not (isclose(PPl1l2.real, PPl2l1.real)):
                            print("deviation in real parts of PP")
                            print(PPl1l2.real, "vs. ", PPl2l1.real)
                            passed = False
                        if not (isclose(PPl1l2.imag, PPl2l1.imag)):
                            print("deviation in imag parts of PP")
                            print(PPl1l2.imag, "vs. ", PPl2l1.imag)
                            passed = False
                        if not (isclose(PQl1l2.real, PQl2l1.real)):
                            passed = False
                            print("deviation in real parts of PQ")
                            print(PQl1l2.real, "vs. ", PQl2l1.real)
                        if not (isclose(PQl1l2.imag, PQl2l1.imag)):
                            passed = False
                            print("deviation in imag parts of PQ")
                            print(PQl1l2.imag, "vs. ", PQl2l1.imag)
    return passed

def test_inverse(A, B):
    """ Test if A and B are inverse matrices
        If matrices are rectangular use a
        pseudo-inverse
    """
    passed = True
    i, j = A.shape
    dim = i
    if not i == j:
        dim = max(i, j)
        Asq = np.eye(dim)
        Asq[0:i, 0:j] = A
        A = Asq
        Bsq = np.eye(dim)
        Bsq[0:i, 0:j] = B
        B = Bsq
    AB = np.dot(A, B)
    for i in range(dim):
        for j in range(dim):
            elem = AB[i,j]
            if not (abs(elem.imag) < 1e-2):
                print("nonvanishing imaginary part at (i,j)=(%i,%i), %.2f" %(i,j, elem.imag))
                passed = False
                break
            if i == j and not (abs(elem.real-1) < 1e-2):
                passed = False
                print("non-1 real part at (i,j)=(%i,%i), %.2f" %(i, j, elem.real))
                break
    return passed



if __name__ == "__main__":
    A1 = translation_matrix(10, 10, 1, 1, 1, 1)
    A2 = translation_matrix(10, 10, 1, 1, 1, 1)
    A = np.dot(A1, A2)
    B = translation_matrix(10, 10, 1, 2, -1, 1) 
    tests = [test_symmetry(10, 3), test_inverse(A, B)]
    for i, test in enumerate(tests):
        print("Test %i passed: " %(i+1), test)
