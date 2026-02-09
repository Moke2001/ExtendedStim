import random

import galois
import numpy as np

from extendedstim.tools.GaloisTools import random_distance_caculator, random_distance_caculator_css


class ConstructiveReflectionCode:
    __slots__=['l','m','A','B','H_X','H_Z','H']

    def __init__(self, l:int, m:int, y_list_a:list[int], p_list_y_a:list[int], x_list_b:list[int], p_list_x_b:list[int]) -> None:
        self.l=l
        self.m=m
        GF = galois.GF(2)
        term_y = GF(np.zeros((l*m, l*m), dtype=int))
        for i, r in zip(y_list_a, p_list_y_a):
            S_m = shift(m, i)
            if np.mod(r,2)==0:
                P_m=np.eye(m, dtype=int)
            else:
                P_m=np.zeros((m,m), dtype=int)
                for j in range(m):
                    P_m[j,m-j-1]=1
            op = GF(np.kron(np.eye(l, dtype=int), S_m @ P_m))
            term_y += op
        term_x = GF(np.zeros((l*m, l*m), dtype=int))
        for i, r in zip(x_list_b, p_list_x_b):
            S_l = shift(l, i)
            if np.mod(r, 2) == 0:
                P_l = np.eye(l, dtype=int)
            else:
                P_l = np.zeros((l, l), dtype=int)
                for j in range(l):
                    P_l[j, l - j - 1] = 1
            op = GF(np.kron(S_l @ P_l, np.eye(m, dtype=int)))
            term_x += op
        self.A = term_y
        self.B = term_x
        H_X = np.hstack((self.A, self.B))
        H_Z = np.hstack((self.B.T, self.A.T))
        self.H_X = GF(H_X)
        self.H_Z = GF(H_Z)
        self.H = GF(np.vstack((H_X, H_Z)))


def shift(n, k):
    """
    Returns the shift matrix of size n x n with shift k.
    S_ij = 1 if j = i + k mod n
    """
    S = np.zeros((n, n), dtype=int)
    for i in range(n):
        S[i, (i + k) % n] = 1
    return S


if __name__ == '__main__':
    l=4
    m=5
    for _ in range(2000):
        y_list_a=np.random.choice(range(l),2)
        p_list_y_a=np.random.choice([0,1],2)
        x_list_b=np.random.choice(range(m),2)
        p_list_x_b=np.random.choice([0,1],2)
        code=ConstructiveReflectionCode(l,m,y_list_a,p_list_y_a,x_list_b,p_list_x_b)
        assert np.all(code.H_X@code.H_Z.T==0)
        physical_number=l*m*2
        logical_number=code.H_X.shape[1]-np.linalg.matrix_rank(code.H_X)-np.linalg.matrix_rank(code.H_Z)
        code_distance=random_distance_caculator_css(code.H_X, code.H_Z,1000)
        print("physical_number:",physical_number,"logical_number:",logical_number,"code_distance:",code_distance)
