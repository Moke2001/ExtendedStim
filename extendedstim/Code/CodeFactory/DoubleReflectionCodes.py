import numpy as np
import galois
from extendedstim.Code.PrimitiveCode.MajoranaCSSCode import MajoranaCSSCode
from extendedstim.Code.PrimitiveCode.PauliCSSCode import PauliCSSCode
from extendedstim.Physics.MajoranaOperator import MajoranaOperator
from extendedstim.Physics.PauliOperator import PauliOperator
from extendedstim.tools.GaloisTools import random_distance_caculator


##  CHAPTER: ====Details====
class DoubleReflectionCode:
    """
    self.l：Dimension of the loop matrix A
    self.m：Dimension of the loop matrix B
    self.H：Check matrix
    """

    __slots__=['l','m','x_list_a','y_list_a','x_list_b','y_list_b','p_list_x_a','p_list_y_a','p_list_x_b','p_list_y_b','H']

    def __init__(self,l:int,m:int,x_list_a:list[int],y_list_a:list[int],x_list_b:list[int],y_list_b:list[int],p_list_x_a:list[int],p_list_y_a:list[int],p_list_x_b:list[int],p_list_y_b:list[int]) -> None:
        """""
        input.l：循环矩阵的大小
        input.m：循环矩阵的大小
        input.x_list_a：A中x项的移位幂次
        input.y_list_a：A中y项的移位幂次
        input.x_list_b：B中x项的移位幂次
        input.y_list_b：B中y项的移位幂次
        input.r_list_a：A中对应项的P矩阵幂次，可以是[r_x_list, r_y_list]或平铺列表
        input.r_list_b：B中对应项的P矩阵幂次，可以是[r_x_list, r_y_list]或平铺列表
        output：无
        influence：构造校验矩阵self.H
        """""
        
        # PART: Details
        self.l=l
        self.m=m
        self.x_list_a=x_list_a
        self.y_list_a=y_list_a
        self.x_list_b=x_list_b
        self.y_list_b=y_list_b
        self.p_list_x_a=p_list_x_a
        self.p_list_y_a=p_list_y_a
        self.p_list_x_b=p_list_x_b
        self.p_list_y_b=p_list_y_b

        # PART: Details
        A = self._construct_matrix(l, m, x_list_a, y_list_a, p_list_x_a, p_list_y_a)
        B = self._construct_matrix(l, m, x_list_b, y_list_b, p_list_x_b, p_list_y_b)
        
        # PART: Details
        if self._check_conditions(A, B):
            A_T = A.T
            B_T = B.T
            H_top = np.hstack((A, B, A_T, B_T))
            H_bottom = np.hstack((B_T, A_T, B, A))
            self.H = galois.GF(2)(np.vstack((H_top, H_bottom)))
        else:
            A_T = A.T
            B_T = B.T
            H_top = np.hstack((A, B, A_T, B_T))
            H_bottom = np.hstack((B_T, A_T, B, A))
            self.H = galois.GF(2)(np.vstack((H_top, H_bottom)))
            #self.H = None

    # SECTION: Details
    def _parse_r_list(self, r_list, len_x, len_y):
        if len(r_list) == 2 and isinstance(r_list[0], (list, tuple, np.ndarray)):
             return r_list[0], r_list[1]
        elif len(r_list) == len_x + len_y:
             return r_list[:len_x], r_list[len_x:]
        else:
             raise ValueError(f"r_list length {len(r_list)} does not match expected length {len_x}+{len_y} or structure [list, list]")

    # SECTION: Details
    def _construct_matrix(self, l, m, x_list, y_list, r_x, r_y):
        GF = galois.GF(2)
        # Term 1: sum x_i S_l^i P_l^a_ix \otimes I_m
        term1 = GF(np.zeros((l*m, l*m), dtype=int))
        for i, r in zip(x_list, r_x):
            S_l = shift(l, i)
            P_l = matrix_power_P(l, r)
            op = GF(np.kron(S_l @ P_l, np.eye(m, dtype=int)))
            term1 += op
            
        # Term 2: sum y_i I_l \otimes S_m^i P_m^a_iy
        term2 = GF(np.zeros((l*m, l*m), dtype=int))
        for i, r in zip(y_list, r_y):
            S_m = shift(m, i)
            P_m = matrix_power_P(m, r)
            op = GF(np.kron(np.eye(l, dtype=int), S_m @ P_m))
            term2 += op
            
        return term1 + term2

    # SECTION: Details
    def _check_conditions(self, A, B):
        # 1. A A^T = A^T A
        cond1 = np.array_equal(A @ A.T, A.T @ A)
        # 2. B B^T = B^T B
        cond2 = np.array_equal(B @ B.T, B.T @ B)
        # 3. A B = B A
        cond3 = np.array_equal(A @ B, B @ A)
            
        return cond1 and cond2 and cond3
        
    # SECTION: Details
    def __str__(self):
        if self.H is None:
            return "DoubleUnamedCode(Construction Failed)"
        A_str='A='
        B_str='B='
        for i in range(len(self.x_list_a)):
            if i==0:
                A_str+=f'x^{self.x_list_a[i]}p^{self.p_list_x_a[i]}y^{self.y_list_a[i]}q^{self.p_list_y_a[i]}'
                B_str+=f'x^{self.x_list_b[i]}p^{self.p_list_x_b[i]}y^{self.y_list_b[i]}q^{self.p_list_y_b[i]}'
            else:
                A_str+=f'+x^{self.x_list_a[i]}p^{self.p_list_x_a[i]}y^{self.y_list_a[i]}q^{self.p_list_y_a[i]}'
                B_str+=f'+x^{self.x_list_b[i]}p^{self.p_list_x_b[i]}y^{self.y_list_b[i]}q^{self.p_list_y_b[i]}'
        return f"{self.l}|{self.m}|"+A_str+'|'+B_str


# CHAPTER: Details
class MajoranaDoubleReflectionCode(MajoranaCSSCode, DoubleReflectionCode):
    
    ##  SECTION: ----Constructor----
    def __init__(self,l:int,m:int,x_list_a:list[int],y_list_a:list[int],x_list_b:list[int],y_list_b:list[int],p_list_x_a:list[int],p_list_y_a:list[int],p_list_x_b:list[int],p_list_y_b:list[int]) -> None:
        DoubleReflectionCode.__init__(self, l, m, x_list_a, y_list_a, x_list_b, y_list_b, p_list_x_a, p_list_y_a, p_list_x_b, p_list_y_b)
        if self.H is None:
            raise ValueError("Code construction failed due to unsatisfied conditions.")
            
        generators_x=[]
        generators_z=[]
        for i in range(len(self.H)):
            generators_x.append(MajoranaOperator.HermitianOperatorFromOccupy(np.where(self.H[i]!=0)[0],[]))
            generators_z.append(MajoranaOperator.HermitianOperatorFromOccupy([],np.where(self.H[i]!=0)[0]))
        MajoranaCSSCode.__init__(self,generators_x,generators_z,self.l*self.m*4)


##  CHAPTER: ====Details====
class PauliDoubleReflectionCode(PauliCSSCode, DoubleReflectionCode):
    
    ##  SECTION: ----Constructor----
    def __init__(self,l:int,m:int,x_list_a:list[int],y_list_a:list[int],x_list_b:list[int],y_list_b:list[int],p_list_x_a:list[int],p_list_y_a:list[int],p_list_x_b:list[int],p_list_y_b:list[int]) -> None:
        DoubleReflectionCode.__init__(self, l, m, x_list_a, y_list_a, x_list_b, y_list_b, p_list_x_a, p_list_y_a, p_list_x_b, p_list_y_b)
        if self.H is None:
            raise ValueError("Code construction failed due to unsatisfied conditions.")
            
        generators_x=[]
        generators_z=[]
        for i in range(len(self.H)):
            generators_x.append(PauliOperator.HermitianOperatorFromOccupy(np.where(self.H[i]!=0)[0],[]))
            generators_z.append(PauliOperator.HermitianOperatorFromOccupy([],np.where(self.H[i]!=0)[0]))
        PauliCSSCode.__init__(self,generators_x,generators_z,self.l*self.m*4)


##  CHAPTER: ====Shift====
def shift(number:int, shift_amount:int):
    """""
    input.number：矩阵维度
    input.shift_amount：循环移位步长
    output：GF(2)矩阵，置换矩阵
    """""
    S = np.zeros((number, number), dtype=int)
    for i in range(number):
        S[i, (i + shift_amount) % number] = 1
    return galois.GF(2)(S)


##  SECTION: ----Power----
def matrix_power_P(number:int, r:int):
    """""
    input.number：矩阵维度
    input.r：幂次
    output：GF(2)矩阵，P^r
    """""
    if r % 2 == 0:
        return galois.GF(2)(np.eye(number, dtype=int))
    else:
        # P is anti-diagonal
        return galois.GF(2)(np.fliplr(np.eye(number, dtype=int)))

if __name__ == "__main__":
    code=DoubleReflectionCode(4,4,[1,0,2],[3,3,2],[3,3,2],[1,0,1],[0,1,1],[1,1,1],[1,0,1],[0,0,1])
    #print(np.all(code.H@code.H.T==0))
    """
x_a=[1, 0, 2], y_a=[3, 3, 2], px_a=[0, 1, 1], py_a=[1, 1, 1]
x_b=[3, 3, 2], y_b=[1, 0, 1], px_b=[1, 0, 1], py_b=[0, 0, 1]
    """
    print(random_distance_caculator(code.H,[],2000))
