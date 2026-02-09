import numpy as np
from extendedstim.Code.PrimitiveCode.MajoranaCSSCode import MajoranaCSSCode
from extendedstim.Code.PrimitiveCode.PauliCSSCode import PauliCSSCode
from extendedstim.Physics.MajoranaOperator import MajoranaOperator
from extendedstim.Physics.PauliOperator import PauliOperator
import galois


##  CHAPTER: ====Twisted Double Bivariant Bicycle Code====
class TwistedDoubleBivariantBicycleCode:
    __slots__=['gamma','l','m','a_x_list','a_y_list','b_x_list','b_y_list','H']

    ##  SECTION: ----Constructor----
    def __init__(self,l:int,m:int,gamma:int,a_x_list:list[int],a_y_list:list[int],b_x_list:list[int],b_y_list:list[int]):
        """""
        input.l：even dimension of the code, i.e. the number of qubits in each cycle.
        input.m：even dimension of the code, i.e. the number of qubits in each cycle.
        input.a_x_list：list of exponents of x in A, i.e. the power of x in A.
        input.a_y_list：list of exponents of y in A, i.e. the power of y in A.
        input.b_x_list：list of exponents of x in B, i.e. the power of x in B.
        input.b_y_list：list of exponents of y in B, i.e. the power of y in B.
        """""

        self.l=l
        self.m=m
        self.gamma=np.mod(gamma,m)
        self.a_x_list=a_x_list
        self.a_y_list=a_y_list
        self.b_x_list=b_x_list
        self.b_y_list=b_y_list
        GF=galois.GF(2)
        x=GF(np.kron(shift(l, 1), np.eye(m, dtype=int)))
        x[m*(l-1):,0:m]=shift(m,self.gamma)
        y=GF(np.kron(np.eye(l, dtype=int), shift(m, 1)))
        A=None
        B=None
        for i in range(len(a_x_list)):
            if A is None:
                A=np.linalg.matrix_power(x, a_x_list[i])@np.linalg.matrix_power(y, a_y_list[i])
            else:
                A=A+np.linalg.matrix_power(x, a_x_list[i])@np.linalg.matrix_power(y, a_y_list[i])
        for i in range(len(b_x_list)):
            if B is None:
                B=np.linalg.matrix_power(x, b_x_list[i])@np.linalg.matrix_power(y, b_y_list[i])
            else:
                B=B+np.linalg.matrix_power(x, b_x_list[i])@np.linalg.matrix_power(y, b_y_list[i])
        A_T=A.T
        B_T=B.T
        H_top=np.hstack((A, B, A_T, B_T))
        H_bottom=np.hstack((B_T, A_T, B, A))
        H=GF(np.vstack((H_top, H_bottom)))
        self.H=H


    ##  SECTION: ----String representation----
    def __str__(self):
        A_str=''
        B_str=''
        for i in range(len(self.a_x_list)):
            if i==0:
                A_str=A_str+f"x^{self.a_x_list[i]}y^{self.a_y_list[i]}"
            else:
                A_str=A_str+f"+x^{self.a_x_list[i]}y^{self.a_y_list[i]}"
        for i in range(len(self.b_x_list)):
            if i==0:
                B_str=B_str+f"x^{self.b_x_list[i]}y^{self.b_y_list[i]}"
            else:
                B_str=B_str+f"+x^{self.b_x_list[i]}y^{self.b_y_list[i]}"
        return f"|{self.l*self.m*4}|{self.l}|{self.m}|{self.gamma}|"+A_str+"|"+B_str+"|"


##  CHAPTER: ====Majorana Twisted Double Bivariant Bicycle Code====
class MajoranaTwistedDoubleBivariantBicycleCode(MajoranaCSSCode,TwistedDoubleBivariantBicycleCode):

    ##  SECTION: ----Constructor----
    def __init__(self,l:int,m:int,gamma:int,a_x_list:list[int],a_y_list:list[int],b_x_list:list[int],b_y_list:list[int]):
        TwistedDoubleBivariantBicycleCode.__init__(self,l,m,gamma,a_x_list,a_y_list,b_x_list,b_y_list)
        generators_x=[]
        generators_z=[]
        for i in range(len(self.H)):
            generators_x.append(MajoranaOperator.HermitianOperatorFromOccupy(np.where(self.H[i]!=0)[0],[]))
            generators_z.append(MajoranaOperator.HermitianOperatorFromOccupy([],np.where(self.H[i]!=0)[0]))
        MajoranaCSSCode.__init__(self,generators_x,generators_z,self.l*self.m*4)


##  CHAPTER: ====Pauli Twisted Double Bivariant Bicycle Code====
class PauliTwistedDoubleBivariantBicycleCode(PauliCSSCode,TwistedDoubleBivariantBicycleCode):

    ##  SECTION: ----Constructor----
    def __init__(self,l:int,m:int,gamma:int,a_x_list:list[int],a_y_list:list[int],b_x_list:list[int],b_y_list:list[int]):
        TwistedDoubleBivariantBicycleCode.__init__(self,l,m,gamma,a_x_list,a_y_list,b_x_list,b_y_list)
        generators_x=[]
        generators_z=[]
        for i in range(len(self.H)):
            generators_x.append(PauliOperator.HermitianOperatorFromOccupy(np.where(self.H[i]!=0)[0],[]))
            generators_z.append(PauliOperator.HermitianOperatorFromOccupy([],np.where(self.H[i]!=0)[0]))
        PauliCSSCode.__init__(self,generators_x,generators_z,self.l*self.m*4)


##  CHAPTER: ====Shift====
def shift(number:int, shift_number:int):
    """""
    input.number：dimension of the matrix, i.e. the number of qubits.
    input.shift：loop shift step, i.e. the number of qubits to shift.
    output：GF(2) matrix, the shift matrix.
    """""

    S = np.zeros((number, number), dtype=int)
    for i in range(number):
        S[i, (i + shift_number) % number] = 1
    return galois.GF(2)(S)
