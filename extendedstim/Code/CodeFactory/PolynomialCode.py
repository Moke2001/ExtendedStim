"""""
Module: PolynomialCode, construct polynomial code from polynomials.
"""""
import copy
import galois
from galois import GF2
from extendedstim.Code.PrimitiveCode.MajoranaCSSCode import MajoranaCSSCode
from extendedstim.Code.PrimitiveCode.PauliCSSCode import PauliCSSCode
from extendedstim.Physics.MajoranaOperator import MajoranaOperator
from extendedstim.Physics.PauliOperator import PauliOperator
from extendedstim.tools.GaloisTools import mip_distance_caculator, minus
from extendedstim.tools.TypingTools import islist, isinteger
import numpy as np


##  CHAPTER: ====Element of finite polynomial ring====
class Element:
    """
    An example:
    $x_0x_1x_2x_3+2x_0x_1x_2x_3+3x_0x_1x_2x_3$
    element = Element(4,[[0,0,0,1],[0,0,1,1],[0,1,1,1]],[2,2,2,2])
    """

    ##  SECTION: ----Constructor----
    def __init__(self,number:int,polynomial,lengths):
        assert isinteger(number)
        assert islist(polynomial)
        assert all([len(temp)==number for temp in polynomial])
        self.number = number
        x = np.array(polynomial,dtype=int)
        self.lengths = np.array(lengths)
        for i in range(x.shape[1]):
            x[:,i]=np.mod(x[:,i],lengths[i])
        row_dtype = np.dtype([("f", x.dtype, (x.shape[1],))])
        x_view = x.view(row_dtype).ravel()
        uniq, counts = np.unique(x_view, return_counts=True)
        self.polynomial = uniq[counts == 1].view(x.dtype).reshape(-1, x.shape[1])

    ##  SECTION: ----Add----
    def __add__(self,other:'Element'):
        if other==0:
            return self.copy()
        elif other==1:
            return Element(self.number,np.vstack([np.array([0]*self.number,dtype=int), other.polynomial]),self.lengths)
        else:
            assert len(self.lengths)==len(other.lengths)
            assert self.number == other.number and all([self.lengths[i] == other.lengths[i] for i in range(len(self.lengths))])
            return Element(self.number, np.vstack([self.polynomial, other.polynomial]), self.lengths)

    ##  SECTION: ----Multiply----
    def __mul__(self,other:'Element'):
        assert other==1 or self.number==other.number
        if other==1:
            return self.copy()
        elif other==0:
            return 0
        x = []
        for i in range(len(self.polynomial)):
            for j in range(len(other.polynomial)):
                x.append(self.polynomial[i]+other.polynomial[j])
        return Element(self.number, np.array(x,dtype=int), self.lengths)

    ##  SECTION: ----Power----
    def __pow__(self,other:int):
        assert isinteger(other)
        temp=self
        if other==0:
            return Element(self.number,np.array([[0]*self.number,],dtype=int),self.lengths)
        else:
            for i in range(other-1):
                temp=temp*self
            return temp

    ##  SECTION: ----Print----
    def __str__(self):
        poly_str=""
        for q in range(len(self.polynomial)):
            if q==0:
                for p in range(len(self.polynomial[q])):
                    if np.all(self.polynomial[q,p]==0):
                        poly_str+="1"
                    else:
                        if self.polynomial[q][p]==0:
                            continue
                        elif self.polynomial[q][p]==1:
                            poly_str=poly_str+f"x_{q}"
                        else:
                            poly_str=poly_str+f"x_{q}^{self.polynomial[q][p]}"
            else:
                poly_str = poly_str + f"+"
                for p in range(len(self.polynomial[q])):
                    if np.all(self.polynomial[q,p]==0):
                        poly_str+="1"
                    else:
                        if self.polynomial[q][p]==0:
                            continue
                        elif self.polynomial[q][p]==1:
                            poly_str=poly_str+f"x_{q}"
                        else:
                            poly_str=poly_str+f"x_{q}^{self.polynomial[q][p]}"

        return poly_str

    ##  SECTION: ----Copy----
    def copy(self):
        return Element(self.number,copy.deepcopy(self.polynomial),self.lengths)

    ##  SECTION: ----Return matrix representation----
    def matrix(self):
        result=0
        for i in range(len(self.polynomial)):
            temp=1
            for j in range(len(self.polynomial[i])):
                temp=np.kron(temp,shift(self.lengths[j],self.polynomial[i][j]))
            if i==0:
                result=GF2(temp)
            else:
                result=result+GF2(temp)
        return result

    ##  SECTION: ----Transpose----
    @property
    def T(self):
        result=self.copy()
        for i in range(result.polynomial.shape[1]):
            result.polynomial[:,i]=self.lengths[i]-result.polynomial[:,i]
        return result


##  CHAPTER: ====Polynomial====
class Polynomial:

    ##  SECTION: ----Constructor----
    def __init__(self,number,lengths):
        self.number=number
        self.lengths=lengths

    ##  SECTION: Return element of some variable of polynomial
    def element(self,index):
        temp=[0]*self.number
        temp[index]=1
        return Element(self.number,np.array([temp,],dtype=int),self.lengths)


##  CHAPTER: ====Check matrix of polynomial code====
class PolynomialCode:

    ##  SECTION: ----Constructor----
    def __init__(self, lengths:list[int]|np.ndarray,cell_size:int,poly_x_list:list[Element],poly_z_list:list[Element]):
        """
        lengths: Lengths of each dimension
        cell_size: Number of qubits in each unit cell
        polynomials: Polynomial description, structure A[q][i][j]
        q: Represents the qth qubit in the unit cell (0 <= q < cell_size)
        i: Represents the ith position of the qubit in the stabilizer (i.e., the ith term in the polynomial)
        j: Represents the dimension index of the offset vector [x, y, z...]
        """

        self.lengths = lengths
        self.cell_size = cell_size
        assert cell_size==len(poly_x_list)==len(poly_z_list)
        assert np.all([temp.lengths==self.lengths for temp in poly_x_list+poly_z_list])
        self.polynomials_x = poly_x_list
        self.polynomials_z = poly_z_list
        matrix_list_x = []
        matrix_list_z = []
        for i in range(self.cell_size):
            matrix_list_x.append(self.polynomials_x[i].matrix())
            matrix_list_z.append(self.polynomials_z[i].matrix())
        self.H_X=np.hstack(matrix_list_x)
        self.H_Z=np.hstack(matrix_list_z)


    ##  SECTION: ----Print----
    def __str__(self)->str:
        """
        input：None
        output：string, to represent the polynomial code
        influence：print the polynomial code
        """

        x_str="("
        for i in range(self.cell_size):
            if i==0:
                x_str=x_str+f"{self.polynomials_x[i]}"
            else:
                x_str=x_str+f",{self.polynomials_x[i]}"
        x_str=x_str+"|0,...,0)^T"

        z_str = "(0,...,0|"
        for i in range(self.cell_size):
            if i == 0:
                z_str = z_str + f"{self.polynomials_z[i]}"
            else:
                z_str = z_str + f",{self.polynomials_z[i]}"
        z_str = z_str + ")^T"

        return f"$H_X={x_str}, H_Z={z_str}$"


##  CHAPTER: ====Majorana polynomial Code====
class MajoranaPolynomialCode(MajoranaCSSCode,PolynomialCode):

    ##  SECTION: ----Constructor----
    def __init__(self,lengths:list[int]|np.ndarray,cell_size:int,poly_x_list:list[Element],poly_z_list:list[Element]):
        PolynomialCode.__init__(self,lengths,cell_size,poly_x_list,poly_z_list)
        generators_x=[]
        generators_z=[]
        for i in range(len(self.H_X)):
            generators_x.append(MajoranaOperator.HermitianOperatorFromOccupy(np.where(self.H_X[i]!=0)[0],[]))
        for i in range(len(self.H_Z)):
            generators_z.append(MajoranaOperator.HermitianOperatorFromOccupy([],np.where(self.H_Z[i]!=0)[0]))
        MajoranaCSSCode.__init__(self,generators_x,generators_z,self.H_X.shape[1])


##  CHAPTER: ====Pauli polynomial Code====
class PauliPolynomialCode(PauliCSSCode,PolynomialCode):

    ##  SECTION: ----Constructor----
    def __init__(self,lengths:list[int]|np.ndarray,cell_size:int,poly_x_list:list[Element],poly_z_list:list[Element]):
        PolynomialCode.__init__(self,lengths,cell_size,poly_x_list,poly_z_list)
        generators_x=[]
        generators_z=[]
        for i in range(len(self.H_X)):
            generators_x.append(PauliOperator.HermitianOperatorFromOccupy(np.where(self.H_X[i]!=0)[0],[]))
        for i in range(len(self.H_Z)):
            generators_z.append(PauliOperator.HermitianOperatorFromOccupy([],np.where(self.H_Z[i]!=0)[0]))
        PauliCSSCode.__init__(self,generators_x,generators_z,self.H_X.shape[1])


##  CHAPTER: ====Tools for constructing polynomial code====
def shift(number:int, shift:int):
    """
    input.number: Matrix dimension
    input.shift: Cyclic shift step length
    output: GF(2) matrix, permutation matrix
    """

    S = np.zeros((number, number), dtype=int)
    for i in range(number):
        S[i, (i + shift) % number] = 1
    return galois.GF(2)(S)