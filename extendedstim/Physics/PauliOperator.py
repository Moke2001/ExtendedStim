"""
Module: PauliOperator
A class for quantum Pauli operators, each of which has a coefficient and a string composed of local operators X and Z.
"""
import galois
import numpy as np
from extendedstim.Physics.Operator import Operator


class PauliOperator(Operator):

    ##  SECTION: ----Constructor----
    def __init__(self, occupy_x, occupy_z, coff):
        """
        input.occupy_x: Which index has local X operator in the operator.
        input.occupy_z: Which index has local Z operator in the operator.
        input.coff: Coefficient 1, -1, i, or -i.
        output：None
        """
        super().__init__(occupy_x, occupy_z, coff)


    ##  SECTION: ----Matrix product----
    def __matmul__(self, other:'PauliOperator')->'PauliOperator':
        """
        input.other: PauliOperator, the operator on the right.
        output: PauliOperator, the product of self and other.
        """

        ##  Check input
        assert isinstance(other, PauliOperator)

        ##  Computing occupy of new operator
        occupy_x = np.setxor1d(self.occupy_x, other.occupy_x, assume_unique=False)
        occupy_z = np.setxor1d(self.occupy_z, other.occupy_z, assume_unique=False)
        exchange_times = np.sum([np.count_nonzero(self.occupy_z == temp) for temp in other.occupy_x])
        
        ##  Computing coefficient of new operator
        if exchange_times % 2 == 1:
            factor = -1
        else:
            factor = 1

        ##  Return new operator
        return PauliOperator(occupy_x, occupy_z, self.coff * other.coff * factor)


    ##  SECTION: ----Left matrix multiplication----
    def __rmatmul__(self, other:'PauliOperator')->'PauliOperator':
        """
        input.other: PauliOperator, the operator on the left.
        output: PauliOperator, the product of self and other.
        """

        ##  Check input
        assert isinstance(other, PauliOperator)
        
        ##  Computing and return new operator
        return other.__matmul__(self)


    ##  SECTION: ----Multiplication with scalar----
    def __mul__(self, other:complex|float|int)->'PauliOperator':
        """
        input.other: complex|float|int, the scalar to multiply.
        output: PauliOperator, the product of self and other.
        """

        ##  Check input
        assert other == 1 or other == -1 or other == 1j or other == -1j
        
        ##  Computing and return new operator
        return PauliOperator(self.occupy_x, self.occupy_z, self.coff * other)


    ##  SECTION: ----Right multiplication with scalar----
    def __rmul__(self, other:complex|float|int)->'PauliOperator':
        """
        input.other: complex|float|int, the scalar to multiply.
        output: PauliOperator, the product of self and other.
        """

        ##  Computing and return new operator
        return self.__mul__(other)


    ##  SECTION: ----Equality----
    def __eq__(self, other:'PauliOperator')->bool:
        """
        input.other: PauliOperator, the operator to compare.
        output: bool, whether self and other are equal.
        """

        ##  Check input
        assert isinstance(other, PauliOperator)
        
        ##  Check equality and return result
        return np.array_equal(self.occupy_x, other.occupy_x) and np.array_equal(self.occupy_z, other.occupy_z) and self.coff == other.coff


    ##  SECTION: ----Negate----
    def __neg__(self)->'PauliOperator':
        """
        output: PauliOperator, the negation of self.
        """

        ##  Computing and return new operator
        return PauliOperator(self.occupy_x, self.occupy_z, -self.coff)


    ##  SECTION: ----Print----
    def __str__(self)->str:
        """
        output: str, the string representation of self.
        """

        ##  Preprocess
        op_str="{self.coff}"  # Initialize operator string
        flag_z=0  # Record current index of Z operator
        flag_x=0  # Record current index of X operator

        ##  Fill operator string by order from left to right, from X to Z
        while flag_x<len(self.occupy_x) and flag_z<len(self.occupy_z):
            if self.occupy_x[flag_x]>=self.occupy_z[flag_z]:
                op_str+=f"\\hat X_{self.occupy_x[flag_x]}"
                flag_x+=1
            else:
                op_str+=f"\\hat Z_{self.occupy_z[flag_z]}"
                flag_z+=1
        
        ##  Return operator string
        return op_str


    ##  SECTION: ----Check hermitian----
    @property
    def is_hermitian(self)->bool:
        """
        output：bool, whether self is hermitian.
        """

        ##  Check hermitian property, i.e. how many overlaps between one's X and the other's Z
        if len(np.intersect1d(self.occupy_x, self.occupy_z, assume_unique=False)) % 2 == 0:
            if not (self.coff == 1 or self.coff == -1):
                return False
        else:
            if not (self.coff == 1j or self.coff == -1j):
                return False

        ##  Return result
        return True


    ##  SECTION: ----Dual operator----
    @property
    def dual(self)->'PauliOperator':
        """
        output：PauliOperator, the dual operator of self.
        """

        ##  Compute and return dual operator
        return PauliOperator(self.occupy_z.copy(), self.occupy_x.copy(), self.coff)


    ##  SECTION: ----Hermitian Pauli operator constructing from specified occupy----
    @staticmethod
    def HermitianOperatorFromOccupy(occupy_x:list[int]|np.ndarray,occupy_z:list[int]|np.ndarray)->'PauliOperator':
        """
        input.occupy_x：which qubits have X operator
        input.occupy_z：which qubits have Z operator
        output：PauliOperator, the Hermitian operator with specified occupy.
        """

        ##  Check input
        try:
            occupy_x=np.array([temp for temp in occupy_x],dtype=int)
            occupy_z=np.array([temp for temp in occupy_z],dtype=int)
        except Exception as e:
            raise ValueError("occupy_x and occupy_z must be a list of integers.") from e

        ##  Compute coefficient
        weight=len(np.intersect1d(occupy_x,occupy_z, assume_unique=False))
        if weight % 2 == 0:
            coff=1
        else:
            coff=1j

        ##  Return result
        return PauliOperator(occupy_x,occupy_z,coff)


    ##  SECTION: ----Hermitian Pauli operator constructing from vector----
    @staticmethod
    def HermitianOperatorFromVector(vector)->'PauliOperator':
        """
        input.vector：vector representation of Pauli operator, even indices are X operators, odd indices are Z operators.
        output：PauliOperator, the Hermitian operator with specified vector.
        """

        ##  Check input
        try:
            GF2=galois.GF(2)
            vector=GF2(np.array([temp for temp in vector],dtype=int))
        except Exception as e:
            raise ValueError("vector must be a list of integers.") from e

        ##  Compute Pauli occupy
        pauli_x = np.where(vector[0::2] == 1)[0]
        pauli_z = np.where(vector[1::2] == 1)[0]

        ##  Return result
        return PauliOperator.HermitianOperatorFromOccupy(pauli_x, pauli_z)


    ##  SECTION: ----Check commutation----
    @staticmethod
    def commute(A:'PauliOperator',B:'PauliOperator')->bool:
        """
        input.A：first Pauli operator
        input.B：second Pauli operator
        output：bool, whether A and B commute.
        """

        ##  Check input
        assert isinstance(A, PauliOperator) and isinstance(B, PauliOperator)

        ##  Compute commutation
        same_time=np.sum([np.count_nonzero(A.occupy_z==temp) for temp in B.occupy_x])
        same_time+=np.sum([np.count_nonzero(A.occupy_x==temp) for temp in B.occupy_z])

        ##  Return result
        return np.mod(same_time, 2) == 0
