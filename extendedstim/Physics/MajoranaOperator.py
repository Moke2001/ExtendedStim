"""
Module: MajoranaOperator
A class for quantum Majorana operators, each of which has a coefficient and a string composed of local operators $\\hat \\gamma$ and $\\hat\\gamma^\\prime$.
"""
import galois
import numpy as np
from extendedstim.Physics.Operator import Operator


class MajoranaOperator(Operator):

    ##  SECTION: ----Constructor----
    def __init__(self, occupy_x, occupy_z, coff):
        """
        input.occupy_x: Which index has local $\\hat\\gamma$ operator in the operator.
        input.occupy_z: Which index has local $\\hat\\gamma^\\prime$ operator in the operator.
        input.coff: Coefficient 1, -1, i, or -i.
        output：None
        """
        super().__init__(occupy_x, occupy_z, coff)


    ##  SECTION: ----Matrix product----
    def __matmul__(self, other:'MajoranaOperator')->'MajoranaOperator':
        """
        input.other: MajoranaOperator, the operator on the right.
        output: MajoranaOperator, the product of self and other.
        """

        ##  Check input
        assert isinstance(other, MajoranaOperator)

        ##  Compute new occupy
        occupy_x = np.setxor1d(self.occupy_x, other.occupy_x, assume_unique=True)
        occupy_z = np.setxor1d(self.occupy_z, other.occupy_z, assume_unique=True)
        self_occupy = np.append(self.occupy_x * 2, self.occupy_z * 2 + 1)
        other_occupy = np.append(other.occupy_x * 2, other.occupy_z * 2 + 1)
        self_occupy = np.sort(self_occupy)
        other_occupy = np.sort(other_occupy)
        
        ##  Compute new coefficient
        exchange_times = np.sum([np.count_nonzero(self_occupy > temp) for temp in other_occupy])
        if exchange_times % 2 == 1:
            factor = -1
        else:
            factor = 1

        ##  Return new operator
        return MajoranaOperator(occupy_x, occupy_z, self.coff * other.coff * factor)


    ##  SECTION: ----Left matrix multiplication----
    def __rmatmul__(self, other:'MajoranaOperator')->'MajoranaOperator':
        """
        input.other: MajoranaOperator, the operator on the left.
        output: MajoranaOperator, the product of other and self.
        """

        ##  Compute and return new operator
        return other.__matmul__(self)


    ##  SECTION: ----Multiplication with scalar----
    def __mul__(self, other:complex|float|int)->'MajoranaOperator':
        """
        input.other: complex|float|int, the scalar to multiply.
        output: MajoranaOperator, the product of self and other.
        """

        ##  Check input
        assert other == 1 or other == -1 or other == 1j or other == -1j, "Scalar must be 1, -1, i, or -i."
        
        ##  Compute and return new operator
        return MajoranaOperator(self.occupy_x, self.occupy_z, self.coff * other)


    ##  SECTION: ----Right multiplication with scalar----
    def __rmul__(self, other:complex|float|int)->'MajoranaOperator':
        """
        input.other: complex|float|int, the scalar to multiply.
        output: MajoranaOperator, the product of other and self.
        """

        ##  Compute and return new operator
        return self.__mul__(other)


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
                op_str+=f"\\hat\\gamma_{self.occupy_x[flag_x]}"
                flag_x+=1
            else:
                op_str+=f"\\hat\\gamma^\\prime_{self.occupy_z[flag_z]}"
                flag_z+=1
        
        ##  Return operator string
        return op_str


    ##  SECTION: ----Equality----
    def __eq__(self, other:'MajoranaOperator')->bool:
        """
        input.other: MajoranaOperator, the operator to compare.
        output: bool, True if self and other are equal.
        """

        ##  Check input
        assert isinstance(other, MajoranaOperator)
        
        ##  Check equality and return result
        return np.array_equal(self.occupy_x, other.occupy_x) and np.array_equal(self.occupy_z, other.occupy_z) and self.coff == other.coff


    ##  SECTION: ----Negation----
    def __neg__(self)->'MajoranaOperator':
        """
        output: MajoranaOperator, the negation of self.
        """

        ##  Compute and return new operator
        return MajoranaOperator(self.occupy_x, self.occupy_z, -self.coff)


    ##  SECTION: ----Check Hermitian----
    @property
    def is_hermitian(self)->bool:
        """
        output：bool, whether self is hermitian.
        """

        ##  Check hermitian property
        if np.mod(self.weight * (self.weight - 1) // 2,2) == 0:
            if not (self.coff == 1 or self.coff == -1):
                return False
        else:
            if not (self.coff == 1j or self.coff == -1j):
                return False

        ##  Return result
        return True


    ##  SECTION: ----Get dual Majorana operator----
    @property
    def dual(self)->'MajoranaOperator':
        """
        output：MajoranaOperator, the dual of self.
        """

        ##  Compute and return new operator
        return MajoranaOperator(self.occupy_z, self.occupy_x, self.coff)


    ##  SECTION: ----Create Hermitian operator from occupy----
    @staticmethod
    def HermitianOperatorFromOccupy(occupy_x:list[int]|np.ndarray,occupy_z:list[int]|np.ndarray)->'MajoranaOperator':
        """
        input.occupy_x：which qubits have $\\hat\\gamma$ operator
        input.occupy_z：which qubits have $\\hat\\gamma^\\prime$ operator
        output：MajoranaOperator, the Hermitian operator with specified occupy.
        """

        ##  Check input
        try:
            occupy_x=np.array([temp for temp in occupy_x],dtype=int)
            occupy_z=np.array([temp for temp in occupy_z],dtype=int)
        except Exception as e:
            raise ValueError("occupy_x and occupy_z must be a list of integers.") from e

        ##  Compute coefficient
        weight = len(occupy_x) + len(occupy_z)
        if (weight * (weight - 1) // 2) % 2 == 0:
            coff = 1
        else:
            coff = 1j

        ##  Return new operator
        return MajoranaOperator(occupy_x, occupy_z, coff)


    ##  SECTION: ----Create Hermitian operator from vector----
    @staticmethod
    def HermitianOperatorFromVector(vector)->'MajoranaOperator':
        """
        input.vector：vector representation of Majorana operator, even indices are $\\hat\\gamma$ operators, odd indices are $\\hat\\gamma^\\prime$ operators.
        output：MajoranaOperator, the Hermitian operator with specified vector.
        """

        ##  Check input
        try:
            GF2=galois.GF(2)
            vector=GF2(np.array([temp for temp in vector],dtype=int))
        except Exception as e:
            raise ValueError("vector must be a list of integers.") from e

        ##  Compute occupy
        occupy_x = np.where(vector[0::2] == 1)[0]
        occupy_z = np.where(vector[1::2] == 1)[0]

        ##  Return new operator
        return MajoranaOperator.HermitianOperatorFromOccupy(occupy_x, occupy_z)


    ##  SECTION: ----Check commute----
    @staticmethod
    def commute(A:'MajoranaOperator', B:'MajoranaOperator')->bool:
        """
        input.A：first Majorana operator
        input.B：second Majorana operator
        output：bool, whether A and B commute.
        """

        ##  Check input
        assert isinstance(A, MajoranaOperator) and isinstance(B, MajoranaOperator)

        ##  Compute overlap
        overlap_x = len(np.intersect1d(A.occupy_x, B.occupy_x))
        overlap_z = len(np.intersect1d(A.occupy_z, B.occupy_z))
        weight = (len(A.occupy_x) + len(A.occupy_z)) * (len(B.occupy_x) + len(B.occupy_z))
        judge = overlap_x + overlap_z + weight

        ##  Check commute property
        return np.mod(judge, 2) == 0
