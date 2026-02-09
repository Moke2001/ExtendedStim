"""
Module: MajoranaCode
A class for Majorana code, each of which has stabilizer generators and physical number of modes as properties.
"""
import galois
import numpy as np
from extendedstim.Code.PrimitiveCode.QuantumCode import QuantumCode
from extendedstim.Physics.MajoranaOperator import MajoranaOperator
from extendedstim.tools.GaloisTools import orthogonalize, solve, mip_distance_caculator, minus


class MajoranaCode(QuantumCode):

    ##  SECTION: ----Constructor----
    def __init__(self,generators,physical_number):
        """
        input.generators：list of MajoranaOperator, stabilizer generators
        input.physical_number: number of physical fermionic sites
        output：None
        """

        ##  Check input
        try:
            generators=[generator if isinstance(generator, MajoranaOperator) else MajoranaOperator(generator) for generator in generators]
        except Exception as e:
            raise ValueError("Invalid generators") from e
        
        try:
            physical_number=int(physical_number)
            assert physical_number>=0, "physical_number must be non-negative"
        except Exception as e:
            raise ValueError("Invalid physical_number") from e

        ##  Initialize super class
        super().__init__(generators, physical_number)


    ##  SECTION: ----Get code distance----
    @property
    def distance(self) -> int:
        """
        output：int, get code distance by mip method
        """

        ##  Calculate code distance
        return mip_distance_caculator(self.check_matrix,minus(self.check_matrix.null_space(),self.check_matrix))


    ##  SECTION: ----Get logical operators----
    @property
    def logical_operators(self):
        """
        output：np.array[MajoranaOperator] independent logical operators
        """

        ##  Get logical operators from check matrix
        matrix = self.check_matrix
        codewords = matrix.null_space()
        independent_null_basis_list = []
        for vec in codewords:
            rank_before = np.linalg.matrix_rank(matrix)
            matrix = np.vstack([matrix, vec])
            if np.linalg.matrix_rank(matrix) == rank_before + 1:
                independent_null_basis_list.append(vec)
        basis_list = orthogonalize(independent_null_basis_list)
        majorana_logical_operators = []
        for i in range(len(basis_list)):
            temp = MajoranaOperator.HermitianOperatorFromVector(basis_list[i])
            majorana_logical_operators.append(temp)
        majorana_logical_operators = np.array(majorana_logical_operators, dtype=MajoranaOperator)
        
        ##  Return logical operators
        return majorana_logical_operators


    ##  SECTION: ----get even or odd----
    @property
    def even_or_odd(self):
        """
        output："even" or "odd", check if there exists a solution H x = 1
        """

        ##  Check if there exists a solution H x = 1
        H=self.check_matrix
        ones=galois.GF2(np.ones(H.shape[1],dtype=int))

        ##  Return "even" or "odd"
        if solve(H,ones) is None:
            return "odd"
        else:
            return "even"


    ##  SECTION: ----Get code from check matrix----
    @staticmethod
    def FromCheckMatrix(check_matrix):
        """
        input.H：GF(2) (m,2n)
        output：MajoranaCode 实例
        """

        ##  Check input
        try:
            check_matrix=galois.GF2(np.array(check_matrix,dtype=int))
        except Exception as e:
            raise ValueError("Invalid check_matrix") from e

        ##  Get Majorana stabilizer generators from checker in check matrix
        generators = []
        for temp in range(check_matrix.shape[0]):
            occupy_x=np.where(check_matrix[temp,0::2]==1)[0]
            occupy_z=np.where(check_matrix[temp,1::2]==1)[0]
            generators.append(MajoranaOperator.HermitianOperatorFromOccupy(occupy_x,occupy_z))
        physical_number=check_matrix.shape[1]
        
        ##  Return MajoranaCode instance
        return MajoranaCode(generators,physical_number)
