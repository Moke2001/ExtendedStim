"""
Module: MajoranaCSSCode
A class for Majorana CSS code, each of which has X and Z stabilizer generators and physical number of modes as properties.
"""
import numpy as np
from extendedstim.Code.PrimitiveCode.MajoranaCode import MajoranaCode
from extendedstim.Code.PrimitiveCode.QuantumCSSCode import QuantumCSSCode
from extendedstim.Physics.MajoranaOperator import MajoranaOperator
from extendedstim.tools.GaloisTools import orthogonalize, mip_distance_caculator, minus


class MajoranaCSSCode(MajoranaCode, QuantumCSSCode):

    ##  SECTION: ----Constructor----
    def __init__(self, generators_x, generators_z, physical_number):
        """
        input.generators_x：list of MajoranaOperator, X Majorana stabilizer
        input.generators_z：list of MajoranaOperator, Z Majorana stabilizer
        input.physical_number：fermionic site number
        """

        ##  Check inputs
        try:
            generators_x: list[MajoranaOperator]=[generator if isinstance(generator, MajoranaOperator) else MajoranaOperator(generator) for generator in generators_x]
            generators_z: list[MajoranaOperator]=[generator if isinstance(generator, MajoranaOperator) else MajoranaOperator(generator) for generator in generators_z]
        except ValueError as e:
            raise ValueError("generators_x must be a list of MajoranaOperator.") from e

        try:
            physical_number=int(physical_number)
            assert physical_number>0, "physical_number must be a positive integer."
        except Exception as e:
            raise ValueError("physical_number must be a positive integer.") from e
        
        ##  Initialize code
        QuantumCSSCode.__init__(self, generators_x, generators_z, physical_number)


    ## SECTION: ----get code distance by mip method----
    @property
    def distance(self):
        """
        output：int, get code distance by mip method
        """

        ##  Return code distance
        return mip_distance_caculator(self.check_matrix,minus(self.check_matrix.null_space(),self.check_matrix))


    ##  SECTION: ----get \\hat gamma code distance by mip method----
    @property
    def distance_x(self):
        """
        output：int, get \\hat gamma code distance by mip method
        """

        ##  Return code distance
        return mip_distance_caculator(self.check_matrix_x,minus(self.check_matrix_x.null_space(),self.check_matrix_x))


    ##  SECTION: ----get \\hat gamma^\\prime code distance by mip method----
    @property
    def distance_z(self):
        """
        output：int, get \\hat gamma^\\prime code distance by mip method
        """

        ##  Return code distance
        return mip_distance_caculator(self.check_matrix_z,minus(self.check_matrix_z.null_space(),self.check_matrix_z))


    ##  SECTION: ----get logical operators----
    @property
    def logical_operators(self)->'list[MajoranaOperator]':
        """
        output：np.array[MajoranaOperator]，independent logical operators
        """

        ##  Return logical operators
        _=self._logical_operators_x
        return self._logical_operators_x+self._logical_operators_z


    ##  SECTION: ----get \\hat gamma logical operators----
    @property
    def logical_operators_x(self)->'list[MajoranaOperator]':
        """
        output：np.array[MajoranaOperator]，independent \\hat gamma logical operators
        """

        ##  Return \\gamma logical operators
        if self._logical_operators_x is None or len(self._logical_operators_x)==0:
            matrix = self.check_matrix_x
            codewords = matrix.null_space()
            independent_null_basis_list = []
            for vec in codewords:
                rank_before = np.linalg.matrix_rank(matrix)
                matrix = np.vstack([matrix, vec])
                if np.linalg.matrix_rank(matrix) == rank_before + 1:
                    independent_null_basis_list.append(vec)
            basis_list = orthogonalize(independent_null_basis_list)
            majorana_logical_operators_x = []
            majorana_logical_operators_z = []
            for i in range(len(basis_list)):
                occupy_temp=np.where(basis_list[i]==1)[0]
                temp = MajoranaOperator.HermitianOperatorFromOccupy(occupy_temp,[])
                majorana_logical_operators_x.append(temp)
                temp = MajoranaOperator.HermitianOperatorFromOccupy([],occupy_temp)
                majorana_logical_operators_z.append(temp)
            self._logical_operators_x = majorana_logical_operators_x
            self._logical_operators_z = majorana_logical_operators_z
            return self._logical_operators_x
        else:
            return self._logical_operators_x
        

    ##  SECTION: ----get \\hat gamma^\\prime logical operators----
    @property
    def logical_operators_z(self)->'list[MajoranaOperator]':
        """
        output：np.array[MajoranaOperator]，independent \\hat gamma^\\prime logical operators
        """

        ##  Return \\gamma^\\prime logical operators
        if self._logical_operators_z is None or len(self._logical_operators_z)==0:
            _=self.logical_operators_x
            return self._logical_operators_z
        else:
            return self._logical_operators_z

        
    ##  SECTION: ----get \\hat gamma^\\prime logical operators----
    @staticmethod
    def FromCheckMatrix(check_matrix):
        """
        input.check_matrix：GF(2) (m,2n)
        output：MajoranaCSSCode
        """
        
        ##  Check input
        try:
            check_matrix=np.array(check_matrix, dtype=int)
            assert check_matrix.shape[1]%2==0, "check_matrix must have even number of columns."
        except Exception as e:
            raise ValueError("check_matrix must be np.ndarray of GF(2).") from e

        ##  Generate generators
        generators_x = []
        generators_z = []
        for i in range(len(check_matrix)):
            generators_x.append(MajoranaOperator.HermitianOperatorFromVector(check_matrix[i]))
            generators_z.append(MajoranaOperator.HermitianOperatorFromVector(check_matrix[i]))
        physical_number=check_matrix.shape[1]

        ##  Initialize code
        return MajoranaCSSCode(generators_x, generators_z, physical_number)
