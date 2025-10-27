import numpy as np
from Code.QuantumCode.MajoranaCode import MajoranaCode
from Code.QuantumCode.QuantumCSSCode import QuantumCSSCode
from Math.BinaryArray import BinaryArray as ba
from Physics.MajoranaOperator import MajoranaOperator


class MajoranaCSSCode(MajoranaCode, QuantumCSSCode):

    # %%  USER：构造方法
    def __init__(self, generators_x, generators_z, physical_number):
        QuantumCSSCode.__init__(self, generators_x, generators_z, physical_number)

    # %%  USER：属性方法
    ##  USER：求码距（x方向）
    @property
    def distance_x(self):
        return ba.distance(self.check_matrix_x,'mip')

    ##  USER：求码距（z方向）
    @property
    def distance_z(self):
        return ba.distance(self.check_matrix_z,'mip')

    ##  USER：求逻辑算符（x方向）
    @property
    def logical_operators_x(self):
        matrix = self.check_matrix_x()
        codewords = matrix.null_space()
        independent_null_basis_list = []
        for vec in codewords:
            rank_before = matrix.rank()
            matrix = ba.vstack(matrix, vec)
            if matrix.rank == rank_before + 1:
                independent_null_basis_list.append(vec)
        basis_list = ba.orthogonalize(independent_null_basis_list)
        majorana_logical_operators = []
        for i in range(len(basis_list)):
            occupy_x=np.where(basis_list[i]==1)[0]
            temp = MajoranaOperator.HermitianOperatorFromOccupy(occupy_x,[])
            majorana_logical_operators.append(temp)
        majorana_logical_operators = np.array(majorana_logical_operators, dtype=MajoranaOperator)
        return majorana_logical_operators

    ##  USER：求逻辑算符（z方向）
    @property
    def logical_operators_z(self):
        matrix = self.check_matrix_z()
        codewords = matrix.null_space()
        independent_null_basis_list = []
        for vec in codewords:
            rank_before = matrix.rank()
            matrix = ba.vstack(matrix, vec)
            if matrix.rank == rank_before + 1:
                independent_null_basis_list.append(vec)
        basis_list = ba.orthogonalize(independent_null_basis_list)
        majorana_logical_operators = []
        for i in range(len(basis_list)):
            occupy_z = np.where(basis_list[i] == 1)[0]
            temp = MajoranaOperator.HermitianOperatorFromOccupy([], occupy_z)
            majorana_logical_operators.append(temp)
        majorana_logical_operators = np.array(majorana_logical_operators, dtype=MajoranaOperator)
        return majorana_logical_operators
