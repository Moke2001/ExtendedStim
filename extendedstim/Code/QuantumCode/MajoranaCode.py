import numpy as np

from .QuantumCode import QuantumCode
from ...Math.BinaryArray import BinaryArray as ba
from ...Physics.MajoranaOperator import MajoranaOperator


class MajoranaCode(QuantumCode):

    #%%  USER：构造方法
    def __init__(self,generators,physical_number):
        super().__init__(generators, physical_number)

    #%%  USER：属性方法
    ##  USER：求码距
    @property
    def distance(self):
        return ba.distance(self.check_matrix,'mip')

    ##  USER：求逻辑算符
    @property
    def logical_operators(self):
        matrix = self.check_matrix()
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
            temp = MajoranaOperator.HermitianOperatorFromVector(basis_list[i])
            majorana_logical_operators.append(temp)
        majorana_logical_operators = np.array(majorana_logical_operators, dtype=MajoranaOperator)
        return majorana_logical_operators

    ##  USER：判断是否为偶数码
    @property
    def even_or_odd(self):
        H=self.check_matrix
        ones=ba.ones(H.shape[1])
        if ba.solve(H,ones) is None:
            return "odd"
        else:
            return "even"

    #%%  USER：静态方法
    ##  USER：基于校验矩阵构造code
    @staticmethod
    def FromCheckMatrix(check_matrix):
        generators = np.empty(check_matrix.shape[0],dtype=MajoranaOperator)
        for temp in range(check_matrix.shape[0]):
            occupy_x=np.where(check_matrix[temp,0::2]==1)
            occupy_z=np.where(check_matrix[temp,1::2]==1)
            generators[temp]=MajoranaOperator.HermitianOperatorFromOccupy(occupy_x,occupy_z)
        physical_number=check_matrix.shape[1]//2
        return MajoranaCode(generators,physical_number)
