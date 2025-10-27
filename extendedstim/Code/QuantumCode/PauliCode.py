import numpy as np
from .QuantumCode import QuantumCode
from ...Physics.MajoranaOperator import MajoranaOperator


class PauliCode(QuantumCode):

    # %%  USER：构造方法
    def __init__(self, generators, physical_number):
        super().__init__(generators, physical_number)

    # %%  USER：属性方法
    ##  TODO：求Pauli code的距离
    @property
    def distance(self):
        return 1

    ##  TODO：求Pauli code的逻辑算子
    @property
    def logical_operators(self):
        return []

    # %%  USER：静态方法
    ##  USER：基于校验矩阵构造code
    @staticmethod
    def FromCheckMatrix(check_matrix):
        generators = np.empty(check_matrix.shape[0], dtype=MajoranaOperator)
        for temp in range(check_matrix.shape[0]):
            occupy_x = np.where(check_matrix[temp, 0::2] == 1)
            occupy_z = np.where(check_matrix[temp, 1::2] == 1)
            generators[temp] = MajoranaOperator.HermitianOperatorFromOccupy(occupy_x, occupy_z)
        physical_number = check_matrix.shape[1] // 2
        return PauliCode(generators, physical_number)
