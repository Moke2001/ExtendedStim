from abc import abstractmethod
import numpy as np
from QuantumCode import QuantumCode
from Physics.Operator import Operator


class QuantumCSSCode(QuantumCode):

    #%%  USER：构造方法
    def __init__(self, generators_x, generators_z, physical_number):
        self.generators_x = np.array(generators_x,dtype=Operator)
        self.generators_z = np.array(generators_z,dtype=Operator)
        self.checker_number_x=len(generators_x)
        self.checker_number_z = len(generators_z)
        super().__init__(generators_x+generators_z, physical_number)

    #%%  USER：属性方法
    ##  USER：求逻辑算符（X方向）
    @property
    @abstractmethod
    def logical_operators_x(self):
        pass

    ##  USER：求逻辑算符（Z方向）
    @property
    @abstractmethod
    def logical_operators_z(self):
        pass

    ##  USER：求校验矩阵的秩（X方向）
    @property
    def rank_x(self):
        return self.check_matrix_x.rank

    ##  USER：求校验矩阵的秩（Z方向）
    @property
    def rank_z(self):
        return self.check_matrix_z.rank

    ##  USER：求校验矩阵（X方向）
    @property
    def check_matrix_x(self):
        matrix=Operator.get_matrix(self.generators_x,self.physical_number)
        return matrix

    ##  USER：求校验矩阵（Z方向）
    @property
    def check_matrix_z(self):
        matrix=Operator.get_matrix(self.generators_z,self.physical_number)
        return matrix

    ##  USER：求码距（X方向）
    @property
    @abstractmethod
    def distance_x(self):
        pass

    ##  USER：求码距（Z方向）
    @property
    @abstractmethod
    def distance_z(self):
        return 1
