import copy
from abc import ABC, abstractmethod
import numpy as np
from Physics.Operator import Operator


class QuantumCode(ABC):

    #%%  USER：构造方法
    def __init__(self,generators,physical_number):
        assert isinstance(generators, list) or isinstance(generators, np.ndarray)
        self.generators = np.array(generators,dtype=Operator)
        self.physical_number=physical_number
        self.checker_number=self.generators.shape[0]

    #%%  USER：属性方法
    ##  USER：求校验矩阵
    @property
    def check_matrix(self):
        matrix=Operator.get_matrix(self.generators,self.physical_number)
        return matrix

    ##  USER：求校验矩阵的秩
    @property
    def rank(self):
        return self.check_matrix.rank

    ##  USER：求logical number
    @property
    def logical_number(self):
        return self.physical_number-self.rank

    ##  USER：求码距
    @property
    @abstractmethod
    def distance(self):
        pass

    ##  USER：求逻辑算符
    @property
    @abstractmethod
    def logical_operators(self):
        pass

    #%%  USER：对象方法
    ##  USER：修改索引
    def index_map(self, index_map):
        for generator in self.generators:
            generator.index_map(index_map)

    ##  USER：复制代码
    def copy(self):
        return copy.deepcopy(self)