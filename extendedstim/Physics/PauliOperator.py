import numpy as np
from Physics.Operator import Operator


class PauliOperator(Operator):

    #%%  USER：构造方法
    def __init__(self, occupy_x, occupy_y, coff):
        super().__init__(occupy_x, occupy_y, coff)

    #%%  USER：重载运算符
    ##  USER：定义算符的乘积
    def __matmul__(self, other):
        assert isinstance(other, PauliOperator)
        occupy_x = np.setxor1d(self.occupy_x, other.occupy_x, assume_unique=False)
        occupy_z = np.setxor1d(self.occupy_z, other.occupy_z, assume_unique=False)
        exchange_times = np.sum([np.count_nonzero(self.occupy_z == temp) for temp in other.occupy_x])
        if exchange_times % 2 == 1:
            factor = -1
        else:
            factor = 1
        return PauliOperator(occupy_x, occupy_z, self.coff * other.coff * factor)

    ##  USER：定义算符的乘积
    def __rmatmul__(self, other):
        assert isinstance(other, PauliOperator)
        return other.__matmul__(self)

    ##  USER：定义算符的乘积
    def __mul__(self, other):
        assert other == 1 or other == -1 or other == 1j or other == -1j
        return PauliOperator(self.occupy_x, self.occupy_z, self.coff * other)

    ##  USER：定义算符的乘积
    def __rmul__(self, other):
        return self.__mul__(other)

    ##  USER：定义算符的字符串表示
    def __str__(self):
        return "PauliOperator(occupy_x={},occupy_z={},coff={})".format(self.occupy_x, self.occupy_z, self.coff)

    ##  USER：定义算符的相等性判断
    def __eq__(self, other):
        assert isinstance(other, PauliOperator)
        return np.array_equal(self.occupy_x, other.occupy_x) and np.array_equal(self.occupy_z, other.occupy_z) and self.coff == other.coff

    ##  USER：定义算符的负号
    def __neg__(self):
        return PauliOperator(self.occupy_x, self.occupy_z, -self.coff)

    #%%  USER：属性方法
    ##  USER：算符的类型
    @property
    def type(self):
        return "PauliOperator"

    ##  USER：算符是否是厄米算符
    @property
    def is_hermitian(self):
        if len(np.intersect1d(self.occupy_x, self.occupy_z, assume_unique=False)) % 2 == 0:
            if not (self.coff == 1 or self.coff == -1):
                return False
        else:
            if not (self.coff == 1j or self.coff == -1j):
                return False
        return True

    ##  USER：算符的对偶算符
    @property
    def dual(self):
        return PauliOperator(self.occupy_z, self.occupy_x, self.coff)

    #%%  USER：对象方法
    ##  USER：将算符在index处分裂
    def split(self, index):
        left_x,left_z,middle_x,middle_z,right_x,right_z = super().split(index)
        return PauliOperator(left_x, left_z, self.coff), PauliOperator(middle_x, middle_z, 1), PauliOperator(right_x, right_z, 1)

    def index_map(self, index):
        x,z,coff = super().index_map(index)
        return PauliOperator(x, z, coff)

    def copy(self):
        return PauliOperator(self.occupy_x.copy(), self.occupy_z.copy(), self.coff)

    #%%  USER：静态方法
    ##  USER：定义一个厄米算符
    @staticmethod
    def HermitianOperatorFromOccupy(occupy_x, occupy_z):
        weight=len(np.intersect1d(occupy_x,occupy_z, assume_unique=False))
        if weight % 2 == 0:
            coff=1
        else:
            coff=1j
        return PauliOperator(occupy_x,occupy_z,coff)

    ##  USER：定义一个厄米算符
    @staticmethod
    def HermitianOperatorFromVector(vector):
        pauli_x = np.where(vector[0::2] == 1)[0]
        pauli_z = np.where(vector[1::2] == 1)[0]
        return PauliOperator.HermitianOperatorFromOccupy(pauli_x, pauli_z)

    ##  USER：判断两个厄米算符是否对易
    @staticmethod
    def commute(A,B):
        assert isinstance(A, PauliOperator) and isinstance(B, PauliOperator)
        same_time=np.sum([np.count_nonzero(A.occupy_z==temp) for temp in B.occupy_x])
        same_time+=np.sum([np.count_nonzero(A.occupy_x==temp) for temp in B.occupy_z])
        return np.mod(same_time, 2) == 0
