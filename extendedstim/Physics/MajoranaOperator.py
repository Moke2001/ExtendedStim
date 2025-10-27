import numpy as np
from Physics.Operator import Operator


class MajoranaOperator(Operator):

    #%%  USER：构造方法
    def __init__(self, occupy_x, occupy_y, coff):
        super().__init__(occupy_x, occupy_y, coff)

    #%%  USER：重载运算符
    ##  USER：矩阵乘法
    def __matmul__(self, other):
        assert isinstance(other, MajoranaOperator)
        occupy_x = np.setxor1d(self.occupy_x, other.occupy_x, assume_unique=False)
        occupy_z = np.setxor1d(self.occupy_z, other.occupy_z, assume_unique=False)
        self_occupy = np.append(self.occupy_x * 2, self.occupy_z * 2 + 1)
        other_occupy = np.append(other.occupy_x * 2, other.occupy_z * 2 + 1)
        self_occupy = np.sort(self_occupy)
        other_occupy = np.sort(other_occupy)
        exchange_times = np.sum([np.count_nonzero(self_occupy < temp) for temp in other_occupy])
        if exchange_times % 2 == 1:
            factor = -1
        else:
            factor = 1
        return MajoranaOperator(occupy_x, occupy_z, self.coff * other.coff * factor)

    ##  USER：右矩阵乘法
    def __rmatmul__(self, other):
        assert isinstance(other, MajoranaOperator)
        return other.__matmul__(self)

    ##  USER：标量乘法
    def __mul__(self, other):
        assert other == 1 or other == -1 or other == 1j or other == -1j
        return MajoranaOperator(self.occupy_x, self.occupy_z, self.coff * other)

    ##  USER：右标量乘法
    def __rmul__(self, other):
        return self.__mul__(other)

    ##  USER：字符串表示
    def __str__(self):
        return "MajoranaOperator(occupy_x={},occupy_z={},coff={})".format(self.occupy_x, self.occupy_z, self.coff)

    ##  USER：相等判断
    def __eq__(self, other):
        assert isinstance(other, MajoranaOperator)
        return np.array_equal(self.occupy_x, other.occupy_x) and np.array_equal(self.occupy_z, other.occupy_z) and self.coff == other.coff

    ##  USER：取负
    def __neg__(self):
        return MajoranaOperator(self.occupy_x, self.occupy_z, -self.coff)

    # %%  USER：属性方法
    ##  USER：算符的类型
    @property
    def type(self):
        return "MajoranaOperator"

    ##  USER：算符是否是厄米算符
    @property
    def is_hermitian(self):
        if self.weight * (self.weight - 1) // 2 == 0:
            if not (self.coff == 1 or self.coff == -1):
                return False
        else:
            if not (self.coff == 1j or self.coff == -1j):
                return False
        return True

    ##  USER：求算符的对偶算符
    @property
    def dual(self):
        return MajoranaOperator(self.occupy_z, self.occupy_x, self.coff)

    # %%  USER：对象方法
    ##  USER：将算符在index位置切分
    def split(self, index):
        left_x, left_z, middle_x, middle_z, right_x, right_z = super().split(index)
        return MajoranaOperator(left_x, left_z, self.coff), MajoranaOperator(middle_x, middle_z, 1), MajoranaOperator(right_x, right_z, 1)

    def index_map(self, index):
        x, z, coff = super().index_map(index)
        return MajoranaOperator(x, z, coff)

    def copy(self):
        return MajoranaOperator(self.occupy_x.copy(), self.occupy_z.copy(), self.coff)

    # %%  USER：静态方法
    ##  USER：定义一个厄米算符
    @staticmethod
    def HermitianOperatorFromOccupy(occupy_x, occupy_z):
        weight = len(occupy_x) + len(occupy_z)
        if (weight * (weight - 1) // 2) % 2 == 0:
            coff = 1
        else:
            coff = 1j
        return MajoranaOperator(occupy_x, occupy_z, coff)

    ##  USER：定义一个厄米算符
    @staticmethod
    def HermitianOperatorFromVector(vector):
        occupy_x = np.where(vector[0::2] == 1)[0]
        occupy_z = np.where(vector[1::2] == 1)[0]
        return MajoranaOperator.HermitianOperatorFromOccupy(occupy_x, occupy_z)

    ##  USER：检查两个厄米算符是否对易
    @staticmethod
    def commute(A, B):
        assert isinstance(A, MajoranaOperator) and isinstance(B, MajoranaOperator)
        overlap_x = len(np.intersect1d(A.occupy_x, B.occupy_x))
        overlap_z = len(np.intersect1d(A.occupy_z, B.occupy_z))
        weight = (len(A.occupy_x) + len(A.occupy_z)) * (len(B.occupy_x) + len(B.occupy_z))
        judge = overlap_x + overlap_z + weight
        return np.mod(judge, 2) == 0
