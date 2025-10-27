import numpy as np

from .PauliCode import PauliCode
from .QuantumCSSCode import QuantumCSSCode
from ...Math.BinaryArray import BinaryArray as ba
from ...Physics.MajoranaOperator import MajoranaOperator


class PauliCSSCode(PauliCode, QuantumCSSCode):

    # %%  USER：构造方法
    def __init__(self, generators_x, generators_z, physical_number):
        QuantumCSSCode.__init__(self, generators_x, generators_z, physical_number)

    # %%  USER：属性方法
    ##  TODO：求Pauli CSS code的距离（x方向）
    @property
    def distance_x(self):
        return 1

    ##  TODO：求Pauli CSS code的距离（z方向）
    @property
    def distance_z(self):
        return 1

    ##  TODO：求Pauli CSS code的逻辑算子（x方向）
    @property
    def logical_operators_x(self):
        return None

    ##  TODO：求Pauli CSS code的逻辑算子（z方向）
    @property
    def logical_operators_z(self):
        return None
