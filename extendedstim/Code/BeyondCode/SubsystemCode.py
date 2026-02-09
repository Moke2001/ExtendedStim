"""
Module: SubsystemCode
"""
from extendedstim.Code.PrimitiveCode.MajoranaCode import MajoranaCode
from extendedstim.Code.PrimitiveCode.PauliCode import PauliCode
from extendedstim.Physics.MajoranaOperator import MajoranaOperator
from extendedstim.Physics.PauliOperator import PauliOperator
from extendedstim.tools.GaloisTools import mip_distance_caculator


class SubsystemCode:

    ##  SECTION: ----Constructor----
    def __init__(self, code: MajoranaCode | PauliCode, bare_logical_operators: list[MajoranaOperator | PauliOperator]):
        """
        input.code: MajoranaCode or PauliCode instance.
        input.bare_logical_operators: List of MajoranaOperator or PauliOperator, serving as measurement operators on the code space.
        output: None.
        influence: None
        """
        self.code = code
        self.bare_logical_operators = bare_logical_operators
        if isinstance(code, MajoranaCode):
            self.code_type = MajoranaCode
            self.operator_type = MajoranaOperator
        elif isinstance(code, PauliCode):
            self.code_type = PauliCode
            self.operator_type = PauliOperator
        else:
            raise ValueError("code must be MajoranaCode or PauliCode")


    ##  SECTION: ----Details----
    @property
    def logical_operators(self) -> list[MajoranaOperator | PauliOperator]:
        """
        output: list[MajoranaOperator|PauliOperator], set of independent logical operators.
        """
        return self.bare_logical_operators


    ##  SECTION: ----Details----
    @property
    def distance(self) -> int:
        """
        output: int, code distance.
        """
        return mip_distance_caculator(self.code.check_matrix, self.operator_type.get_matrix(self.bare_logical_operators, self.code.physical_number))
