"""""
Module: Operator
An abstract class for quantum operators, each of which has a coefficient and a string composed of local operators.
"""""
from copy import copy
from abc import abstractmethod, ABC
import galois
import numpy as np
from numpy._typing import NDArray


class Operator(ABC):
    """
    An example:
    $-i\\hat A_0\\hat A_2\\hat B_3\\hat A_3\\hat B_3$
    op = Operator([0, 2, 3], [2, 3], -1j)
    """

    __slots__ = ['occupy_x', 'occupy_z', 'coff']

    ##  SECTION: ----Constructor----
    def __init__(self, occupy_x: list[int] | NDArray[np.integer], occupy_z: list[int] | NDArray[np.integer], coff: int | float | complex):
        """
        input.occupy_x: Which index has local A operator in the operator.
        input.occupy_z: Which index has local B operator in the operator.
        input.coff: Coefficient 1, -1, i, or -i.
        output：None
        """

        ##  Check input
        try:
            occupy_x = np.array([temp for temp in occupy_x], dtype=int)
            occupy_z = np.array([temp for temp in occupy_z], dtype=int)
        except Exception as e:
            raise ValueError("occupy_x and occupy_z must be a list of integers.") from e

        assert coff == -1 or coff == 1 or coff == 1j or coff == -1j, "coff must be 1, -1, 1j, or -1j."

        ##  Assign value
        self.occupy_x = occupy_x
        self.occupy_z = occupy_z
        self.occupy_x.sort()
        self.occupy_z.sort()
        self.coff = coff

    ##  SECTION: ----Abstract matrix multiply----
    @abstractmethod
    def __matmul__(self, other: 'Operator') -> 'Operator':
        """
        input.other：Another operator of the same type.
        output：Product operator.
        """

        pass

    ##  SECTION: ----Right matrix multiply----
    @abstractmethod
    def __rmatmul__(self, other: 'Operator') -> 'Operator':
        """
        input.other：Left operand.
        output：Product operator.
        """

        pass

    ##  SECTION: ----Multiply with scalar----
    @abstractmethod
    def __mul__(self, other: int | float | complex) -> 'Operator':
        """
        input.other：Coefficient (±1, ±i)
        output：New operator (coefficient is scaled)
        """

        pass

    ##  SECTION: ----Left Multiply----
    @abstractmethod
    def __rmul__(self, other: int | float | complex) -> 'Operator':
        """
        input.other：Coefficient (±1, ±i)
        output：New operator (coefficient is scaled)
        """

        pass

    ##  SECTION: ----Print----
    @abstractmethod
    def __str__(self) -> str:
        """
        output：str，可读表示
        """

        pass

    ##  SECTION: ----Compare----
    @abstractmethod
    def __eq__(self, other) -> bool:
        """
        input.other：Another operator of the same type.
        output：bool, True if the two operators have the same structure and coefficient.
        """

        pass

    ##  SECTION: ----Negate----
    @abstractmethod
    def __neg__(self):
        """
        output：New operator with negated coefficient.
        """

        pass

    ##  SECTION: ----Weight----
    @property
    def weight(self) -> int:
        """
        output: int, number of local operators in this operator
        """
        return len(self.occupy_x) + len(self.occupy_z)

    ##  SECTION: Details
    @property
    @abstractmethod
    def is_hermitian(self) -> bool:
        """
        output: bool, True if the operator is hermitian.
        """

        pass

    ##  SECTION: ----Dual----
    @property
    @abstractmethod
    def dual(self) -> None:
        """
        output：Dual operator (X/Z swapped, keep coefficient rule).
        """

        pass

    ##  SECTION: ----Index Map----
    def index_map(self, index) -> None:
        """
        input.index: New index mapping (array or list).
        output：None (modify object in place)
        An example:
        Map $-i\\hat X_0\\hat X_2\\hat Z_3\\hat X_3\\hat Z_3$ on 4 qubits from [0,2,3,3] to [1,3,4,4]:
        op.index_map([1,3,4,4])
        """

        ##  Check input
        assert isinstance(index, np.ndarray) or isinstance(index, list)
        if len(self.occupy_x) > 0:
            assert len(index) >= self.occupy_x.max() + 1
        if len(self.occupy_z) > 0:
            assert len(index) >= self.occupy_z.max() + 1

        ##  Deal with index
        x = np.array([index[i] for i in self.occupy_x], dtype=int)
        z = np.array([index[i] for i in self.occupy_z], dtype=int)
        self.occupy_x = x
        self.occupy_z = z

    ##  SECTION: ----Get vector represent of this operator----
    def get_vector(self, number: int):
        """
        input.number: Total number of qubits.
        output：Vector representation of this operator (GF(2) length 2*number)
        """

        ##  Check input
        try:
            number = int(number)
        except Exception as e:
            raise ValueError("number must be an integer.") from e

        ##  Create vector
        vector = galois.GF2(np.zeros(number * 2, dtype=int))
        vector[self.occupy_x * 2] = 1
        vector[self.occupy_z * 2 + 1] = 1

        ##  Return vector on GF(2)
        return vector

    ##  SECTION: ----Deepcopy----
    def copy(self):
        """
        output：deepcopy of new instance
        """
        return copy.deepcopy(self)

    ##  SECTION: ----Get matrix represent of these operators----
    @staticmethod
    def get_matrix(ops, number: int):
        """
        input.ops: List of Operators, each operator is a subclass of Operator.
        input.number: Total number of qubits.
        output：Matrix representation of these operators (GF(2) shape=(len(ops),2*number))
        """

        ##  Check input
        try:
            number = int(number)
        except Exception as e:
            raise ValueError("number must be an integer.") from e
        if not isinstance(ops, list):
            raise ValueError("ops must be a list of Operator.")
        for op in ops:
            if not isinstance(op, Operator):
                raise ValueError("ops must be a list of Operator.")

        ##  Vertically stack each vector as a matrix
        matrix = None
        for i, op in enumerate(ops):
            vector = op.get_vector(number)
            if i == 0:
                matrix = vector
            else:
                matrix = np.vstack([matrix, vector])

        if matrix is not None and matrix.ndim == 1:
            matrix = matrix.reshape(1, -1)

        ##  Return matrix on GF(2)
        return matrix

    ##  SECTION: ----Create Hermitian operator from occupy----
    @staticmethod
    @abstractmethod
    def HermitianOperatorFromOccupy(occupy_x, occupy_z) -> 'Operator':
        """
        input.occupy_x：Which sites have X operator.
        input.occupy_z：Which sites have Z operator.
        output：Same type new object, hermitian operator with given occupy.
        """

        pass

    ##  SECTION: ----Create Hermitian operator from vector----
    @staticmethod
    @abstractmethod
    def HermitianOperatorFromVector(vector) -> 'Operator':
        """
        input.vector：GF(2) vector (2n length)
        output：Same type new object, hermitian operator with given vector.
        """

        pass

    ##  SECTION: ----Check commute----
    @staticmethod
    @abstractmethod
    def commute(A, B) -> bool:
        """
        input.A：Operator A
        input.B：Operator B
        output：bool, True if A and B commute.
        """

        pass
