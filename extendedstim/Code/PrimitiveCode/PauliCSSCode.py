"""
Module: PauliCode
A class for Pauli CSS code, each of which has X and Z stabilizer generators and physical number of modes as properties.
"""
import numpy as np
from extendedstim.Code.PrimitiveCode.PauliCode import PauliCode
from extendedstim.Code.PrimitiveCode.QuantumCSSCode import QuantumCSSCode
from extendedstim.Physics.PauliOperator import PauliOperator
from extendedstim.tools.GaloisTools import orthogonalize, mip_distance_caculator


class PauliCSSCode(PauliCode, QuantumCSSCode):
    """
    An example:
    $S_0=i\\hat X_0\\hat X_1\\hat X_2\\hat X_3$
    $S_1=i\\hat X_1\\hat X_2\\hat X_5\\hat X_6$
    $S_2=i\\hat X_0\\hat X_2\\hat X_4\\hat X_6$
    $S_3=i\\hat Z_0\\hat Z_1\\hat Z_2\\hat Z_3$
    $S_4=i\\hat Z_1\\hat Z_2\\hat Z_5\\hat Z_6$
    $S_5=i\\hat Z_0\\hat Z_2\\hat Z_4\\hat Z_6$
    n=7
    """

    ##  SECTION: ----Constructor----
    def __init__(self, generators_x, generators_z, physical_number: int)->None:
        """
        input.generators_x：list of PauliOperator, X stabilizers
        input.generators_z：list of PauliOperator, Z stabilizers
        input.physical_number：number of physical qubits
        """

        ##  Check inputs
        try:
            generators_x: list[PauliOperator]=[generator if isinstance(generator, PauliOperator) else PauliOperator(generator) for generator in generators_x]
            generators_z: list[PauliOperator]=[generator if isinstance(generator, PauliOperator) else PauliOperator(generator) for generator in generators_z]
        except ValueError as e:
            raise ValueError("generators_x must be a list of PauliOperator.") from e

        try:
            physical_number=int(physical_number)
            assert physical_number>0, "physical_number must be a positive integer."
        except Exception as e:
            raise ValueError("physical_number must be a positive integer.") from e
        
        ##  Initialize parent classes
        QuantumCSSCode.__init__(self, generators_x, generators_z, physical_number)


    ##  SECTION: ----Get X code distance----
    @property
    def distance_x(self) -> int:
        """
        output：int, code distance by mip algorithm for X stabilizers
        """
        return mip_distance_caculator(self.check_matrix_x, self.logical_operators_x)


    ##  SECTION: ----Get Z code distance----
    @property
    def distance_z(self) -> int:
        """
        output：int, code distance by mip algorithm for Z stabilizers
        """
        return mip_distance_caculator(self.check_matrix_z, self.logical_operators_z)


    ##  SECTION: ----Get X logical operators----
    @property
    def logical_operators_x(self):
        if self._logical_operators_x is not None:
            return self._logical_operators_x
        else:
            matrix = self.check_matrix_x
            codewords = matrix.null_space()
            independent_null_basis_list = []
            for vec in codewords:
                rank_before = np.linalg.matrix_rank(matrix)
                matrix = np.vstack([matrix, vec])
                if np.linalg.matrix_rank(matrix) == rank_before + 1:
                    independent_null_basis_list.append(vec)
            basis_list = orthogonalize(independent_null_basis_list)
            pauli_logical_operators_x = []
            pauli_logical_operators_z = []
            for i in range(len(basis_list)):
                occupy_temp = np.where(basis_list[i] == 1)[0]
                temp = PauliOperator.HermitianOperatorFromOccupy(occupy_temp, [])
                pauli_logical_operators_x.append(temp)
                temp = PauliOperator.HermitianOperatorFromOccupy([], occupy_temp)
                pauli_logical_operators_z.append(temp)
            self._logical_operators_x = np.array(pauli_logical_operators_x, dtype=PauliOperator)
            self._logical_operators_z = np.array(pauli_logical_operators_z, dtype=PauliOperator)
            return self._logical_operators_x
        

    ##  SECTION: ----get Z logical operators----
    @property
    def logical_operators_z(self):
        if self._logical_operators_z is not None:
            return self._logical_operators_z
        else:
            _=self.logical_operators_x
            return self._logical_operators_z


    # SECTION: ----get \\hat gamma^\\prime logical operators----
    @staticmethod
    def FromCheckMatrix(H_x,H_z):
        """
        input.H_x: np.ndarray of GF(2), check matrix of X stabilizers
        input.H_z: np.ndarray of GF(2), check matrix of Z stabilizers
        output：PauliCSSCode, code generated from check matrices
        """

        ##  Check input
        try:
            H_x=np.array(H_x, dtype=int)
            H_z=np.array(H_z, dtype=int)
            assert H_x.shape[1]==H_z.shape[1], "H_x and H_z must have the same number of columns."
        except Exception as e:
            raise ValueError("H_x and H_z must be np.ndarray of GF(2).") from e
        
        ##  Generate generators
        generators_x = []
        generators_z = []
        for i in range(len(H_x)):
            generators_x.append(PauliOperator.HermitianOperatorFromOccupy(np.where(H_x[i]==1)[0], []))
        for i in range(len(H_z)):
            generators_z.append(PauliOperator.HermitianOperatorFromOccupy([], np.where(H_z[i]==1)[0]))
        physical_number: int=H_x.shape[1]

        ##  Initialize code
        return PauliCSSCode(generators_x, generators_z, physical_number)
