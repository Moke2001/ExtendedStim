"""
Module: Frame
Defines the Pauli Frame model for tracking error propagation in Clifford circuit simulation.
It is a simplified/faster version of `Platform`.
"""
import galois
import numpy as np
from extendedstim.Physics.MajoranaOperator import MajoranaOperator
from extendedstim.Physics.PauliOperator import PauliOperator
from extendedstim.tools.TypingTools import isinteger


class Frame:
    GF = galois.GF(2)

    ##  SECTION: ----Constructor----
    def __init__(self):
        """
        Initialize an empty frame with 0 qubits and 0 fermionic sites.
        output: None.
        influence: Initialize an empty frame.
        """
        self.pauli_number = 0
        self.majorana_number = 0
        self.frame = None
        self.pauli_frame = None
        self.majorana_frame = None

    ##  SECTION: ----Initialization----
    def trap(self, majorana_number: int, pauli_number: int):
        """
        input.majorana_number: Number of fermionic sites.
        input.pauli_number: Number of qubits.
        influence: Randomly initialize the frame based on stabilizers (simulating projection to random +1/-1 eigenstates).
        """
        self.pauli_number = pauli_number
        self.majorana_number = majorana_number
        self.frame = self.GF.Zeros(2 * pauli_number + 2 * majorana_number)
        for i in range(majorana_number):
            if np.random.rand() > 0.5:
                temp = self.GF(np.zeros(2 * majorana_number, dtype=int))
                temp[i * 2] = 1
                temp[i * 2 + 1] = 1
                self.frame[0:self.majorana_number * 2] += temp
        for i in range(pauli_number):
            if np.random.rand() > 0.5:
                temp = self.GF(np.zeros(2 * pauli_number, dtype=int))
                temp[i * 2 + 1] = 1
                self.frame[self.majorana_number * 2:] += temp

    ##  SECTION: ----Measurement----
    def measure(self, op: MajoranaOperator | PauliOperator, reference_value: int) -> int:
        """
        input.op: Hermitian operator to measure.
        input.reference_value: Noiseless reference measurement result (+1 or -1).
        output: int, Actual measurement result (+1 or -1).
        influence: Flip reference value based on commutation/anticommutation with current frame, and update frame with 50% probability.
        """

        ##  Check inputs
        assert op.is_hermitian, "Input operator must be Hermitian."
        assert reference_value == 1 or reference_value == -1, "Reference value must be +1 or -1."

        ##  Computing final state
        if isinstance(op, MajoranaOperator):
            v0 = op.get_vector(self.majorana_number)
            v1 = self.frame[0:self.majorana_number * 2]
            overlap_number = np.dot(v0, v1)
            weight_0 = np.dot(v0, v0)
            weight_1 = np.dot(v1, v1)
            if overlap_number + weight_0 * weight_1 == 0:
                result = reference_value
            else:
                result = -reference_value
            if np.random.rand() > 0.5:
                self.frame[0:self.majorana_number * 2] += v0
            return result
        elif isinstance(op, PauliOperator):
            v0 = op.get_vector(self.pauli_number)
            v0_x = v0[0::2]
            v0_z = v0[1::2]
            v1 = self.frame[self.majorana_number * 2:]
            v1_x = v1[0::2]
            v1_z = v1[1::2]
            overlap_number = np.dot(v0_x, v1_z) + np.dot(v0_z, v1_x)
            if overlap_number == 0:
                result = reference_value
            else:
                result = -reference_value
            if np.random.rand() > 0.5:
                self.frame[self.majorana_number * 2:] += v0
            return result
        else:
            raise ValueError

    ##  SECTION: ----Details----
    def x(self, qubit_index: int):
        """
        input.qubit_index: Target qubit.
        influence: Ideal single-qubit gate does not change error frame, so no operation.
        """
        pass

    ##  SECTION: ----Details----
    def y(self, qubit_index: int):
        """
        input.qubit_index: Target qubit.
        influence: Ideal gate, no operation.
        """
        pass

    ##  SECTION: ----Details----
    def z(self, qubit_index: int):
        """
        input.qubit_index: Target qubit.
        influence: Ideal gate, no operation.
        """
        pass

    ##  SECTION: ----Details----
    def h(self, qubit_index: int):
        """
        input.qubit_index: Target qubit.
        influence: Swap X and Z components of the corresponding qubit in Pauli frame.
        """

        ##  Check inputs
        assert isinteger(qubit_index) and 0 <= qubit_index < self.pauli_number
        qubit_index_x = self.majorana_number * 2 + qubit_index * 2
        qubit_index_z = self.majorana_number * 2 + qubit_index * 2 + 1

        ##  Computing final state
        cache = self.frame[qubit_index_x].copy()
        self.frame[qubit_index_x] = self.frame[qubit_index_z]
        self.frame[qubit_index_z] = cache

    ##  SECTION: ----S gate----
    def s(self, pauli_index: int):
        """
        input.pauli_index: Target qubit.
        influence: Add X component to Z component in Pauli frame (Z -> ZX).
        """

        ##  Check inputs
        assert isinteger(pauli_index) and 0 <= pauli_index < self.pauli_number
        qubit_index_x = self.majorana_number * 2 + pauli_index * 2
        qubit_index_z = self.majorana_number * 2 + pauli_index * 2 + 1

        ##  Computing final state
        self.frame[qubit_index_z] += self.frame[qubit_index_x]

    ##  SECTION: ----Details----
    def u(self, majorana_index: int):
        """
        input.majorana_index: Target fermionic site.
        influence: Ideal gate, no operation.
        """
        pass

    ##  SECTION: ----Details----
    def v(self, majorana_index: int):
        """
        input.majorana_index: Target fermionic site.
        influence: Ideal gate, no operation.
        """
        pass

    ##  SECTION: ----Details----
    def n(self, majorana_index: int):
        """
        input.majorana_index: Target fermionic site.
        influence: Ideal gate, no operation.
        """
        pass

    ##  SECTION: ----Details----
    def p(self, majorana_index: int):
        """
        input.majorana_index: Target fermionic site.
        influence: Swap X and Z components of corresponding fermionic site in frame.
        """

        ##  Check inputs
        assert isinteger(majorana_index) and 0 <= majorana_index < self.majorana_number
        majorana_index_x = majorana_index * 2
        majorana_index_z = majorana_index * 2 + 1

        ##  Computing fianl state
        cache = self.frame[majorana_index_x].copy()
        self.frame[majorana_index_x] = self.frame[majorana_index_z]
        self.frame[majorana_index_z] = cache

    ##  SECTION: ----Control X gate----
    def cx(self, control_index: int, target_index: int):
        """
        input.control_index: Control qubit index.
        input.target_index: Target qubit index.
        influence: Update Pauli frame: target_X += control_X, control_Z += target_Z.
        """

        ##  Check inputs
        assert isinteger(control_index) and 0 <= control_index < self.pauli_number
        assert isinteger(target_index) and 0 <= target_index < self.pauli_number
        control_qubit_index_x = self.majorana_number * 2 + control_index * 2
        control_qubit_index_z = self.majorana_number * 2 + control_index * 2 + 1
        target_qubit_index_x = self.majorana_number * 2 + target_index * 2
        target_qubit_index_z = self.majorana_number * 2 + target_index * 2 + 1

        ##  Computing final state
        self.frame[target_qubit_index_x] += self.frame[control_qubit_index_x]
        self.frame[control_qubit_index_z] += self.frame[target_qubit_index_z]


    ##  SECTION: ----Fermionic control X gate----
    def cnx(self, control_index: int, target_index: int):
        """
        input.control_index: Fermionic control site index.
        input.target_index: Target qubit index.
        influence: Update fermion-qubit frame based on mixed gate rules.
        """

        ##  Check inputs
        assert isinteger(control_index) and 0 <= control_index < self.majorana_number
        assert isinteger(target_index) and 0 <= target_index < self.pauli_number
        control_majorana_index_x = control_index * 2
        control_majorana_index_z = control_index * 2 + 1
        target_qubit_index_x = self.majorana_number * 2 + target_index * 2
        target_qubit_index_z = self.majorana_number * 2 + target_index * 2 + 1

        ##  Computing final state
        target_x = self.frame[target_qubit_index_x] + self.frame[control_majorana_index_x] + self.frame[control_majorana_index_z]
        control_x = self.frame[control_majorana_index_x] + self.frame[target_qubit_index_z]
        control_z = self.frame[control_majorana_index_z] + self.frame[target_qubit_index_z]
        self.frame[target_qubit_index_x] = target_x
        self.frame[control_majorana_index_x] = control_x
        self.frame[control_majorana_index_z] = control_z


    ##  SECTION: ----Fermionic control N gate----
    def cnn(self, control_index: int, target_index: int):
        """
        input.control_index: Control fermionic site index.
        input.target_index: Target fermionic site index.
        influence: Update frame based on fermionic gate rules.
        """

        ##  Check inputs
        assert isinteger(control_index) and 0 <= control_index < self.majorana_number
        assert isinteger(target_index) and 0 <= target_index < self.majorana_number
        control_majorana_index_x = control_index * 2
        control_majorana_index_z = control_index * 2 + 1
        target_majorana_index_x = target_index * 2
        target_majorana_index_z = target_index * 2 + 1

        ##  Computing final state
        control_x = self.frame[control_majorana_index_x] + self.frame[target_majorana_index_x] + self.frame[target_majorana_index_z]
        control_z = self.frame[control_majorana_index_z] + self.frame[target_majorana_index_x] + self.frame[target_majorana_index_z]
        target_x = self.frame[target_majorana_index_x] + self.frame[control_majorana_index_x] + self.frame[control_majorana_index_z]
        target_z = self.frame[target_majorana_index_z] + self.frame[control_majorana_index_x] + self.frame[control_majorana_index_z]
        self.frame[control_majorana_index_x] = control_x
        self.frame[control_majorana_index_z] = control_z
        self.frame[target_majorana_index_x] = target_x
        self.frame[target_majorana_index_z] = target_z


    ##  SECTION: ----Fermionic braid gate----
    def braid(self, control_index: int, target_index: int, *args):
        """
        input.control_index: Fermionic site index.
        input.target_index: Fermionic site index.
        influence: Swap specific components in fermionic frame.
        """

        ##  Check inputs
        assert isinteger(control_index) and 0 <= control_index < self.majorana_number
        assert isinteger(target_index) and 0 <= target_index < self.majorana_number
        control_majorana_index_z = control_index * 2 + 1
        target_majorana_index_x = target_index * 2

        ##  Computing final state
        cache = self.frame[control_majorana_index_z].copy()
        self.frame[control_majorana_index_z] = self.frame[target_majorana_index_x]
        self.frame[target_majorana_index_x] = cache


    ##  SECTION: ----Error gate----
    def error(self, p: float, error_pattern: str, index_0: int, index_1: int = None):
        """
        input.p: Error probability.
        input.error_pattern: Error pattern string (e.g. 'X', 'Y', 'Z'...).
        input.index_0: Target index 0.
        input.index_1: Target index 1.
        influence: Apply an error with probability p.
        """

        ##  Compute
        if np.random.rand() < p:
            string_0 = error_pattern[0]
            string_1 = error_pattern[1]
            if string_0 == 'X':
                self.frame[index_0 * 2 + self.majorana_number * 2] += self.GF(1)
            elif string_0 == 'Y':
                self.frame[index_0 * 2 + self.majorana_number * 2] += self.GF(1)
                self.frame[index_0 * 2 + self.majorana_number * 2 + 1] += self.GF(1)
            elif string_0 == 'Z':
                self.frame[index_0 * 2 + self.majorana_number * 2 + 1] += self.GF(1)
            elif string_0 == 'U':
                self.frame[index_0 * 2] += self.GF(1)
            elif string_0 == 'V':
                self.frame[index_0 * 2 + 1] += self.GF(1)
            elif string_0 == 'N':
                self.frame[index_0 * 2] += self.GF(1)
                self.frame[index_0 * 2 + 1] += self.GF(1)
            else:
                raise ValueError(f"string_0={string_0} is not valid.")
            if string_1 == 'X':
                self.frame[index_1 * 2 + self.majorana_number * 2] += self.GF(1)
            elif string_1 == 'Y':
                self.frame[index_1 * 2 + self.majorana_number * 2] += self.GF(1)
                self.frame[index_1 * 2 + self.majorana_number * 2 + 1] += self.GF(1)
            elif string_1 == 'Z':
                self.frame[index_1 * 2 + self.majorana_number * 2 + 1] += self.GF(1)
            elif string_1 == 'U':
                self.frame[index_1 * 2] += self.GF(1)
            elif string_1 == 'V':
                self.frame[index_1 * 2 + 1] += self.GF(1)
            elif string_1 == 'N':
                self.frame[index_1 * 2] += self.GF(1)
                self.frame[index_1 * 2 + 1] += self.GF(1)
            elif string_1 == '_':
                pass
            else:
                raise ValueError(f"string_0={string_0} or string_1={string_1} is not valid.")


    ##  SECTION: ----Reset gate----
    def reset(self, pauli_index: int):
        """
        input.pauli_index: Target qubit.
        influence: Clear corresponding qubit components in frame and randomly introduce a Z error (simulating reset after measurement).
        """

        ##  Check inputs
        assert isinteger(pauli_index) and 0 <= pauli_index < self.pauli_number
        qubit_index_x = self.majorana_number * 2 + pauli_index * 2
        qubit_index_z = self.majorana_number * 2 + pauli_index * 2 + 1

        ##  Compute
        self.frame[qubit_index_x] = 0
        self.frame[qubit_index_z] = 0
        if np.random.rand() < 0.5:
            self.frame[qubit_index_z] = self.GF(1)
