"""
Module: Platform
Defines platform state evolution and measurement model based on stabilizer formalism, supporting Pauli and Majorana gates, noise, measurement, and reset.
"""
import galois
import numpy as np
from galois import GF2
from extendedstim.Physics.MajoranaOperator import MajoranaOperator
from extendedstim.Physics.PauliOperator import PauliOperator
from extendedstim.tools.GaloisTools import solve
from extendedstim.tools.TypingTools import isinteger


class Platform:
    GF = galois.GF(2)

    ##  SECTION: ----Constructor----
    def __init__(self):
        """
        Initialize an empty platform.
        output: None.
        influence: Initialize empty internal state.
        """
        self.pauli_number = 0
        self.majorana_number = 0
        self.stabilizers_pauli = []
        self.stabilizers_majorana = []
        self.stabilizers = None
        self.coffs = None


    ##  SECTION: ----Initialization----
    def trap(self, majorana_number: int, pauli_number: int):
        """
        input.majorana_number: Number of fermionic sites.
        input.pauli_number: Number of qubits.
        influence: Initialize stabilizers for the given number of sites/qubits.
        """

        ##  Check inputs
        assert isinteger(majorana_number) and majorana_number >= 0
        assert isinteger(pauli_number) and pauli_number >= 0

        ##  qubitsfermionic sites
        self.pauli_number = pauli_number
        self.majorana_number = majorana_number
        self.stabilizers = self.GF(np.zeros((majorana_number + pauli_number, 2 * majorana_number + 2 * pauli_number), dtype=int))
        self.coffs = np.ones(majorana_number + pauli_number, dtype=complex)

        ##  0
        for i in range(majorana_number):
            self.stabilizers[i, 2 * i] = 1
            self.stabilizers[i, 2 * i + 1] = 1
            self.coffs[i] = 1j
        for i in range(pauli_number):
            self.stabilizers[i + majorana_number, 2 * i + 1 + 2 * majorana_number] = 1
            self.coffs[i + majorana_number] = 1

        ##  SECTION: ----Measurement----
    def measure(self, op: MajoranaOperator | PauliOperator) -> int:
        """
        input.op: PauliOperator or MajoranaOperator (must be Hermitian).
        output: int, +1 or -1 measurement result.
        influence: Update stabilizer group or check consistency.
        """

        ##  Check inputs
        assert op.is_hermitian
        if isinstance(op, MajoranaOperator):
            vector_op = np.append(op.get_vector(self.majorana_number), self.GF.Zeros(self.pauli_number * 2))
        else:
            vector_op = np.append(self.GF.Zeros(self.majorana_number * 2), op.get_vector(self.pauli_number))

        ##  Measurement
        first_index = None
        for i in range(len(self.stabilizers)):
            pauli_commute = (np.dot(self.stabilizers[i][self.majorana_number * 2::2], vector_op[self.majorana_number * 2 + 1::2]) +
                             np.dot(self.stabilizers[i][self.majorana_number * 2 + 1::2], vector_op[self.majorana_number * 2::2]))
            majorana_commute = (np.dot(self.stabilizers[i][0:self.majorana_number * 2], vector_op[0:self.majorana_number * 2]) +
                                np.sum(vector_op[0:self.majorana_number * 2]) * np.sum(self.stabilizers[i][0:self.majorana_number * 2]))

            if majorana_commute + pauli_commute == 1 and first_index is None:
                first_index = i
            elif majorana_commute + pauli_commute == 1 and first_index is not None:
                factor = np.sum([np.sum(self.stabilizers[i][temp + 1:2 * self.majorana_number]) for temp in range(2 * self.majorana_number) if
                                 self.stabilizers[first_index][temp] == 1])
                self.stabilizers[i] += self.stabilizers[first_index]
                if np.mod(factor, 2) == 0:
                    factor = 1
                else:
                    factor = -1
                self.coffs[i] = self.coffs[i] * self.coffs[first_index] * factor
            else:
                pass

        ##  Measurement
        if first_index is not None:
            if np.random.rand() < 0.5:
                self.stabilizers[first_index] = vector_op
                self.coffs[first_index] = op.coff
                return 1
            else:
                self.stabilizers[first_index] = vector_op
                self.coffs[first_index] = -op.coff
                return -1

        ##  Measurement
        else:
            solution = solve(self.stabilizers, vector_op)
            coff = np.prod([self.coffs[i] for i in range(len(self.stabilizers)) if solution[i] == 1])
            vector_now, factor = majorana_factor(self.majorana_number,
                                                 [self.stabilizers[i] for i in range(len(self.stabilizers)) if solution[i] == 1])
            coff = coff * factor
            if coff == op.coff:
                return 1
            else:
                return -1


    ##  SECTION: ----X gate----
    def x(self, qubit_index: int):
        """
        input.qubit_index: Target qubit.
        influence: Flip phase of stabilizers overlapping with Z.
        """

        ##  Check inputs
        assert isinteger(qubit_index) and 0 <= qubit_index < self.pauli_number

        ##  Compute final state
        indices = np.where(self.stabilizers[:, self.majorana_number * 2 + qubit_index * 2 + 1] == 1)[0]
        self.coffs[indices] = -self.coffs[indices]


        ##  SECTION: ----Y gate----
    def y(self, qubit_index: int):
        """
        input.qubit_index: Target qubit.
        influence: Flip phase of stabilizers overlapping with both X and Z.
        """

        ##  Check inputs
        assert isinteger(qubit_index) and 0 <= qubit_index < self.pauli_number

        ##  Compute final state
        indices_x = np.where(self.stabilizers[:, self.majorana_number * 2 + qubit_index * 2 + 1] == 1)[0]
        self.coffs[indices_x] = -self.coffs[indices_x]
        indices_z = np.where(self.stabilizers[:, self.majorana_number * 2 + qubit_index * 2] == 1)[0]
        self.coffs[indices_z] = -self.coffs[indices_z]


        ##  SECTION: ----Z gate----
    def z(self, qubit_index: int):
        """
        input.qubit_index: Target qubit.
        influence: Flip phase of stabilizers overlapping with X.
        """

        ##  Check inputs
        assert isinteger(qubit_index) and 0 <= qubit_index < self.pauli_number

        ##  Compute final state
        indices = np.where(self.stabilizers[:, self.majorana_number * 2 + qubit_index * 2] == 1)[0]
        self.coffs[indices] = -self.coffs[indices]


    ##  SECTION: ----H gate----
    def h(self, qubit_index: int):
        """
        input.qubit_index: Target qubit.
        influence: Swap X/Z support and update phase.
        """

        ##  Check inputs
        assert isinteger(qubit_index) and 0 <= qubit_index < self.pauli_number

        ##  Compute final state
        indices = np.where(np.logical_and(self.stabilizers[:, self.majorana_number * 2 + qubit_index * 2 + 1] == 1,
                                          self.stabilizers[:, self.majorana_number * 2 + qubit_index * 2] == 1))[0]
        self.coffs[indices] = -self.coffs[indices]

        caches = self.stabilizers[:, self.majorana_number * 2 + qubit_index * 2].copy()
        self.stabilizers[:, self.majorana_number * 2 + qubit_index * 2] = self.stabilizers[:, self.majorana_number * 2 + qubit_index * 2 + 1]
        self.stabilizers[:, self.majorana_number * 2 + qubit_index * 2 + 1] = caches


    ##  SECTION: ----S gate----
    def s(self, pauli_index: int):
        """
        input.pauli_index: Target qubit.
        influence: Z += X (Phase gate).
        """

        ##  Check inputs
        assert isinteger(pauli_index) and 0 <= pauli_index < self.pauli_number

        ##  Compute final state
        indices = np.where(self.stabilizers[:, self.majorana_number * 2 + pauli_index * 2] == 1)[0]
        self.coffs[indices] = 1j * self.coffs[indices]
        self.stabilizers[:, self.majorana_number * 2 + pauli_index * 2 + 1] += self.stabilizers[:, self.majorana_number * 2 + pauli_index * 2]


    ##  SECTION: ----U gate----
    def u(self, majorana_index: int):
        """
        input.majorana_index: Target fermionic site.
        influence: Flip phase based on overlap weights.
        """

        ##  Check inputs
        assert isinteger(majorana_index) and 0 <= majorana_index < self.majorana_number

        ##  Compute final state
        weights = np.sum(self.stabilizers[:, 0:2 * self.majorana_number], axis=1)
        overlaps = self.GF(np.where(self.stabilizers[:, majorana_index * 2] == 1, 1, 0))
        indices = np.where(weights + overlaps == 0)[0]
        self.coffs[indices] = -self.coffs[indices]


    ##  SECTION: ----V gate----
    def v(self, majorana_index: int):
        """
        input.majorana_index: Target fermionic site.
        influence: Flip phase based on Z overlap.
        """

        ##  Check inputs
        assert isinteger(majorana_index) and 0 <= majorana_index < self.majorana_number

        ##  Compute final state
        weights = np.sum(self.stabilizers[:, 0:2 * self.majorana_number], axis=1)
        overlaps = self.GF(np.where(self.stabilizers[:, majorana_index * 2 + 1] == 1, 1, 0))
        indices = np.where(weights + overlaps == 0)[0]
        self.coffs[indices] = -self.coffs[indices]


    ##  SECTION: ----N gate----
    def n(self, majorana_index: int):
        """
        input.majorana_index: Target fermionic site.
        influence: Flip phase based on X or Z parity.
        """

        ##  Check inputs
        assert isinteger(majorana_index) and 0 <= majorana_index < self.majorana_number

        ##  Compute final state
        weights = np.sum(self.stabilizers[:, 0:2 * self.majorana_number], axis=1)
        overlaps = self.GF(
            np.where(np.logical_xor(self.stabilizers[:, majorana_index * 2] == 1, self.stabilizers[:, majorana_index * 2 + 1] == 1), 1, 0))
        indices = np.where(weights + overlaps == 1)[0]
        self.coffs[indices] = -self.coffs[indices]


    ##  SECTION: ----P gate----
    def p(self, majorana_index: int):
        """
        input.majorana_index: Target fermionic site.
        influence: Swap X/Z support.
        """

        ##  Check inputs
        assert isinteger(majorana_index) and 0 <= majorana_index < self.majorana_number

        ##  Compute final state
        indices = np.where(np.logical_and(self.stabilizers[:, majorana_index * 2 + 1] == 1, self.stabilizers[:, majorana_index * 2] == 0))[0]
        self.coffs[indices] = -self.coffs[indices]
        caches = self.stabilizers[:, majorana_index * 2].copy()
        self.stabilizers[:, majorana_index * 2] = self.stabilizers[:, majorana_index * 2 + 1]
        self.stabilizers[:, majorana_index * 2 + 1] = caches


    ##  SECTION: ----CX gate----
    def cx(self, control_index: int, target_index: int):
        """
        input.control_index: Control qubit index.
        input.target_index: Target qubit index.
        influence: Linear transformation of stabilizers (X propagates to target, Z to control).
        """

        ##  Check inputs
        assert isinteger(control_index) and 0 <= control_index < self.pauli_number
        assert isinteger(target_index) and 0 <= target_index < self.pauli_number

        ##  Compute final state
        control_qubit_index_x = control_index * 2 + self.majorana_number * 2
        control_qubit_index_z = control_index * 2 + 1 + self.majorana_number * 2
        target_qubit_index_x = target_index * 2 + self.majorana_number * 2
        target_qubit_index_z = target_index * 2 + 1 + self.majorana_number * 2
        targets_x = self.stabilizers[:, control_qubit_index_x] + self.stabilizers[:, target_qubit_index_x]
        control_z = self.stabilizers[:, control_qubit_index_z] + self.stabilizers[:, target_qubit_index_z]
        self.stabilizers[:, target_qubit_index_x] = targets_x
        self.stabilizers[:, control_qubit_index_z] = control_z


    ##  SECTION: ----CNX gate----
    def cnx(self, control_index: int, target_index: int):
        """
        input.control_index: Fermionic control site index.
        input.target_index: Target qubit index.
        influence: Update stabilizers according to derived rules (see docs).
        """

        ##  Check inputs
        assert isinteger(control_index) and 0 <= control_index < self.majorana_number
        assert isinteger(target_index) and 0 <= target_index < self.pauli_number

        ##  Compute final state
        control_majorana_index_x = control_index * 2
        control_majorana_index_z = control_index * 2 + 1
        target_qubit_index_x = self.majorana_number * 2 + target_index * 2
        target_qubit_index_z = self.majorana_number * 2 + target_index * 2 + 1
        factor = [1j if self.stabilizers[temp, target_qubit_index_z] == 1 else 1 for temp in range(len(self.stabilizers))]
        indices = np.where(np.logical_and(self.stabilizers[:, target_qubit_index_z] == 1, self.stabilizers[:, control_majorana_index_z] == 1))[0]
        self.coffs = self.coffs * np.array(factor)
        self.coffs[indices] = -self.coffs[indices]
        targets_x = self.stabilizers[:, target_qubit_index_x] + self.stabilizers[:, control_majorana_index_x] + self.stabilizers[
            :, control_majorana_index_z]
        controls_x = self.stabilizers[:, control_majorana_index_x] + self.stabilizers[:, target_qubit_index_z]
        controls_z = self.stabilizers[:, control_majorana_index_z] + self.stabilizers[:, target_qubit_index_z]
        self.stabilizers[:, target_qubit_index_x] = targets_x
        self.stabilizers[:, control_majorana_index_x] = controls_x
        self.stabilizers[:, control_majorana_index_z] = controls_z

    ##  SECTION: ----CNN gate----
    def cnn(self, control_index: int, target_index: int):
        """
        input.control_index: Control fermionic site index.
        input.target_index: Target fermionic site index.
        influence: Update stabilizers according to derived rules (see docs).
        """

        ##  Check inputs
        assert isinteger(control_index) and 0 <= control_index < self.majorana_number
        assert isinteger(target_index) and 0 <= target_index < self.majorana_number
        if target_index < control_index:
            cache = target_index
            target_index = control_index
            control_index = cache

        ##  Compute final state
        control_majorana_index_x = control_index * 2
        control_majorana_index_z = control_index * 2 + 1
        target_majorana_index_x = target_index * 2
        target_majorana_index_z = target_index * 2 + 1
        temp = self.stabilizers[:, (control_majorana_index_x, control_majorana_index_z, target_majorana_index_x, target_majorana_index_z)]
        indices_1j = np.where(np.logical_or(np.logical_or(np.all(temp == GF2([1, 0, 0, 0]), axis=1), np.all(temp == GF2([0, 1, 0, 0]), axis=1)),
                                            np.logical_or(np.all(temp == GF2([0, 0, 1, 0]), axis=1), np.all(temp == GF2([0, 0, 0, 1]), axis=1))))[0]
        indices_minus_1 = np.where(np.logical_or(np.all(temp == GF2([1, 0, 0, 1]), axis=1), np.all(temp == GF2([0, 1, 1, 0]), axis=1)))[0]
        indices_minus_1j = np.where(np.logical_or(np.logical_or(np.all(temp == GF2([1, 1, 1, 0]), axis=1), np.all(temp == GF2([1, 1, 0, 1]), axis=1)),
                                                  np.logical_or(np.all(temp == GF2([1, 0, 1, 1]), axis=1),
                                                                np.all(temp == GF2([0, 1, 1, 1]), axis=1))))[0]
        self.coffs[indices_1j] = self.coffs[indices_1j] * 1j
        self.coffs[indices_minus_1] = -self.coffs[indices_minus_1]
        self.coffs[indices_minus_1j] = self.coffs[indices_minus_1j] * (-1j)
        control_x = self.stabilizers[:, control_majorana_index_x] + self.stabilizers[:, target_majorana_index_x] + self.stabilizers[
            :, target_majorana_index_z]
        control_z = self.stabilizers[:, control_majorana_index_z] + self.stabilizers[:, target_majorana_index_x] + self.stabilizers[
            :, target_majorana_index_z]
        target_x = self.stabilizers[:, target_majorana_index_x] + self.stabilizers[:, control_majorana_index_x] + self.stabilizers[
            :, control_majorana_index_z]
        target_z = self.stabilizers[:, target_majorana_index_z] + self.stabilizers[:, control_majorana_index_x] + self.stabilizers[
            :, control_majorana_index_z]
        self.stabilizers[:, control_majorana_index_x] = control_x
        self.stabilizers[:, control_majorana_index_z] = control_z
        self.stabilizers[:, target_majorana_index_x] = target_x
        self.stabilizers[:, target_majorana_index_z] = target_z


    ##  SECTION: ----BRAID gate----
    def braid(self, control_index: int, target_index: int):
        """
        input.control_index: Fermionic site index.
        input.target_index: Fermionic site index.
        influence: Swap specific support and flip phase.
        """

        ##  Check inputs
        assert isinteger(control_index) and 0 <= control_index < self.majorana_number
        assert isinteger(target_index) and 0 <= target_index < self.majorana_number
        control_majorana_index_z = control_index * 2 + 1
        target_majorana_index_x = target_index * 2
        if target_index > control_index:
            numbers_of_mid = np.sum(self.stabilizers[:, control_majorana_index_z + 1:target_majorana_index_x], axis=1)
        else:
            numbers_of_mid = np.sum(self.stabilizers[:, target_majorana_index_x + 1:control_majorana_index_z], axis=1)

        ##  Compute final state
        factors_0 = self.stabilizers[:, control_majorana_index_z] * numbers_of_mid
        factors_1 = self.stabilizers[:, target_majorana_index_x] * numbers_of_mid
        factors_2 = self.stabilizers[:, control_majorana_index_z] * self.stabilizers[:, target_majorana_index_x]
        factors_4 = self.stabilizers[:, target_majorana_index_x]
        judge = factors_0 + factors_1 + factors_2 + factors_4
        factors = [-1 if judge[temp] == 1 else 1 for temp in range(len(self.stabilizers))]
        self.coffs = self.coffs * np.array(factors)
        caches = self.stabilizers[:, control_majorana_index_z].copy()
        self.stabilizers[:, control_majorana_index_z] = self.stabilizers[:, target_majorana_index_x]
        self.stabilizers[:, target_majorana_index_x] = caches


    ##  SECTION: ----ERROR gate----
    def error(self, p: float, error_pattern: str, index_0: int, index_1: int = None):
        """
        input.p: Error probability.
        input.error_pattern: Error pattern string.
        input.index_0: Target index 0.
        input.index_1: Target index 1.
        influence: Apply an error with probability p.
        """

        ##  Check inputs
        assert isinteger(index_0) and 0 <= index_0 < self.majorana_number
        if index_1 is not None:
            assert isinteger(index_1) and 0 <= index_1 < self.majorana_number

        if np.random.rand() < p:
            string_0 = error_pattern[0]
            string_1 = error_pattern[1]
            if string_0 == 'X':
                self.x(index_0)
            elif string_0 == 'Y':
                self.y(index_0)
            elif string_0 == 'Z':
                self.z(index_0)
            elif string_0 == 'U':
                self.u(index_0)
            elif string_0 == 'V':
                self.v(index_0)
            elif string_0 == 'N':
                self.n(index_0)
            else:
                raise ValueError(f"string_0={string_0} is not valid.")
            if string_1 == 'X':
                self.x(index_1)
            elif string_1 == 'Y':
                self.y(index_1)
            elif string_1 == 'Z':
                self.z(index_1)
            elif string_1 == 'U':
                self.u(index_1)
            elif string_1 == 'V':
                self.v(index_1)
            elif string_1 == 'N':
                self.n(index_1)
            elif string_1 == '_':
                pass
            else:
                raise ValueError(f"string_0={string_0} or string_1={string_1} is not valid.")


    ##  SECTION: ----RESET gate----
    def reset(self, pauli_index: int):
        """
        input.pauli_index: Target qubit.
        influence: Set corresponding stabilizer row to Z=+1.
        """

        ##  Check inputs
        assert isinteger(pauli_index) and 0 <= pauli_index < self.pauli_number

        ##  0
        vector_op = np.append(self.GF.Zeros(self.majorana_number * 2), PauliOperator([], [pauli_index], 1).get_vector(self.pauli_number))
        first_index = None
        for i in range(len(self.stabilizers)):
            pauli_commute = (np.dot(self.stabilizers[i][self.majorana_number * 2::2], vector_op[self.majorana_number * 2 + 1::2]) +
                             np.dot(self.stabilizers[i][self.majorana_number * 2 + 1::2], vector_op[self.majorana_number * 2::2]))
            majorana_commute = (np.dot(self.stabilizers[i][0:self.majorana_number * 2], vector_op[0:self.majorana_number * 2]) +
                                np.sum(vector_op[0:self.majorana_number * 2]) * np.sum(self.stabilizers[i][0:self.majorana_number * 2]))

            if majorana_commute + pauli_commute == 1 and first_index is None:
                first_index = i
            elif majorana_commute + pauli_commute == 1 and first_index is not None:
                factor = majorana_factor(self.majorana_number, [self.stabilizers[i], self.stabilizers[first_index]])
                self.stabilizers[i] += self.stabilizers[first_index]
                if np.mod(factor, 2) == 0:
                    factor = 1
                else:
                    factor = -1
                self.coffs[i] = self.coffs[i] * self.coffs[first_index] * factor
            else:
                pass
        if first_index is not None:
            self.stabilizers[first_index] = vector_op
            self.coffs[first_index] = 1
        else:
            solution = solve(self.stabilizers, vector_op)
            coff = np.prod([self.coffs[i] for i in range(len(self.stabilizers)) if solution[i] == 1])
            vector_now, factor = majorana_factor(self.majorana_number,
                                                 [self.stabilizers[i] for i in range(len(self.stabilizers)) if solution[i] == 1])
            coff = coff * factor
            if coff == -1:
                self.x(pauli_index)

    ##  SECTION: ----DETECTOR----
    def detect(self, op: MajoranaOperator | PauliOperator) -> int | None:
        """
        input.op: MajoranaOperator or PauliOperator.
        output: int or None, +1/-1 if deterministic, None if random.
        influence: None.
        """

        ##  Check inputs
        assert op.is_hermitian
        if isinstance(op, MajoranaOperator):
            vector_op = np.append(op.get_vector(self.majorana_number), self.GF.Zeros(self.pauli_number * 2))
        else:
            vector_op = np.append(self.GF.Zeros(self.majorana_number * 2), op.get_vector(self.pauli_number))

        ##  Compute final state
        solution = solve(self.stabilizers, vector_op)
        if solution is None:
            return None
        coff = np.prod([self.coffs[i] for i in range(len(self.stabilizers)) if solution[i] == 1])
        vector_now, factor = majorana_factor(self.majorana_number,
                                             [vector_op] + [self.stabilizers[i] for i in range(len(self.stabilizers)) if solution[i] == 1])
        coff = coff * factor
        if coff == op.coff:
            return 1
        else:
            return -1


##  CHAPTER: ====Majorana factor====
def majorana_factor(majorana_number: int, args: list) -> tuple[np.ndarray, int]:
    """
    input.majorana_number: Number of fermionic sites.
    input.args: List of vectors.
    output: tuple, (resultant vector, phase factor).
    influence: Calculate phase factor from Majorana operator multiplication.
    """
    first_vector = args[0]
    vector_now = first_vector
    factor = 0
    for vector in args[1:]:
        factor = factor + np.sum(
            [np.count_nonzero(vector_now[temp + 1:2 * majorana_number]) for temp in range(2 * majorana_number) if vector[temp] == 1])
        vector_now = vector_now + vector
    if np.mod(factor, 2) == 0:
        factor = 1
    elif np.mod(factor, 2) == 1:
        factor = -1
    return vector_now, factor
