import numpy as np
from Math.BinaryArray import BinaryArray as ba
from Physics.MajoranaOperator import MajoranaOperator
from Physics.PauliOperator import PauliOperator


class Platform:

    # %%  USER：构造方法
    def __init__(self):
        self.pauli_number = 0
        self.majorana_number = 0
        self.stabilizers_pauli = []
        self.stabilizers_majorana = []

    # %%  USER：对象方法
    ##  USER：初始化平台，定义fermionic sites和qubits数目
    def initialize(self, majorana_number, pauli_number, *args):
        self.pauli_number = pauli_number
        self.majorana_number = majorana_number
        if len(args) == 0:
            for i in range(majorana_number):
                self.stabilizers_majorana.append(MajoranaOperator([i], [i], 1j))
                self.stabilizers_pauli.append(PauliOperator([], [], 1))
            for i in range(pauli_number):
                self.stabilizers_majorana.append(MajoranaOperator([], [], 1))
                self.stabilizers_pauli.append(PauliOperator([], [i], 1))
        elif len(args) == 2:
            self.stabilizers_majorana = args[0]
            self.stabilizers_pauli = args[1]

    ##  USER：测量算符op，返回测量结果
    def measure(self, op):
        if isinstance(op, PauliOperator):
            stabilizers_now = self.stabilizers_pauli
        elif isinstance(op, MajoranaOperator):
            stabilizers_now = self.stabilizers_majorana
        else:
            raise NotImplementedError
        assert op.is_hermitian

        first_pauli = None
        first_index = -1
        for i in range(len(stabilizers_now)):
            if isinstance(op, PauliOperator):
                commute_flag = PauliOperator.commute(op, stabilizers_now[i])
            elif isinstance(op, MajoranaOperator):
                commute_flag = MajoranaOperator.commute(op, stabilizers_now[i])
            else:
                raise NotImplementedError

            if not commute_flag:
                if first_index == -1:
                    first_pauli = stabilizers_now[i]
                    first_index = i
                else:
                    stabilizers_now[i] = stabilizers_now[i] @ first_pauli
        if first_index == -1:
            matrix_pauli = PauliOperator.get_matrix(self.stabilizers_pauli, self.pauli_number)
            matrix_majorana = MajoranaOperator.get_matrix(self.stabilizers_majorana, self.majorana_number)
            matrix = ba.hstack(matrix_majorana, matrix_pauli)
            if isinstance(op, MajoranaOperator):
                vector_majorana = op.get_vector(self.majorana_number)
                vector_pauli = ba.zeros(self.pauli_number * 2)
                vector = ba.hstack(vector_majorana, vector_pauli)
            else:
                vector_majorana = ba.zeros(self.majorana_number * 2)
                vector_pauli = op.get_vector(self.pauli_number)
                vector = ba.hstack(vector_majorana, vector_pauli)
            result = ba.solve(matrix, vector)
            op_mul_pauli = PauliOperator([], [], 1)
            op_mul_majorana = MajoranaOperator([], [], 1)
            for i in range(len(result)):
                if result[i] == 1:
                    op_mul_majorana = op_mul_majorana @ self.stabilizers_majorana[i]
                    op_mul_pauli = op_mul_pauli @ self.stabilizers_pauli[i]
            coff = op_mul_pauli.coff * op_mul_majorana.coff
            if coff == op.coff:
                return 1
            else:
                return -1
        else:
            if np.random.rand() < 0.5:
                stabilizers_now[first_index] = op.copy()
                return 1
            else:
                stabilizers_now[first_index] = -op.copy()
                return -1

    ##  USER：X门，作用于qubit_index
    def X(self, qubit_index: int):
        for i in range(len(self.stabilizers_pauli)):
            if self.stabilizers_pauli[i].is_exist_occupy_z(qubit_index) is not None:
                self.stabilizers_pauli[i].coff = -self.stabilizers_pauli[i].coff

    ##  USER：Y门，作用于qubit_index
    def Y(self, qubit_index: int):
        for i in range(len(self.stabilizers_pauli)):
            if self.stabilizers_pauli[i].is_exist_occupy_x(qubit_index) is not None:
                self.stabilizers_pauli[i].coff = -self.stabilizers_pauli[i].coff
            if self.stabilizers_pauli[i].is_exist_occupy_z(qubit_index) is not None:
                self.stabilizers_pauli[i].coff = -self.stabilizers_pauli[i].coff

    ##  USER：Z门，作用于qubit_index
    def Z(self, qubit_index: int):
        for i in range(len(self.stabilizers_pauli)):
            if self.stabilizers_pauli[i].is_exist_occupy_x(qubit_index) is not None:
                self.stabilizers_pauli[i].coff = -self.stabilizers_pauli[i].coff

    ##  USER：Hadamard gate，作用于qubit_index
    def H(self, qubit_index: int):
        for i, stabilizer in enumerate(self.stabilizers_pauli):
            assert isinstance(stabilizer, PauliOperator)
            left, middle, right = stabilizer.split(qubit_index)
            if len(middle.occupy_x) == 0 and len(middle.occupy_z) == 0:
                continue
            elif len(middle.occupy_x) == 1 and len(middle.occupy_z) == 1:
                continue
            elif len(middle.occupy_x) == 0 and len(middle.occupy_z) == 1:
                middle.occupy_x = [qubit_index]
                middle.occupy_z = []
            elif len(middle.occupy_x) == 1 and len(middle.occupy_z) == 0:
                middle.occupy_x = []
                middle.occupy_z = [qubit_index]
            else:
                raise ValueError
            self.stabilizers_pauli[i] = left @ middle @ right

    ##  USER：gamma门，作用于majorana_index
    def U(self, majorana_index: int):
        op = MajoranaOperator([majorana_index], [], 1)
        for i, stabilizer in enumerate(self.stabilizers_majorana):
            if not MajoranaOperator.commute(stabilizer, op):
                stabilizer.coff = -stabilizer.coff

    ##  USER：gamma_prime门，作用于majorana_index
    def V(self, majorana_index: int):
        op = MajoranaOperator([], [majorana_index], 1)
        for i, stabilizer in enumerate(self.stabilizers_majorana):
            if not MajoranaOperator.commute(stabilizer, op):
                stabilizer.coff = -stabilizer.coff

    ##  USER：i*gamma*gamma_prime门，作用于majorana_index
    def N(self, majorana_index: int):
        op = MajoranaOperator([majorana_index], [majorana_index], 1j)
        for i, stabilizer in enumerate(self.stabilizers_majorana):
            if not MajoranaOperator.commute(stabilizer, op):
                stabilizer.coff = -stabilizer.coff

    ##  USER：fermionic phase gate，作用于majorana_index
    def P(self, majorana_index: int):
        pass

    ##  USER：S门，作用于pauli_index
    def S(self, pauli_index: int):
        for i, stabilizer in enumerate(self.stabilizers_pauli):
            assert isinstance(stabilizer, PauliOperator)
            left, middle, right = stabilizer.split(pauli_index)
            if len(middle.occupy_x) == 0:
                continue
            elif len(middle.occupy_x) == 1 and len(middle.occupy_z) == 0:
                middle.occupy_z = [pauli_index]
                middle.coff = 1j
            elif len(middle.occupy_x) == 1 and len(middle.occupy_z) == 1:
                middle.occupy_z = []
                middle.coff = 1j
            else:
                raise ValueError
            self.stabilizers_pauli[i] = left @ middle @ right

    ##  USER：CNOT门，作用于control_index,target_index，两者是qubits，前者是控制位
    def CX(self, control_index, target_index):
        for i, stabilizer in enumerate(self.stabilizers_pauli):
            assert isinstance(stabilizer, PauliOperator)
            left, middle, right = stabilizer.split(control_index)
            if len(middle.occupy_x) == 1:
                middle = PauliOperator([control_index, target_index], [], 1)
            if target_index > control_index:
                right_left, right_middle, right_right = right.split(target_index)
                if len(right_middle.occupy_z) == 1:
                    right_middle = PauliOperator([], [control_index, target_index], 1)
                self.stabilizers_pauli[i] = left @ middle @ right_left @ right_middle @ right_right
            elif target_index < control_index:
                left_left, left_middle, left_right = left.split(target_index)
                if len(left_middle.occupy_z) == 1:
                    left_middle = PauliOperator([], [control_index, target_index], 1)
                self.stabilizers_pauli[i] = left_left @ left_middle @ left_right @ middle @ right
            else:
                raise ValueError

    ##  USER：CN-NOT门，作用于control_index,target_index，前者是fermionic site控制位，后者是qubit目标位
    def CNX(self, control_index, target_index):
        for i in range(len(self.stabilizers_pauli)):
            stabilizer_pauli = self.stabilizers_pauli[i]
            stabilizer_majorana = self.stabilizers_majorana[i]
            left_control, middle_control, right_control = stabilizer_majorana.split(control_index)
            left_target, middle_target, right_target = stabilizer_pauli.split(target_index)
            majorana_product = PauliOperator([], [], 1)
            pauli_product = MajoranaOperator([], [], 1)
            if len(middle_control.occupy_x) == 1:
                majorana_product = PauliOperator([target_index], [], 1)
            if len(middle_control.occupy_z) == 1:
                majorana_product = majorana_product @ PauliOperator([target_index], [], 1)
            if len(middle_target.occupy_z) == 1:
                pauli_product = MajoranaOperator([target_index], [target_index], 1j)
            self.stabilizers_pauli[i] = majorana_product @ left_target @ middle_target @ right_target
            self.stabilizers_majorana[i] = left_control @ middle_control @ right_control @ pauli_product

    ##  USER：CU-NOT门，作用于control_index,target_index，前者是fermionic site控制位，后者是qubit目标位
    def CUX(self, control_index, target_index):
        for i in range(len(self.stabilizers_pauli)):
            stabilizer_pauli = self.stabilizers_pauli[i]
            stabilizer_majorana = self.stabilizers_majorana[i]
            left_control, middle_control, right_control = stabilizer_majorana.split(control_index)
            left_target, middle_target, right_target = stabilizer_pauli.split(target_index)
            majorana_product = PauliOperator([], [], 1)
            pauli_product = MajoranaOperator([], [], 1)
            if len(middle_control.occupy_z) == 1:
                majorana_product = PauliOperator([target_index], [], 1)
            if len(middle_target.occupy_z) == 1:
                pauli_product = MajoranaOperator([control_index], [], 1)
            self.stabilizers_pauli[i] = majorana_product @ left_target @ middle_target @ right_target
            self.stabilizers_majorana[i] = left_control @ middle_control @ right_control @ pauli_product

    ##  USER：CV-NOT门，作用于control_index,target_index，前者是fermionic site控制位，后者是qubit目标位
    def CVX(self, control_index, target_index):
        for i in range(len(self.stabilizers_pauli)):
            stabilizer_pauli = self.stabilizers_pauli[i]
            stabilizer_majorana = self.stabilizers_majorana[i]
            left_control, middle_control, right_control = stabilizer_majorana.split(control_index)
            left_target, middle_target, right_target = stabilizer_pauli.split(target_index)
            majorana_product = PauliOperator([], [], 1)
            pauli_product = MajoranaOperator([], [], 1)
            if len(middle_control.occupy_x) == 1:
                majorana_product = PauliOperator([target_index], [], 1)
            if len(middle_target.occupy_z) == 1:
                pauli_product = MajoranaOperator([], [control_index], 1)
            self.stabilizers_pauli[i] = majorana_product @ left_target @ middle_target @ right_target
            self.stabilizers_majorana[i] = left_control @ middle_control @ right_control @ pauli_product

    ##  USER：执行pauli_index上的X-error
    def x_error(self, pauli_index, p):
        if np.random.rand() < p:
            self.X(pauli_index)

    ##  USER：执行pauli_index上的Y-error
    def y_error(self, pauli_index, p):
        if np.random.rand() < p:
            self.Y(pauli_index)

    ##  USER：执行pauli_index上的Z-error
    def z_error(self, pauli_index, p):
        if np.random.rand() < p:
            self.Z(pauli_index)

    ##  USER：执行majorana_index上的U-error
    def u_error(self, majorana_index, p):
        if np.random.rand() < p:
            self.U(majorana_index)

    ##  USER：执行majorana_index上的V-error
    def v_error(self, majorana_index, p):
        if np.random.rand() < p:
            self.V(majorana_index)

    ##  USER：执行majorana_index上的N-error
    def n_error(self, majorana_index, p):
        if np.random.rand() < p:
            self.N(majorana_index)

    ##  USER：将系统在op上重置
    def clear(self, op):
        if isinstance(op, PauliOperator):
            stabilizers_now = self.stabilizers_pauli
        elif isinstance(op, MajoranaOperator):
            stabilizers_now = self.stabilizers_majorana
        else:
            raise NotImplementedError
        assert op.is_hermitian

        first_pauli = None
        first_index = -1
        for i in range(len(stabilizers_now)):
            if isinstance(op, PauliOperator):
                commute_flag = PauliOperator.commute(op, stabilizers_now[i])
            elif isinstance(op, MajoranaOperator):
                commute_flag = MajoranaOperator.commute(op, stabilizers_now[i])
            else:
                raise NotImplementedError

            if not commute_flag:
                if first_index == -1:
                    first_pauli = stabilizers_now[i]
                    first_index = i
                else:
                    stabilizers_now[i] = stabilizers_now[i] @ first_pauli
        if first_index == -1:
            matrix_pauli = PauliOperator.get_matrix(self.stabilizers_pauli, self.pauli_number)
            matrix_majorana = MajoranaOperator.get_matrix(self.stabilizers_majorana, self.majorana_number)
            if self.pauli_number == 0 and self.majorana_number != 0:
                matrix = matrix_majorana
            elif self.majorana_number == 0 and self.pauli_number != 0:
                matrix = matrix_pauli
            elif self.pauli_number != 0 and self.majorana_number != 0:
                matrix = ba.hstack(matrix_majorana, matrix_pauli)
            else:
                raise NotImplementedError

            if isinstance(op, MajoranaOperator):
                vector_majorana = op.get_vector(self.majorana_number)
                vector_pauli = ba.zeros(self.pauli_number * 2)
                vector = ba.hstack(vector_majorana, vector_pauli)
            else:
                vector_majorana = ba.zeros(self.majorana_number * 2)
                vector_pauli = op.get_vector(self.pauli_number)
                vector = ba.hstack(vector_majorana, vector_pauli)
            result = ba.solve(matrix, vector)
            op_mul_pauli = None
            op_mul_majorana = None
            flag = None
            for i in range(len(result)):
                if result[i] == 1:
                    if op_mul_pauli is None:
                        flag = i
                        op_mul_pauli = self.stabilizers_pauli[i]
                        op_mul_majorana = self.stabilizers_majorana[i]
                    else:
                        op_mul_pauli = op_mul_pauli @ self.stabilizers_pauli[i]
                        op_mul_majorana = op_mul_majorana @ self.stabilizers_majorana[i]
            if op_mul_pauli.coff * op_mul_majorana.coff == op.coff:
                pass
            else:
                assert flag is not None
                self.stabilizers_pauli[flag].coff = -self.stabilizers_pauli[flag].coff
        else:
            stabilizers_now[first_index] = op

    ##  USER：将系统在pauli_index上重置为0态
    def reset(self, pauli_index):
        op = PauliOperator([], [pauli_index], 1)
        self.clear(op)

    ##  USER：将系统在majorana_index上重置为空态
    def fermionic_reset(self, majorana_index):
        op = MajoranaOperator([majorana_index], [majorana_index], 1j)
        self.clear(op)


if __name__ == '__main__':
    from qutip import *
    gamma_0 = fcreate(2, 0) + fdestroy(2, 0)
    gamma_prime_0 = 1j * fdestroy(2, 0) - 1j * fcreate(2, 0)
    gamma_1 = fcreate(2, 1) + fdestroy(2, 1)
    gamma_prime_1 = 1j * fdestroy(2, 1) - 1j * fcreate(2, 1)

    T = (1j * (np.pi/2) * fcreate(2, 0)@fdestroy(2,1 ) +1j * (np.pi/2)*fcreate(2, 1) @fdestroy(2,0)).expm()

    all_list = [gamma_0 @ gamma_0,
                gamma_0, gamma_prime_0, gamma_1, gamma_prime_1,
                1j*gamma_0 @ gamma_prime_0, 1j*gamma_0 @ gamma_1, 1j*gamma_0 @gamma_prime_1, 1j*gamma_prime_0 @ gamma_1, 1j*gamma_prime_0@gamma_prime_1,
                1j*gamma_1 @ gamma_prime_1,
                1j*gamma_0 @ gamma_prime_0 @ gamma_1, 1j*gamma_0 @ gamma_prime_0 @ gamma_prime_1, 1*gamma_0 @ gamma_1 @ gamma_prime_1,
                1j*gamma_prime_0 @ gamma_1 @ gamma_prime_1,
                gamma_0 @ gamma_prime_0 @ gamma_1 @ gamma_prime_1,
                ]
    name = ['I',
            'gamma_0', 'gamma_prime_0', 'gamma_1', 'gamma_prime_1',
            'gamma_0@gamma_prime_0', 'gamma_0@gamma_1', 'gamma_0@gamma_prime_1','gamma_prime_0@gamma_1', 'gamma_prime_0@gamma_prime_1', 'gamma_1@gamma_prime_1',
            'gamma_0@gamma_prime_0@gamma_1','gamma_0@gamma_prime_0@gamma_prime_1', 'gamma_0@gamma_1@gamma_prime_1', 'gamma_prime_0@gamma_1@gamma_prime_1',
            'gamma_0@gamma_prime_0@gamma_1@gamma_prime_1']
    target = T @ gamma_prime_0 @ T.dag()
    for ii in range(len(all_list)):
        if (all_list[ii] - target).norm() < 1e-5:
            print(name[ii], 1)
        elif (-all_list[ii] - target).norm() < 1e-5:
            print(name[ii], -1)
