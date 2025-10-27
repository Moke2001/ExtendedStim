import copy
import qiskit
import numpy as np
import stim
import stimbposd
from qiskit.circuit import CircuitError
from qiskit.circuit.library import XGate
from Physics.MajoranaOperator import MajoranaOperator
from Physics.PauliOperator import PauliOperator
from Platform.Platform import Platform


class Circuit:
    #%%  USER：构造方法
    def __init__(self):
        self.majorana_number = 0
        self.pauli_number = 0
        self.sequence = []
        self.noise = []
        self.measurements = []
        self.detectors = []
        self.observables = []
        self._dem=None

    #%%  USER：重载运算符
    ##  USER：获取序列中的元素
    def __getitem__(self, item):
        return self.sequence[item]

    ##  USER：设置序列中的元素
    def __setitem__(self, key, value):
        self.sequence[key] = value

    #%%  USER：对象方法
    ##  USER：添加量子线路组分操作
    def append(self, name, target, *args):

        ##  添加single-qubit gate
        if name == 'X' or name == 'Y' or name == 'Z' or name == 'H' or name == 'S':
            assert len(args) == 0
            self.sequence.append({'name': name, 'target': target})

        ##  添加single-fermionic-site gate
        elif name == 'U' or name == 'V' or name == 'N' or name == 'P':
            assert len(args) == 0
            self.sequence.append({'name': name, 'target': target})

        ##  添加受控非门
        elif name == 'CX' or name == 'CUX' or name == 'CVX' or name == 'CNX':
            assert len(args) == 0
            assert len(target) == 2
            self.sequence.append({'name': name, 'target': target})

        ##  添加single-qubit and single-fermionic-site上的测量
        elif name == 'MZ' or name == 'MN':
            assert len(args) == 0 or len(args) == 1
            if name == 'MZ':
                self.sequence.append({'name': 'MPP', 'target': PauliOperator([], [target], 1)})
            else:
                self.sequence.append({'name': 'MPP', 'target': MajoranaOperator([target], [target], 1j)})
            self.measurements.append(len(self.sequence) - 1)
            if len(args) == 1:
                self.sequence.append({'name': 'M_ERROR', 'p': args[0]})
                self.noise.append(len(self.sequence) - 1)

        ##  添加qubit上的去极化噪声
        elif name == 'DEPOLARIZE1':
            assert len(args) == 1
            fix=27/26
            self.sequence.append({'name': 'X_ERROR', 'target': target, 'p': args[0]*fix / 3})
            self.noise.append(len(self.sequence) - 1)
            self.sequence.append({'name': 'Y_ERROR', 'target': target, 'p': args[0]*fix / 3})
            self.noise.append(len(self.sequence) - 1)
            self.sequence.append({'name': 'Z_ERROR', 'target': target, 'p': args[0]*fix / 3})
            self.noise.append(len(self.sequence) - 1)

        ##  添加fermionic-site上的去极化噪声
        elif name == 'FDEPOLARIZE1':
            assert len(args) == 1
            fix = 27 / 26
            self.sequence.append({'name': 'U_ERROR', 'target': target, 'p': args[0]*fix / 3})
            self.noise.append(len(self.sequence) - 1)
            self.sequence.append({'name': 'V_ERROR', 'target': target, 'p': args[0]*fix / 3})
            self.noise.append(len(self.sequence) - 1)
            self.sequence.append({'name': 'N_ERROR', 'target': target, 'p': args[0]*fix / 3})
            self.noise.append(len(self.sequence) - 1)

        ##  添加string算符的测量
        elif name == 'MPP':
            self.sequence.append({'name': name, 'target': target})
            self.measurements.append(len(self.sequence) - 1)
            if len(args) == 1:
                self.sequence.append({'name': 'M_ERROR', 'p': args[0]})
                self.noise.append(len(self.sequence) - 1)

        ##  添加监视器
        elif name == 'DETECTOR':
            assert all(target[i] < 0 for i in range(len(target)))
            together = [len(self.measurements) + temp for temp in target]
            self.detectors.append(together)

        ##  添加可观测量
        elif name == 'OBSERVABLE_INCLUDE':
            assert all(target[i] < 0 for i in range(len(target)))
            together = [len(self.measurements) + temp for temp in target]
            self.observables.append(together)

        ##  添加qubit重置
        elif name == 'R':
            if isinstance(target, int):
                self.sequence.append({'name': name, 'target': target})
                if target == self.pauli_number:
                    self.pauli_number += 1
                elif target > self.pauli_number or target < 0:
                    raise ValueError("R gate target must be consecutive")
                else:
                    pass

        ##  添加fermionic-site重置
        elif name == 'FR':
            if isinstance(target, int):
                self.sequence.append({'name': name, 'target': target})
                if target == self.majorana_number :
                    self.majorana_number = target+1
                elif target > self.majorana_number or target < 0:
                    raise ValueError("FR gate target must be consecutive")
                else:
                    pass
        else:
            raise NotImplementedError

    ##  USER：生成无噪声的线路
    def ideal_circuit(self):
        sequence = copy.deepcopy(self.sequence)
        for i in range(len(self.noise)):
            gate = sequence[self.noise[i]]
            assert isinstance(gate, dict)
            gate['p'] = 0
        ideal_circuit = Circuit()
        ideal_circuit.majorana_number = self.majorana_number
        ideal_circuit.pauli_number = self.pauli_number
        ideal_circuit.sequence = sequence
        ideal_circuit.measurements = self.measurements
        ideal_circuit.detectors = self.detectors
        ideal_circuit.observables = self.observables
        ideal_circuit.noise = self.noise
        return ideal_circuit

    ##  USER：执行线路并返回测量结果
    def execute(self):
        platform = Platform()
        platform.initialize(self.majorana_number, self.pauli_number)
        measurement_sample = np.empty(len(self.measurements), dtype=int)
        flag_measurement = 0
        for i, gate in enumerate(self.sequence):
            p = None
            target = None
            name = gate['name']
            if 'target' in gate:
                target = gate['target']
            if 'p' in gate:
                p = gate['p']
            if name == 'X':
                assert target is not None
                platform.X(target)
            elif name == 'Y':
                assert target is not None
                platform.Y(target)
            elif name == 'Z':
                assert target is not None
                platform.Z(target)
            elif name == 'H':
                assert target is not None
                platform.H(target)
            elif name == 'S':
                assert target is not None
                platform.S(target)
            elif name == 'U':
                assert target is not None
                platform.U(target)
            elif name == 'V':
                assert target is not None
                platform.V(target)
            elif name == 'N':
                assert target is not None
                platform.N(target)
            elif name == 'R':
                assert target is not None
                platform.reset(target)
            elif name == 'FR':
                assert target is not None
                platform.fermionic_reset(target)
            elif name == 'X_ERROR':
                assert target is not None
                assert p is not None
                platform.x_error(target, p)
            elif name == 'Y_ERROR':
                assert target is not None
                assert p is not None
                platform.y_error(target, p)
            elif name == 'Z_ERROR':
                assert target is not None
                assert p is not None
                platform.z_error(target, p)
            elif name == 'U_ERROR':
                assert target is not None
                assert p is not None
                platform.u_error(target, p)
            elif name == 'V_ERROR':
                assert target is not None
                assert p is not None
                platform.v_error(target, p)
            elif name == 'N_ERROR':
                assert target is not None
                assert p is not None
                platform.n_error(target, p)
            elif name == 'CX':
                assert target is not None
                platform.CX(target[0], target[1])
            elif name == 'CUX':
                assert target is not None
                platform.CUX(target[0], target[1])
            elif name == 'CVX':
                assert target is not None
                platform.CVX(target[0], target[1])
            elif name == 'CNX':
                assert target is not None
                platform.CNX(target[0], target[1])
            elif name == 'MPP':
                assert target is not None
                measurement_sample[flag_measurement] = platform.measure(target)
                flag_measurement += 1
            elif name == 'M_ERROR':
                assert p is not None
                if np.random.rand() < p:
                    measurement_sample[flag_measurement-1] = -measurement_sample[flag_measurement-1]
            else:
                raise NotImplementedError

        detector_sample = np.empty(len(self.detectors), dtype=bool)
        flag_detector = 0
        for i, detector in enumerate(self.detectors):
            value = measurement_sample[detector][0]
            detector_sample[flag_detector] = False
            for temp in measurement_sample[detector]:
                if value==temp:
                    continue
                else:
                    detector_sample[flag_detector]=True
                    break
            flag_detector += 1

        observable_sample = np.empty(len(self.observables), dtype=bool)
        flag_observable = 0
        for i, observable in enumerate(self.observables):
            if len(observable)==1:
                if measurement_sample[observable][0]==1:
                    observable_sample[flag_observable] = False
                else:
                    observable_sample[flag_observable] = True
            else:
                value = measurement_sample[observable][0]
                observable_sample[flag_observable] = False
                for temp in measurement_sample[observable]:
                    if value==temp:
                        continue
                    else:
                        observable_sample[flag_observable]=True
                        break
            flag_observable += 1
        return measurement_sample, detector_sample, observable_sample,platform

    ##  USER：生成检测错误模型
    def detector_error_model(self) -> stim.DetectorErrorModel:
        if self._dem is not None:
            return self._dem
        ideal_circuit = self.ideal_circuit()
        measurement_sample_origin, detector_sample_origin, observable_sample_origin,platform = ideal_circuit.execute()

        ##  执行检验线路的稳定性
        for time in range(5):
            measurement_sample, detector_sample, observable_sample,platform = ideal_circuit.execute()
            assert np.all(detector_sample==detector_sample_origin),f'原始线路的detector不是稳定的'
            assert np.all(observable_sample==observable_sample_origin),f'原始线路的observable不是稳定的'

        errors = []
        dem_str = ''
        for i in range(len(ideal_circuit.noise)):
            order = ideal_circuit.noise[i]
            gate_ideal = ideal_circuit.sequence[order]
            assert isinstance(gate_ideal, dict)
            gate_ideal['p'] = 1.1
            gate = self.sequence[order]
            assert isinstance(gate, dict)
            p = gate['p']
            measurement_sample, detector_sample, observable_sample,platform = ideal_circuit.execute()
            detector_sample_diff = [detector_sample_origin[j] ^ detector_sample[j] for j in range(len(detector_sample))]
            observable_sample_diff = [observable_sample_origin[j] ^ observable_sample[j] for j in range(len(observable_sample))]
            errors.append(len(errors))
            detectors_trigger = np.where(np.array(detector_sample_diff) == True)[0]
            observables_trigger = np.where(np.array(observable_sample_diff) == True)[0]
            if len(detectors_trigger) > 0 or len(observables_trigger) > 0:
                temp = f'error({p})'
                for index in detectors_trigger:
                    temp = temp + f' D{index}'
                for index in observables_trigger:
                    temp = temp + f' L{index}'
                dem_str += ('\n' + temp)
            gate_ideal['p'] = 0
        dem = stim.DetectorErrorModel(dem_str)
        self._dem=dem
        return dem

    ##  USER：生成编译后的采样器
    def compiler_sampler(self):
        dem = self.detector_error_model()
        return dem.compile_sampler()

    ##  USER：生成解码函数
    def decoder(self, method):
        dem = self.detector_error_model()
        if method == 'bposd':
            dec=stimbposd.bp_osd.BPOSD(model=dem,bp_method='minimum_sum')
        else:
            raise NotImplementedError
        return dec

    ##  USER：执行线路并返回错误率
    def sample(self,sample_number):
        sampler=self.compiler_sampler()
        decoder=self.decoder('bposd')
        detector_data, obs_data, error_data = sampler.sample(shots=sample_number)
        predictions = decoder.decode_batch(detector_data)
        num_errors = 0
        for shot in range(sample_number):
            actual_for_shot = obs_data[shot]
            predicted_for_shot = predictions[shot]
            if not np.array_equal(actual_for_shot, predicted_for_shot):
                num_errors += 1
        return num_errors/sample_number

    ##  USER：绘制线路图
    def draw(self, filename):
        # 绘制一个带有barriers和更多寄存器中，绘制一个新的电路
        F = qiskit.QuantumRegister(self.majorana_number, name='F')
        Q = qiskit.QuantumRegister(self.pauli_number, name='Q')
        C = qiskit.ClassicalRegister(1, name='C')
        A = qiskit.QuantumRegister(1, name='A')
        circuit_qiskit = qiskit.QuantumCircuit(F, Q, C,A)
        cux = qiskit.circuit.ControlledGate(name='CUX', num_qubits=2, params=[], label=None, num_ctrl_qubits=1, base_gate=XGate())
        cvx = qiskit.circuit.ControlledGate(name='CVX', num_qubits=2, params=[], label=None, num_ctrl_qubits=1, base_gate=XGate())
        cnx = qiskit.circuit.ControlledGate(name='CNX', num_qubits=2, params=[], label=None, num_ctrl_qubits=1, base_gate=XGate())
        x_error = qiskit.circuit.Gate('X', 1, label=None, params=[])
        y_error = qiskit.circuit.Gate('Y', 1, label=None, params=[])
        z_error = qiskit.circuit.Gate('Z', 1, label=None, params=[])
        u_error = qiskit.circuit.Gate('U', 1, label=None, params=[])
        v_error = qiskit.circuit.Gate('V', 1, label=None, params=[])
        n_error = qiskit.circuit.Gate('N', 1, label=None, params=[])
        m_error = qiskit.circuit.Gate('M', 1, label=None, params=[])
        rreset = qiskit.circuit.Gate('R', 1, label=None, params=[])
        for gate in self.sequence:
            if gate['name'] == 'R':
                circuit_qiskit.append(rreset, [Q[gate['target']]])
            elif gate['name'] == 'FR':
                circuit_qiskit.append(rreset, [F[gate['target']]])
            elif gate['name'] == 'X':
                circuit_qiskit.x(Q[gate['target']])
            elif gate['name'] == 'Y':
                circuit_qiskit.y(Q[gate['target']])
            elif gate['name'] == 'Z':
                circuit_qiskit.z(Q[gate['target']])
            elif gate['name'] == 'H':
                circuit_qiskit.h(Q[gate['target']])
            elif gate['name'] == 'S':
                circuit_qiskit.s(Q[gate['target']])
            elif gate['name'] == 'CX':
                circuit_qiskit.cx(Q[gate['target'][0]], Q[gate['target'][1]])
            elif gate['name'] == 'CUX':
                circuit_qiskit.append(cux, [F[gate['target'][0]],Q[gate['target'][1]]])
            elif gate['name'] == 'CVX':
                circuit_qiskit.append(cvx, [F[gate['target'][0]],Q[gate['target'][1]]])
            elif gate['name'] == 'CNX':
                circuit_qiskit.append(cnx, [F[gate['target'][0]],Q[gate['target'][1]]])
            elif gate['name'] == 'MPP':
                op = gate['target']
                if isinstance(op, MajoranaOperator):
                    f_flag_x = op.occupy_x
                    f_flag_z = op.occupy_z
                    f_flag_n = np.intersect1d(f_flag_x, f_flag_z)
                    f_flag_x = np.setdiff1d(f_flag_x, f_flag_n)
                    f_flag_z = np.setdiff1d(f_flag_z, f_flag_n)
                    f = np.concatenate([f_flag_x, f_flag_z, f_flag_n])
                    if len(f)>1:
                        mppx=qiskit.circuit.ControlledGate(name='MPPX', num_qubits=len(f)+1, params=[], label=None, num_ctrl_qubits=len(f), base_gate=XGate())
                        circuit_qiskit.append(mppx, F[f.tolist()]+[A[0]])
                        circuit_qiskit.measure(A[0], C[0])
                    else:
                        circuit_qiskit.measure(F[f[0]], C[0])
                elif isinstance(op, PauliOperator):
                    p_flag_x = op.occupy_x
                    p_flag_z = op.occupy_z
                    p_flag_y = np.intersect1d(p_flag_x, p_flag_z)
                    p_flag_x = np.setdiff1d(p_flag_x, p_flag_y)
                    p_flag_z = np.setdiff1d(p_flag_z, p_flag_y)
                    p = np.concatenate([p_flag_x, p_flag_y, p_flag_z])
                    if len(p)>1:
                        mppx=qiskit.circuit.ControlledGate(name='MPPX', num_qubits=len(p)+1, params=[], label=None, num_ctrl_qubits=len(p), base_gate=XGate())
                        circuit_qiskit.append(mppx, Q[p.tolist()]+[A[0]])
                        circuit_qiskit.measure(A[0], C[0])
                    else:
                        circuit_qiskit.measure(Q[p[0]], C[0])
                else:
                    raise CircuitError("cannot set parameters on immutable base gate")
            elif gate['name'] == 'M_ERROR':
                circuit_qiskit.append(m_error, [Q[gate['target']]])
            elif gate['name'] == 'X_ERROR':
                circuit_qiskit.append(x_error, [Q[gate['target']]])
            elif gate['name'] == 'Y_ERROR':
                circuit_qiskit.append(y_error, [Q[gate['target']]])
            elif gate['name'] == 'Z_ERROR':
                circuit_qiskit.append(z_error, [Q[gate['target']]])
            elif gate['name'] == 'U_ERROR':
                circuit_qiskit.append(u_error, [Q[gate['target']]])
            elif gate['name'] == 'V_ERROR':
                circuit_qiskit.append(v_error, [Q[gate['target']]])
            elif gate['name'] == 'N_ERROR':
                circuit_qiskit.append(n_error, [Q[gate['target']]])
        red='#E77081'
        blue='#5375CD'
        green='#00857B'
        grey='#8C92AC'
        purple='#5D548C'
        orange='#F15D22'
        pink='#FFACC5'
        cyan='#C9DCC4'
        circuit_qiskit.draw(output='mpl', filename=filename, style={
            'displaycolor': {'X': red, 'Y': red, 'Z': red,
                             'U': red, 'V': red, 'N': red,
                             'M': red, 'R': cyan,'measure':grey,
                             'x':blue, 'y':blue, 'z':blue,'s':blue,
                             'CNX':blue,'CUX':purple,'CVX':pink,
                             'cx':None,'MPPX':grey
                             },
            'fontsize': 15
        })