"""""
Module: Circuit
"""""
import copy
import os
import time
from multiprocessing import Pool, cpu_count
import pymatching
import qiskit
import numpy as np
import stim
import stimbposd
from qiskit.circuit import CircuitError
from qiskit.circuit.library import XGate, ZGate
from extendedstim.Physics.MajoranaOperator import MajoranaOperator
from extendedstim.Physics.PauliOperator import PauliOperator
from extendedstim.Platform.Frame import Frame
from extendedstim.Platform.Platform import Platform
from extendedstim.tools.TypingTools import isinteger, islist


class Circuit:
    __slots__=['majorana_number', 'pauli_number', 'sequence', 'ideal_sequence',
               '_sequence', '_noises', '_measurements', '_detectors', '_observables',
               '_prototype', '_reference_measurements', '_alert']

    ##  CHAPTER：====Constructor====
    def __init__(self):
        """""
        circuit.append({'name': 'TRAP','majorana_number':10,'pauli_number':10})来指定初始的qubits and fermionic sites数目
        circuit.append({'name': 'X','target':0})等语句来添加gates，具体参考append方法
        """""

        ##  Public attributes
        self.majorana_number=0  # fermionic sites的数目
        self.pauli_number=0  # qubits的数目
        self.sequence=[]  # 量子线路的操作序列
        self.ideal_sequence=[]  # 理想量子线路的操作序列

        ##  Private attributes
        self._sequence:list[dict]=[]  # 记录量子线路的操作序列在真实计算使用
        self._noises=[]  # 记录noise在_sequence中的索引
        self._measurements=[]  # 记录measurement在_sequence中的索引
        self._detectors=[]  # 记录detector在_measurements中的索引
        self._observables=[]  # 记录observable在_measurements中的索引
        self._prototype=None  # 量子线路的detector error model的原型，probability可以修改
        self._reference_measurements=None  # 量子线路的参考测量结果，用于后续sample
        self._alert=False  # 记录是否有alert，用于禁用append方法

    ##  SECTION：----Get gate----
    def __getitem__(self, item: int) -> dict:
        return self.sequence[item]

    ##  SECTION：-----Get reference sample-----
    @property
    def _reference_sample(self)->np.ndarray:
        if self._reference_measurements is None:

            ##  PART：----Initialization----
            platform=Platform()  # 生成量子平台
            measurement_sample=np.empty(len(self._measurements), dtype=int)  # 生成测量值的样本数组
            flag_measurement=0

            ##  PART：----Traverse the sequence----
            for i, gate in enumerate(self._sequence):
                name=gate['name']
                if name in ['X', 'Y', 'Z', 'H', 'S', 'U', 'V', 'N']:  # 执行单门
                    eval('platform.'+name.lower()+'(gate[\'target\'])')
                elif name in ['CX', 'CNX', 'BRAID', 'CNN']:  # 执行双门
                    eval('platform.'+name.lower()+'(gate[\'target\'][0],gate[\'target\'][1])')
                elif name=='R':  # 执行重置
                    platform.reset(gate['target'])
                elif name=='MPP':  # 执行测量
                    measurement_sample[flag_measurement]=platform.measure(gate['target'])
                    flag_measurement+=1
                elif name=='TRAP':  # 初始化一定数目的qubits和fermionic sites
                    platform.trap(self.majorana_number,self.pauli_number)
            self._reference_measurements=measurement_sample

        ##  PART：----Return the observable results----
        return self._reference_measurements

    ##  SECTION：----Get the prototype of detector error model----
    @property
    def prototype(self)->list[list]:
        if self._prototype is None:

            ##  PART：----Preprocessing----
            ##  Get noiseless results
            measurement_sample_origin=self._reference_sample
            detector_sample_origin=diff(measurement_sample_origin, self._detectors)
            observable_sample_origin=diff(measurement_sample_origin, self._observables)

            ##  Check the stability of the original circuit
            for _ in range(10):
                measurement_sample, detector_sample, observable_sample=self.sample(noiseless=True)
                assert np.all(detector_sample==detector_sample_origin), f'原始线路的detector不是稳定的'
                assert np.all(observable_sample==observable_sample_origin), f'原始线路的observable不是稳定的'

            ##  PART：----Sample the circuit with each noise, analyze the effect of each noise on detectors and observables----
            ##  Divide the noises into groups
            group_number=cpu_count()//2  # 分组数量
            sample_number=len(self._noises)//group_number  # 每个分组的噪声数量
            noise_group=[self._noises[temp*sample_number:(temp+1)*sample_number]
                         for temp in range(group_number)]  # 噪声组合的列表
            if len(self._noises[group_number*sample_number:])>0:  # 处理最后一组噪声
                noise_group.append(self._noises[group_number*sample_number:])

            ##  Sample the circuit with each noise group, analyze the effect of each noise group on detectors and observables
            with Pool(processes=len(noise_group)) as pool:
                results=[pool.apply_async(self._sample_batch, args=(noise_group[i],)) for i in range(len(noise_group))]
                final_results=[result.get() for result in results]

            ##  PART：----Analyze the effect of each noise group on detectors and observables----
            flag=0  # 记录当前noise的序号
            prototype=[]  # 初始化原型
            for i in range(len(final_results)):
                for j in range(len(final_results[i])):

                    measurement_sample, detector_sample, observable_sample=final_results[i][j]

                    ##  找到触发的detector和observable
                    detector_sample_diff=[detector_sample_origin[j]^detector_sample[j] for j in range(len(detector_sample))]
                    observable_sample_diff=[observable_sample_origin[j]^observable_sample[j] for j in range(len(observable_sample))]

                    ##  记录触发的detector和observable的位置
                    detectors_trigger=np.where(np.array(detector_sample_diff)==True)[0]
                    observables_trigger=np.where(np.array(observable_sample_diff)==True)[0]

                    ##  合成错误语句
                    if len(detectors_trigger)>0 or len(observables_trigger)>0:
                        temp_error=f'error({self._sequence[self._noises[flag]]['p']}) '
                        temp_trigger=''
                        for index in detectors_trigger:
                            temp_trigger=temp_trigger+f' D{index}'
                        for index in observables_trigger:
                            temp_trigger=temp_trigger+f' L{index}'
                        if self._sequence[self._noises[flag]]['name']=='M_ERROR':
                            prototype.append([temp_error, temp_trigger, 'M_ERROR'])
                        elif 'type' not in self._sequence[self._noises[flag]]:
                            prototype.append([temp_error, temp_trigger, 'ONE_ERROR'])
                        elif self._sequence[self._noises[flag]]['type']=='DEP1':
                            prototype.append([temp_error, temp_trigger, 'DEP1'])
                        elif self._sequence[self._noises[flag]]['type']=='DEP2':
                            prototype.append([temp_error, temp_trigger, 'DEP2'])
                        else:
                            raise NotImplementedError(f'未实现的错误类型')
                        flag+=1

        ##  PART：----Return the prototype of detector error model----
        return prototype

    ##  SECTION：----Add circuit component operations----
    def append(self, params):
        """""
        {
        'name': str，线路操作的名称
        'target': list or int，操作作用的对象
        'p': float，操作对应的概率，不一定有
        'pauli_number': int，强制初始化的qubits数目
        'majorana_number': int，强制初始化的fermionic sites数目
        }
        支持的线路操作名称：
        'X', 'Y', 'Z', 'H', 'S'：{'name':str,'target':int or list}, single qubit上的qugate
        'X_ERROR', 'Y_ERROR', 'Z_ERROR', 'DEPOLARIZE1': {'name':str,'target':int or list}, single qubit上的噪声
        'XX_ERROR', 'XY_ERROR', ... , 'DEPOLARIZE2' {'name':str,'target':list}, two qubits上的噪声
        'U', 'V', 'N', 'P'：{'name':str,'target':int or list}, single fermionic site上的fgate
        'U_ERROR', 'V_ERROR', 'N_ERROR', 'FDEPOLARIZE1': {'name':str,'target':int or list}, single fermionic site上的噪声
        'UU_ERROR', 'UV_ERROR', ... , 'FDEPOLARIZE2': {'name':str,'target':int or list}, two fermionic sites上的噪声
        'MIXTURE'：{'name':str,'target':list,'p':float}, single qubit and single fermionic site上的mixture depolarization
        'CX', 'CNX', 'CNN', 'BRAID': {'name':str,'target':list}, two qubit or fermionic sites上的gates
        'R': {'name':str,'target':int or list}, single qubit or single fermionic sites上的重置到空态或0态
        'MZ', 'MN': {'name':str,'target':int or list}, single qubit or single fermionic sites上的measurement
        'MPP': {'name':str,'target':list or Operator}, Pauli string operators or Majorana string operators的measurement
        'TRAP'：{'name':str,'pauli_number':int,'majorana_number':int}, 强制初始化到空态和0态
        'DETECTOR': {'name':str,'target':list of int}，探测器
        'OBSERVABLE_INCLUDE': {'name':str,'target':list of int}，可观测量
        """""

        ##  PART：----Data preprocessing----
        assert isinstance(params, dict)
        assert not self._alert
        assert 'name' in params
        name=params["name"]  # 线路操作的名称

        ##  PART：----Add single-body gate----
        if name in ['X', 'Y', 'Z', 'H', 'S', 'X_ERROR', 'Y_ERROR', 'Z_ERROR', 'U', 'V', 'N', 'P', 'U_ERROR', 'V_ERROR', 'N_ERROR']:

            ##  作用在整数index对应目标上
            if isinteger(params['target']):
                self._sequence.append(params.copy())
                self.sequence.append(self._sequence[-1])
                self.ideal_sequence.append(self._sequence[-1])

            ##  作用在一系列目标上
            elif islist(params['target']):
                for temp in params['target']:
                    params_temp={'name': name, 'target': temp}
                    self.append(params_temp)

        ##  PART：----Add two-body gate----
        elif name in ['CX', 'CNX', 'CNN',
                      'XX_ERROR','XY_ERROR','XZ_ERROR','YX_ERROR','YY_ERROR','YZ_ERROR','ZX_ERROR','ZY_ERROR','ZZ_ERROR',
                      'UU_ERROR','UV_ERROR','UN_ERROR','VU_ERROR','VV_ERROR','VN_ERROR','NU_ERROR','NV_ERROR','NN_ERROR',
                      'UX_ERROR','UY_ERROR','UZ_ERROR','VX_ERROR','VY_ERROR','VZ_ERROR','NX_ERROR','NY_ERROR','NZ_ERROR',]:

            ##  添加单个two gate
            if islist(params['target']) and len(params['target'])==2 and isinteger(params['target'][0]) and isinteger(params['target'][1]):
                self._sequence.append(params.copy())
                self.sequence.append(self._sequence[-1])
                self.ideal_sequence.append(self._sequence[-1])

            ##  用单个列表添加多个two gates
            elif islist(params['target']) and isinteger(params['target'][0]):
                for i in range(len(params['target'])//2):
                    params_temp={'name': name, 'target': [params['target'][2*i], params['target'][2*i+1]]}
                    self.append(params_temp)

            ##  用多个列表添加多个two gates
            elif islist(params['target']):
                for temp in params['target']:
                    self.append({'name': name, 'target': temp})

            ##  其他情况抛出异常
            else:
                raise ValueError("Gate must be applied to two")

        ##  PART：----Add braid gate----
        elif name=='BRAID' or name=='braid':

            if islist(params['target']) and len(params['target'])==2 and isinteger(params['target'][0]) and isinteger(params['target'][1]):
                x=params['target'][0]  # 控制位
                y=params['target'][1]  # 目标位

                ##  控制位等于目标位时，根据verse参数判断是否添加N门
                if x==y and 'verse' not in params:
                    self._sequence.append({'name': 'P', 'target': x})
                    self.sequence.append(self._sequence[-1])
                    self.ideal_sequence.append(self._sequence[-1])
                elif x==y and 'verse' in params and params['verse']==True:
                    self._sequence.append({'name': 'P', 'target': y})
                    self.sequence.append(self._sequence[-1])
                    self.ideal_sequence.append(self._sequence[-1])
                    self._sequence.append({'name': 'N', 'target': y})
                    self.sequence.append(self._sequence[-1])
                    self.ideal_sequence.append(self._sequence[-1])

                ##  一般情况下braid gate作用在两个不同的fermionic sites上
                elif x!=y:
                    self._sequence.append({'name': 'BRAID', 'target': [x, y]})
                    self.sequence.append(self._sequence[-1])
                    self.ideal_sequence.append(self._sequence[-1])

            ##  用单个列表添加多个two gates
            elif islist(params['target']) and isinteger(params['target'][0]):
                for i in range(len(params['target'])//2):
                    params_temp={'name': name, 'target': [params['target'][2*i], params['target'][2*i+1]]}
                    self.append(params_temp)

            ##  用多个列表添加多个two gates
            elif islist(params['target']):
                for temp in params['target']:
                    self.append({'name': name, 'target': temp})

            ##  其他情况抛出异常
            else:
                raise ValueError("CX, CNX, CNN gate must be applied to two")

        ##  PART：----Add single-qubit depolarization----
        elif name=='DEPOLARIZE1':
            assert 'p' in params

            ##  作用在整数index对应目标上
            if isinteger(params['target']):
                self.sequence.append({'name': 'DEPOLARIZE1', 'target': params['target'], 'p': params['p']})
                fix=(1-np.sqrt(1-4*params["p"]/3))/2
                self._sequence.append({'name': 'X_ERROR', 'target': params["target"], 'p': fix,'type':'DEP1'})
                self._noises.append(len(self._sequence)-1)
                self._sequence.append({'name': 'Y_ERROR', 'target': params["target"], 'p': fix,'type':'DEP1'})
                self._noises.append(len(self._sequence)-1)
                self._sequence.append({'name': 'Z_ERROR', 'target': params["target"], 'p': fix,'type':'DEP1'})
                self._noises.append(len(self._sequence)-1)

            ##  作用在一系列目标上
            elif islist(params['target']):
                for temp in params['target']:
                    self.append({'name': params['name'], 'target': temp, 'p': params["p"]})

        ##  PART：----Add two-qubit depolarization----
        elif name=='DEPOLARIZE2':
            assert 'p' in params

            ##  用单个列表添加单个
            if islist(params['target']) and isinteger(params['target'][0]) and isinteger(params['target'][1]) and len (params['target'])==2:
                self.sequence.append({'name': 'DEPOLARIZE2', 'target': params['target'], 'p': params['p'],'type':'DEP2'})
                fix =(1/2)*(1 - (1 - 16 * params["p"] / 15)**(1/8))
                for case in ['X_ERROR', 'Y_ERROR', 'Z_ERROR']:
                    self._sequence.append({'name': case, 'target': params["target"][0], 'p': fix,'type':'DEP2'})
                    self._noises.append(len(self._sequence) - 1)
                for case in ['X_ERROR', 'Y_ERROR', 'Z_ERROR']:
                    self._sequence.append({'name': case, 'target': params["target"][1], 'p': fix,'type':'DEP2'})
                    self._noises.append(len(self._sequence) - 1)
                for case in ['XX_ERROR', 'XY_ERROR', 'XZ_ERROR','YX_ERROR','YY_ERROR','YZ_ERROR','ZX_ERROR','ZY_ERROR','ZZ_ERROR']:
                    self._sequence.append({'name': case, 'target': params["target"], 'p': fix,'type':'DEP2'})
                    self._noises.append(len(self._sequence) - 1)

            ##  用单个列表添加多个
            elif islist(params['target']) and isinteger(params['target'][0]) and len(params['target'])>2:
                for i in range(len(params['target'])//2):
                    self.append({'name': name, 'target':[params['target'][2*i], params['target'][2*i+1]], 'p': params['p']})

            ##  用多个列表添加多个
            elif islist(params['target']) and islist(params['target'][0]) and isinteger(params['target'][0][0]):
                for temp in params['target']:
                    self.append({'name': params['name'], 'target': temp, 'p': params["p"]})

            ##  其他情况抛出异常
            else:
                raise ValueError("Gate must be applied to two")

        ##  PART：----Add single-site depolarization----
        elif name=='FDEPOLARIZE1':
            assert 'p' in params

            ##  作用在整数index对应目标上
            if isinteger(params['target']):
                self.sequence.append({'name': 'FDEPOLARIZE1', 'target': params['target'], 'p': params['p']})
                fix=(1-np.sqrt(1-4*params["p"]/3))/2
                self._sequence.append({'name': 'U_ERROR', 'target': params["target"], 'p': fix,'type':'DEP1'})
                self._noises.append(len(self._sequence)-1)
                self._sequence.append({'name': 'V_ERROR', 'target': params["target"], 'p': fix,'type':'DEP1'})
                self._noises.append(len(self._sequence)-1)
                self._sequence.append({'name': 'N_ERROR', 'target': params["target"], 'p': fix,'type':'DEP1'})
                self._noises.append(len(self._sequence)-1)

            ##  作用在一系列目标上
            elif islist(params['target']):
                for temp in params['target']:
                    self.append({'name': params['name'], 'target': temp, 'p': params["p"]})

        ##  PART：----Add two-site depolarization----
        elif name == 'FDEPOLARIZE2':
            assert 'p' in params

            ##  用一个列表添加一个
            if islist(params['target']) and isinteger(params['target'][0]):
                self.sequence.append({'name': 'FDEPOLARIZE2', 'target': params['target'], 'p': params['p']})
                fix = (1 / 2) * (1 - (1 - 16 * params["p"] / 15) ** (1 / 8))
                for case in ['U_ERROR', 'V_ERROR', 'N_ERROR']:
                    self._sequence.append({'name': case, 'target': params["target"][0], 'p': fix,'type':'DEP2'})
                    self._noises.append(len(self._sequence) - 1)
                for case in ['U_ERROR', 'V_ERROR', 'N_ERROR']:
                    self._sequence.append({'name': case, 'target': params["target"][1], 'p': fix,'type':'DEP2'})
                    self._noises.append(len(self._sequence) - 1)
                for case in ['UU_ERROR', 'UV_ERROR', 'UN_ERROR','VU_ERROR','VV_ERROR','VN_ERROR','NU_ERROR','NV_ERROR','NN_ERROR']:
                    self._sequence.append({'name': case, 'target': params["target"], 'p': fix,'type':'DEP2'})
                    self._noises.append(len(self._sequence) - 1)

            ##  用一个列表添加多个
            elif islist(params['target']) and isinteger(params['target'][0]) and len(params['target'])>2:
                for i in range(len(params['target'])//2):
                    self.append({'name': name, 'target':[params['target'][2*i], params['target'][2*i+1]], 'p': params['p']})

            ##  用多个列表添加多个
            elif islist(params['target']) and islist(params['target'][0]) and isinteger(params['target'][0][0]):
                for temp in params['target']:
                    self.append({'name': params['name'], 'target': temp, 'p': params["p"]})

            ##  其他情况抛出异常
            else:
                raise ValueError("Gate must be applied to two")

        ##  PART：----Add site-qubit depolarization----
        elif name=='MIXTURE':
            assert 'p' in params

            ##  用单个列表添加单个
            if islist(params['target']) and isinteger(params['target'][0]):
                self.sequence.append({'name': 'MIXTURE2', 'target': params['target'], 'p': params['p']})
                fix = (1 / 2) * (1 - (1 - 16 * params["p"] / 15) ** (1 / 8))
                for case in ['U_ERROR', 'V_ERROR', 'N_ERROR']:
                    self._sequence.append({'name': case, 'target': params["target"][0], 'p': fix,'type':'DEP2'})
                    self._noises.append(len(self._sequence) - 1)
                for case in ['X_ERROR', 'Y_ERROR', 'Z_ERROR']:
                    self._sequence.append({'name': case, 'target': params["target"][1], 'p': fix,'type':'DEP2'})
                    self._noises.append(len(self._sequence) - 1)
                for case in ['UX_ERROR', 'UY_ERROR', 'UZ_ERROR', 'VX_ERROR', 'VY_ERROR', 'VZ_ERROR', 'NX_ERROR', 'NY_ERROR', 'NZ_ERROR']:
                    self._sequence.append({'name': case, 'target': params["target"], 'p': fix,'type':'DEP2'})
                    self._noises.append(len(self._sequence) - 1)

            ##  用一个列表添加多个
            elif islist(params['target']) and isinteger(params['target'][0]) and len(params['target'])>2:
                for i in range(len(params['target'])//2):
                    self.append({'name': name, 'target':[params['target'][2*i], params['target'][2*i+1]], 'p': params['p']})

            ##  用多个列表添加多个
            elif islist(params['target']) and islist(params['target'][0]) and isinteger(params['target'][0][0]):
                for temp in params['target']:
                    self.append({'name': params['name'], 'target': temp, 'p': params["p"]})

            ##  其他情况抛出异常
            else:
                raise ValueError("Gate must be applied to two")

        ##  PART：----Add single qubit and single fermionic site measurement----
        elif name=='MZ' or name=='MN':
            target=params['target']

            ##  作用在一系列目标上
            if islist(target):
                for temp in target:
                    if 'p' in params:
                        dict_temp={'name': name, 'target': temp, 'p': params['p']}
                    else:
                        dict_temp={'name': name, 'target': temp}
                    self.append(dict_temp)

            ##  作用在整数index对应目标上
            elif isinteger(target):
                if name=='MZ':
                    if 'p' in params:
                        dict_temp={'name': 'MPP', 'target': PauliOperator([], [target], 1), 'p': params['p']}
                    else:
                        dict_temp={'name': 'MPP', 'target': PauliOperator([], [target], 1)}
                else:
                    if 'p' in params:
                        dict_temp={'name': 'MPP', 'target': MajoranaOperator([target], [target], 1j), 'p': params['p']}
                    else:
                        dict_temp={'name': 'MPP', 'target': MajoranaOperator([target], [target], 1j)}
                self.append(dict_temp)

        ##  PART：----Add qubit reset----
        elif name=='R':
            assert 'target' in params
            target=params['target']

            ##  作用在整数index对应目标上
            if isinteger(target):
                self._sequence.append({'name': name, 'target': target})
                self.sequence.append(self._sequence[-1])
                self.ideal_sequence.append(self._sequence[-1])
                if target==self.pauli_number:
                    self.pauli_number+=1
                elif target>self.pauli_number or target<0:
                    raise ValueError("R gate target must be consecutive")
                else:
                    pass

            ##  作用在一系列目标上
            elif islist(target):
                for temp in target:
                    self.append({'name': name, 'target': temp})

            ##  其他情况抛出异常
            else:
                raise ValueError

        ##  PART：----Add trap----
        elif name=='TRAP':
            assert len(self._sequence)==0
            assert 'pauli_number' in params
            assert 'majorana_number' in params
            self.pauli_number=params["pauli_number"]
            self.majorana_number=params["majorana_number"]
            self._sequence.append({'name': name})
            self.sequence.append(self._sequence[-1])
            self.ideal_sequence.append(self._sequence[-1])

        ##  PART：----Add string operator measurement----
        elif name=='MPP':
            assert 'target' in params

            ##  在列表上
            if islist(params['target']):
                for i, temp in enumerate(params['target']):
                    dict_temp={'name': 'MPP', 'target': temp}
                    if 'index' in params:
                        dict_temp['index']=params['index'][i]
                    self.append(dict_temp)

            ##  单个算符
            elif isinstance(params['target'], (PauliOperator, MajoranaOperator)):
                dict_temp={'name': 'MPP', 'target': params['target']}
                if 'index' in params:
                    dict_temp['index']=params['index']
                self._sequence.append(dict_temp)
                self._measurements.append(len(self._sequence)-1)
                if 'p' in params:
                    self._sequence.append({'name': 'M_ERROR', 'p': params["p"]})
                    self._noises.append(len(self._sequence)-1)
                    self.sequence.append({'name': 'MPP', 'target': params['target'], 'p': params["p"]})
                    self.ideal_sequence.append({'name': 'MPP', 'target': params['target']})
                else:
                    self.sequence.append(self._sequence[-1])
                    self.ideal_sequence.append({'name': 'MPP', 'target': params['target']})

            ##  其他情况抛出异常
            else:
                raise ValueError

        ##  PART：----Add detector----
        elif name=='DETECTOR':
            assert 'target' in params
            target=params['target']

            ##  对于负数的情况
            if all(target[i]=='negative' or target[i]<0 for i in range(len(target))):
                together=[len(self._measurements)+temp for temp in target if temp!='negative']  # 在测量中找到对应索引
                if 'negative' in target:
                    together.append('negative')
                self._detectors.append(together)

            ##  对于正数的情况
            elif all(target[i]=='negative' or target[i]>=0 for i in range(len(target))):
                together=[temp for temp in target if temp!='negative']
                if 'negative' in target:
                    together.append('negative')
                self._detectors.append(together)

            ##  抛出异常
            else:
                raise ValueError("DETECTOR gate target must be consecutive")

        ##  PART：----Add observable----
        elif name=='OBSERVABLE_INCLUDE':
            assert 'target' in params
            target=params['target']

            ##  对于负数的情况
            if all(target[i]=='negative' or target[i]<0 for i in range(len(target))):
                together=[len(self._measurements)+temp for temp in target if temp!='negative']  # 在测量中找到索引
                if 'negative' in target:
                    together.append('negative')
                self._observables.append(together)

            ##  对于正数的情况
            elif all(target[i]=='negative' or target[i]>=0 for i in range(len(target))):
                together=[temp for temp in target if temp!='negative']  # 在测量中找到索引
                if 'negative' in target:
                    together.append('negative')
                self._observables.append(together)

            ##  抛出异常
            else:
                raise ValueError("OBSERVABLE_INCLUDE gate target must be consecutive")

        ##  PART：----Add tick----
        elif name=='TICK':
            self.sequence.append({'name': 'TICK'})

        ##  PART：----Raise error----
        else:
            raise NotImplementedError


    #%%  SECTION：----Get batch sampling----
    def _sample_batch(self, noises=None, number=None,noiseless=False):
        if noises is not None:
            return [self.sample(i,noiseless) for i in noises]
        elif noises is None and isinteger(number):
            return [self.sample(noiseless=noiseless) for i in range(number)]
        else:
            raise ValueError

    ##  SECTION：----Get single sampling----
    def sample(self,designated_noise_index=None,noiseless=False):
        frame=Frame()
        reference_measurement_sample=self._reference_sample
        measurement_sample=np.empty(len(self._measurements), dtype=int)  # 生成测量值的样本数组
        flag_measurement=0

        ##  遍历整个操作序列
        for i, gate in enumerate(self._sequence):
            name=gate['name']

            if name in ['H','S','P']:
                target: int=gate['target']
                if name=='H':
                    frame.h(target)
                elif name=='S':
                    frame.s(target)
                elif name=='P':
                    frame.p(target)

            ##  执行双门
            if name in ['CX', 'CNX', 'BRAID', 'CNN']:
                target: list=gate['target']
                if name=='CX':
                    frame.cx(target[0], target[1])
                elif name=='CNX':
                    frame.cnx(target[0], target[1])
                elif name=='BRAID':
                    frame.braid(target[0], target[1])
                elif name=='CNN':
                    frame.cnn(target[0], target[1])

            ##  执行重置
            elif name=='R':
                target: int=gate['target']
                frame.reset(target)

            ##  执行误差门
            elif name =='M_ERROR':
                if designated_noise_index is None:
                    p: float=gate['p']
                else:
                    p=1.1
                if np.random.rand()<p:
                    measurement_sample[flag_measurement-1]=-measurement_sample[flag_measurement-1]

            elif name in ['X_ERROR', 'Y_ERROR', 'Z_ERROR', 'U_ERROR', 'V_ERROR', 'N_ERROR',
                          'XX_ERROR','XY_ERROR','XZ_ERROR','YX_ERROR','YY_ERROR','YZ_ERROR','ZX_ERROR','ZY_ERROR','ZZ_ERROR',
                          'UU_ERROR','UV_ERROR','UN_ERROR','VU_ERROR','VV_ERROR','VN_ERROR','NU_ERROR','NV_ERROR','NN_ERROR',
                          'UX_ERROR','UY_ERROR','UZ_ERROR','VX_ERROR','VY_ERROR','VZ_ERROR','NX_ERROR','NY_ERROR','NZ_ERROR']:
                if (designated_noise_index is None or i==designated_noise_index) and not noiseless:
                    if designated_noise_index is None:
                        p: float=gate['p']
                    else:
                        p=1.1
                    if isinteger(gate['target']):
                        frame.error(p,name,gate['target'])
                    elif islist(gate['target']) and len(gate['target'])==2:
                        frame.error(p,name,gate['target'][0],gate['target'][1])
                    else:
                        raise ValueError

            ##  执行测量
            elif name=='MPP':
                target=gate['target']
                measurement_sample[flag_measurement]=frame.measure(target, reference_measurement_sample[flag_measurement])
                flag_measurement+=1

            ##  执行初始化
            elif name=='TRAP':
                frame.trap(self.majorana_number, self.pauli_number)

        ##  PART：----Get detector sample----
        detector_sample=diff(measurement_sample, self._detectors)
        observable_sample=diff(measurement_sample, self._observables)
        return measurement_sample, detector_sample, observable_sample


    ##  SECTION：----Copy function----
    def copy(self):
        return copy.deepcopy(self)


##  CHAPTER：====Calculate detector sample====
def diff(measurement_sample, detectors):
    ##  计算探测器的结果
    detector_sample=np.empty(len(detectors), dtype=bool)
    flag_detector=0
    for i, detector in enumerate(detectors):
        value=1
        detector_sample[flag_detector]=False
        for temp in detector:
            if isinteger(temp):
                value=value*measurement_sample[temp]
            elif temp=='negative':
                value=-value
            else:
                raise NotImplementedError
        if value==1:
            detector_sample[flag_detector]=False
        elif value==-1:
            detector_sample[flag_detector]=True
        else:
            raise ValueError("测量结果中包含非元素")
        flag_detector+=1
    return detector_sample
