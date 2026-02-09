"""
Module: Code2Circuit
"""
from extendedstim.Code.BeyondCode.FermionicLatticeSurgery import FermionicLatticeSurgery
from extendedstim.Code.PrimitiveCode.MajoranaCode import MajoranaCode
from extendedstim.Circuit.Circuit import Circuit
from extendedstim.Code.PrimitiveCode.MajoranaCSSCode import MajoranaCSSCode
from extendedstim.Code.PrimitiveCode.PauliCSSCode import PauliCSSCode
from extendedstim.Code.PrimitiveCode.PauliCode import PauliCode
from extendedstim.Physics.MajoranaOperator import MajoranaOperator
from extendedstim.Physics.PauliOperator import PauliOperator
from extendedstim.tools.TypingTools import isinteger
import numpy as np


##  CHAPTER：====Details====


def Code2Circuit(code: MajoranaCode | PauliCode | FermionicLatticeSurgery, noise_model: str, cycle_number: int, pre_circuit: Circuit = None) -> Circuit | tuple[Circuit, Circuit]:
    """
    input.code: The quantum code to convert.
    input.noise_model: Noise model, options are 'phenomenological' or 'circuit-level'.
    input.cycle_number: Number of cycles.
    input.pre_circuit: Circuit to be appended before the main circuit.
    output: The converted circuit.
    influence: None
    """

    ##  SECTION：----Details----
    assert isinteger(cycle_number) and cycle_number >= 0
    assert isinstance(code, MajoranaCode) or isinstance(code, PauliCode) or isinstance(code, FermionicLatticeSurgery)

    ##  SECTION：----Details----
    if noise_model == 'phenomenological':
        if isinstance(code, MajoranaCSSCode):
            return _MajoranaCSSCode2PhenomenologicalCircuit(code, cycle_number)
        elif isinstance(code, PauliCSSCode):
            return _PauliCSSCode2PhenomenologicalCircuit(code, cycle_number)
        elif isinstance(code, MajoranaCode):
            raise NotImplementedError
        elif isinstance(code, PauliCode):
            raise NotImplementedError
        else:
            raise NotImplementedError

    elif noise_model == 'circuit-level':
        if isinstance(code, MajoranaCSSCode):
            return _MajoranaCSSCode2CircuitLevelCircuit(code, cycle_number, pre_circuit)
        elif isinstance(code, PauliCSSCode):
            return _PauliCSSCode2CircuitLevelCircuit(code, cycle_number)
        elif isinstance(code, MajoranaCode):
            raise NotImplementedError
        elif isinstance(code, PauliCode):
            raise NotImplementedError
        elif isinstance(code, FermionicLatticeSurgery):
            return _FermionicLatticeSurgery2CircuitLevelCircuit(code, cycle_number)
        else:
            raise NotImplementedError

    elif noise_model == 'code-capacity':
        if isinstance(code, MajoranaCSSCode):
            return _MajoranaCSSCode2CodeCapacityCircuit(code)
        else:
            raise NotImplementedError

    else:
        raise ValueError('noise_model must be phenomenological, circuit-level, or code-capacity')


def _MajoranaCSSCode2CodeCapacityCircuit(code: MajoranaCSSCode) -> Circuit:
    """
    input.code: A MajoranaCSSCode.
    output: The converted circuit.
    influence: None
    """

    ##  SECTION：----Details----
    stabilizers_x = code.generators_x

    ##  NZ
    stabilizers_n = []
    for i in range(len(code.generators_z)):
        stabilizers_n.append(code.generators_x[i] @ code.generators_z[i])

    ##  N
    logical_x = code.logical_operators_x
    logical_z = code.logical_operators_z
    logical_occupy = [1j * logical_x[temp] @ logical_z[temp] for temp in range(len(logical_x))]  # Particle number operator set as logical operator set

    majorana_number = code.physical_number  # Number of fermionic sites
    stabilizer_number = len(stabilizers_x) + len(stabilizers_n)  # Number of stabilizers

    ##  SECTION：----Details----
    circuit = Circuit()
    circuit.append({'name': 'TICK'})
    circuit.append({'name': 'TRAP', 'majorana_number': majorana_number, 'pauli_number': 0})
    circuit.append({'name': 'TICK'})

    observable_include = []  # Record indices of observables
    for i, logical_operator in enumerate(logical_occupy):
        circuit.append({"name": "MPP", "target": logical_operator})
        observable_include.append(len(circuit._measurements) - 1)

    ##  X
    for i, stabilizer in enumerate(stabilizers_x):
        circuit.append({"name": "MPP", "target": stabilizer})

    ##  N
    for i, stabilizer in enumerate(stabilizers_n):
        circuit.append({"name": "MPP", "target": stabilizer})

    for i in range(majorana_number):
        circuit.append({"name": "FDEPOLARIZE1", "target": i, "p": 0})

    ##  X
    for i, stabilizer in enumerate(stabilizers_x):
        circuit.append({"name": "MPP", "target": stabilizer, 'p': 0})

    ##  Z
    for i, stabilizer in enumerate(stabilizers_n):
        circuit.append({"name": "MPP", "target": stabilizer, 'p': 0})

    for i in range(stabilizer_number):
        circuit.append({"name": "DETECTOR", "target": [i - stabilizer_number, i - 2 * stabilizer_number]})

    for i, logical_operator in enumerate(logical_occupy):
        circuit.append({"name": "MPP", "target": logical_operator})
        circuit.append({"name": "OBSERVABLE_INCLUDE", "target": [len(circuit._measurements) - 1, observable_include[i]]})

    ##  SECTION：----Details----
    return circuit


##  CHAPTER：====Details====


def _MajoranaCSSCode2PhenomenologicalCircuit(code: MajoranaCSSCode, cycle_number: int) -> Circuit:
    """
    input.code: A MajoranaCSSCode.
    input.cycle_number: Number of syndrome measurements.
    output: The converted circuit.
    influence: None
    """

    ##  SECTION：----Details----
    stabilizers_x = code.generators_x

    ##  NZ
    stabilizers_n = []
    for i in range(len(code.generators_z)):
        stabilizers_n.append(code.generators_x[i] @ code.generators_z[i])

    ##  N
    logical_x = code.logical_operators_x
    logical_z = code.logical_operators_z
    logical_occupy = [1j * logical_x[temp] @ logical_z[temp] for temp in range(len(logical_x))]  # Particle number operator set as logical operator set

    majorana_number = code.physical_number  # Number of fermionic sites
    stabilizer_number = len(stabilizers_x) + len(stabilizers_n)  # Number of stabilizers

    ##  SECTION：----Details----
    circuit = Circuit()
    circuit.append({'name': 'TICK'})
    circuit.append({'name': 'TRAP', 'majorana_number': majorana_number, 'pauli_number': 0})
    circuit.append({'name': 'TICK'})

    observable_include = []  # Record indices of observables
    for i, logical_operator in enumerate(logical_occupy):
        circuit.append({"name": "MPP", "target": logical_operator})
        observable_include.append(len(circuit._measurements) - 1)

    ##  N
    for i, stabilizer in enumerate(stabilizers_n):
        circuit.append({"name": "MPP", "target": stabilizer})

    for _ in range(cycle_number):

        for i in range(majorana_number):
            circuit.append({"name": "FDEPOLARIZE1", "target": i, "p": 0})

        ##  X
        for i, stabilizer in enumerate(stabilizers_x):
            circuit.append({"name": "MPP", "target": stabilizer, 'p': 0})

        ##  Z
        for i, stabilizer in enumerate(stabilizers_n):
            circuit.append({"name": "MPP", "target": stabilizer, 'p': 0})

        if _ == 0:
            for i in range(stabilizer_number // 2):
                circuit.append({"name": "DETECTOR", "target": [i - stabilizer_number // 2, i - stabilizer_number - stabilizer_number // 2]})
        else:
            for i in range(stabilizer_number):
                circuit.append({"name": "DETECTOR", "target": [i - stabilizer_number, i - 2 * stabilizer_number]})

    for i in range(code.physical_number):
        circuit.append({"name": "MN", "target": i, 'p': 0})

    ##  N
    for i, stabilizer in enumerate(stabilizers_n):
        target = []  # Determine the N measurement and previous measurement to detect
        stabilizer_now = MajoranaOperator([], [], 1)  # Record the operator corresponding to the actual measurement product
        for index in stabilizer.occupy_x:
            target.append(-code.physical_number + index)
            stabilizer_now = stabilizer_now @ MajoranaOperator([index], [index], 1j)
        target.append(-code.physical_number - stabilizer_number // 2 + i)  # Add previous N measurement
        if stabilizer_now.coff != stabilizer.coff:
            target.append('negative')
        circuit.append({'name': 'DETECTOR', 'target': target})

    for i, logical_operator in enumerate(logical_occupy):
        circuit.append({"name": "MPP", "target": logical_operator})
        circuit.append({"name": "OBSERVABLE_INCLUDE", "target": [len(circuit._measurements) - 1, observable_include[i]]})

    ##  SECTION：----Details----
    return circuit


##  CHAPTER：====Details====
def _MajoranaCSSCode2CircuitLevelCircuit(code: MajoranaCSSCode, cycle_number: int, pre_circuit: Circuit = None) -> Circuit:
    """
    input.code: A MajoranaCSSCode.
    input.cycle_number: Number of syndrome measurements.
    input.pre_circuit: Circuit to be appended before the main circuit.
    output: The converted circuit.
    influence: None
    """

    ##  SECTION：----Details----
    stabilizers_x = code.generators_x
    stabilizers_z = code.generators_z
    stabilizers_n = []
    for i in range(len(code.generators_z)):
        stabilizers_n.append(code.generators_x[i] @ code.generators_z[i])

    logical_x = code.logical_operators_x
    logical_z = code.logical_operators_z
    logical_occupy = [1j * logical_x[temp] @ logical_z[temp] for temp in range(len(logical_x))]  # Particle number operator set as logical operator set

    majorana_number = code.physical_number
    stabilizer_number = len(stabilizers_x) + len(stabilizers_z)
    pauli_number = stabilizer_number

    ##  SECTION：----Details----
    circuit = Circuit()
    circuit.append({'name': 'TRAP', 'majorana_number': majorana_number, 'pauli_number': pauli_number})
    circuit.append({'name': 'TICK'})

    observable_include = []  # Record indices of observables
    for i, logical_operator in enumerate(logical_occupy):
        circuit.append({"name": "MPP", "target": logical_operator})
        observable_include.append(len(circuit._measurements) - 1)

    for i, stabilizer in enumerate(stabilizers_n):
        circuit.append({"name": "MPP", "target": stabilizer})
    circuit.append({'name': 'TICK'})

    for i in range(majorana_number):
        circuit.append({"name": "FDEPOLARIZE1", "target": i, "p": 0})
    for i in range(pauli_number):
        circuit.append({"name": "DEPOLARIZE1", "target": i, "p": 0})
    circuit.append({'name': 'TICK'})

    if pre_circuit is not None:
        for temp in pre_circuit._sequence:
            circuit.append(temp)
        circuit.append({'name': 'TICK'})

    for _ in range(cycle_number):
        for i, stabilizer in enumerate(stabilizers_x):
            sequence_temp = _syndrome_majorana_css_measurement_circuit(stabilizer, i, 'x')
            for temp in sequence_temp:
                circuit.append(temp)
            circuit.append({'name': 'TICK'})

        for i, stabilizer in enumerate(stabilizers_z):
            sequence_temp = _syndrome_majorana_css_measurement_circuit(stabilizer, i + len(stabilizers_x), 'z')
            for temp in sequence_temp:
                circuit.append(temp)
            circuit.append({'name': 'TICK'})

        ##  N
        if _ == 0:
            for i in range(stabilizer_number // 2):
                circuit.append({"name": "DETECTOR", "target": [i - stabilizer_number, i - stabilizer_number // 2, i - 3 * stabilizer_number // 2]})
        else:

            ##  X
            for i in range(stabilizer_number // 2):
                circuit.append({"name": "DETECTOR", "target": [i - stabilizer_number, i - 2 * stabilizer_number]})

            ##  N
            for i in range(stabilizer_number // 2):
                circuit.append({"name": "DETECTOR", "target": [i - stabilizer_number, i - stabilizer_number // 2, i - 3 * stabilizer_number // 2, i - 2 * stabilizer_number]})

    for i in range(code.physical_number):
        circuit.append({"name": "MN", "target": i, 'p': 0})

    for i, stabilizer in enumerate(stabilizers_n):
        target = []
        stabilizer_now = MajoranaOperator([], [], 1)
        for index in stabilizer.occupy_x:
            target.append(-code.physical_number + index)
            stabilizer_now = stabilizer_now @ MajoranaOperator([index], [index], 1j)
        target.append(-code.physical_number - stabilizer_number // 2 + i)
        target.append(-code.physical_number - stabilizer_number + i)
        if stabilizer_now.coff != stabilizer.coff:
            target.append('negative')
        circuit.append({'name': 'DETECTOR', 'target': target})

    for i, logical_operator in enumerate(logical_occupy):
        circuit.append({"name": "MPP", "target": logical_operator})
        circuit.append({"name": "OBSERVABLE_INCLUDE", "target": [len(circuit._measurements) - 1, observable_include[i]]})

    ##  SECTION：----Details----
    return circuit


##  CHAPTER：====Details====
def _FermionicLatticeSurgery2CircuitLevelCircuit(lattice_surgery: FermionicLatticeSurgery, cycle_number: int) -> Circuit:
    """
    input.lattice_surgery: A FermionicLatticeSurgeryCode.
    input.cycle_number: Number of syndrome measurements.
    output: The converted circuit.
    influence: None
    """

    ##  SECTION：----Details----
    irrelated_A_x = lattice_surgery.irrelated_stabilizers_A_x
    irrelated_A_z = lattice_surgery.irrelated_stabilizers_A_z
    irrelated_B_x = lattice_surgery.irrelated_stabilizers_B_x
    irrelated_B_z = lattice_surgery.irrelated_stabilizers_B_z
    related_A_z = lattice_surgery.related_stabilizers_A_z
    related_B_z = lattice_surgery.related_stabilizers_B_z
    modifiy_A_x = lattice_surgery.modify_stabilizers_A
    modifiy_B_x = lattice_surgery.modify_stabilizers_B
    gauge = lattice_surgery.gauge_stabilizers

    stabilizers_n = ([irrelated_A_x[temp] @ irrelated_A_z[temp] for temp in range(len(irrelated_A_x))] +
                     [irrelated_B_x[temp] @ irrelated_B_z[temp] for temp in range(len(irrelated_B_x))] +
                     [modifiy_A_x[temp] @ related_A_z[temp] for temp in range(len(modifiy_A_x))] +
                     [modifiy_B_x[temp] @ related_B_z[temp] for temp in range(len(modifiy_B_x))] +
                     gauge)
    stabilizers_z = irrelated_A_z + irrelated_B_z + related_A_z + related_B_z

    logicals = lattice_surgery.bare_logical_operators_x
    logical_occupy = [
        MajoranaOperator.HermitianOperatorFromOccupy(np.concatenate((temp.occupy_x, temp.occupy_z)), np.concatenate((temp.occupy_x, temp.occupy_z)))
        for temp in logicals]
    logical_occupy.append(MajoranaOperator.HermitianOperatorFromOccupy(lattice_surgery.remain_bare_logical_operator.occupy_z, lattice_surgery.remain_bare_logical_operator.occupy_z))

    majorana_number = lattice_surgery.physical_number
    stabilizer_number_z = len(stabilizers_z)
    stabilizer_number_n = len(stabilizers_n)
    pauli_number = stabilizer_number_z + stabilizer_number_n

    ##  SECTION：----Details----
    circuit = Circuit()
    circuit.append({'name': 'TRAP', 'majorana_number': majorana_number, 'pauli_number': pauli_number})
    circuit.append({'name': 'TICK'})

    observable_include = []  # Record indices of observables
    for i, logical_operator in enumerate(logical_occupy):
        circuit.append({"name": "MPP", "target": logical_operator})
        observable_include.append(len(circuit._measurements) - 1)

    for i in range(majorana_number):
        circuit.append({"name": "FDEPOLARIZE1", "target": i, "p": 0})
    for i in range(pauli_number):
        circuit.append({"name": "DEPOLARIZE1", "target": i, "p": 0})
    circuit.append({'name': 'TICK'})

    for _ in range(cycle_number):

        for i, stabilizer in enumerate(stabilizers_z):
            sequence_temp = _syndrome_majorana_css_measurement_circuit(stabilizer, i, 'z')
            for temp in sequence_temp:
                circuit.append(temp)
            circuit.append({'name': 'TICK'})

        for i, stabilizer in enumerate(stabilizers_n):
            sequence_temp = _syndrome_majorana_css_measurement_circuit(stabilizer, i + len(stabilizers_z), 'n')
            for temp in sequence_temp:
                circuit.append(temp)
            circuit.append({'name': 'TICK'})

        ##  N
        if _ == 0:
            for i in range(stabilizer_number_n):
                circuit.append({"name": "DETECTOR", "target": [i - stabilizer_number_n]})
        else:

            ##  X
            for i in range(stabilizer_number_z):
                circuit.append({"name": "DETECTOR",
                                "target": [i - stabilizer_number_n - stabilizer_number_z, i - 2 * (stabilizer_number_z + stabilizer_number_n)]})

            ##  N
            for i in range(stabilizer_number_n):
                circuit.append({"name": "DETECTOR", "target": [i - stabilizer_number_n, i - 2 * stabilizer_number_n - stabilizer_number_z]})

    for i in range(lattice_surgery.physical_number):
        circuit.append({"name": "MN", "target": i, 'p': 0})

    for i, stabilizer in enumerate(stabilizers_n):
        target = []
        for index in stabilizer.occupy_x:
            target.append(-lattice_surgery.physical_number + index)
        target.append(-lattice_surgery.physical_number - stabilizer_number_n + i)
        circuit.append({'name': 'DETECTOR', 'target': target})

    for i, logical_operator in enumerate(logical_occupy):
        circuit.append({"name": "MPP", "target": logical_operator})
        circuit.append({"name": "OBSERVABLE_INCLUDE", "target": [len(circuit._measurements) - 1, observable_include[i]]})

    ##  SECTION：----Details----
    return circuit


##  CHAPTER：====Details====
def _PauliCSSCode2PhenomenologicalCircuit(code: PauliCSSCode, cycle_number: int) -> tuple[Circuit, Circuit]:
    """
    input.code: A PauliCSSCode.
    input.cycle_number: Number of syndrome measurements.
    output: The converted circuit (x and z).
    influence: None
    """

    ##  SECTION：----Details----
    stabilizers_x = code.generators_x
    stabilizers_z = code.generators_z

    logical_x = code.logical_operators_x
    logical_z = code.logical_operators_z

    stabilizer_number = len(stabilizers_x) + len(stabilizers_z)
    pauli_number = code.physical_number

    ##  SECTION：----Details----
    circuit_x = Circuit()
    circuit_z = Circuit()
    circuit_x.append({'name': 'TRAP', 'majorana_number': 0, 'pauli_number': pauli_number})
    circuit_x.append({'name': 'H', 'target': range(pauli_number)})
    circuit_z.append({'name': 'TRAP', 'majorana_number': 0, 'pauli_number': pauli_number})

    for _ in range(cycle_number):

        for i in range(pauli_number):
            circuit_x.append({"name": "DEPOLARIZE1", "target": i, "p": 0})
            circuit_z.append({"name": "DEPOLARIZE1", "target": i, "p": 0})

        for i, stabilizer in enumerate(stabilizers_x):
            circuit_x.append({"name": "MPP", "target": stabilizer, 'p': 0})
            circuit_z.append({"name": "MPP", "target": stabilizer, 'p': 0})
        for i, stabilizer in enumerate(stabilizers_z):
            circuit_z.append({"name": "MPP", "target": stabilizer, 'p': 0})
            circuit_x.append({"name": "MPP", "target": stabilizer, 'p': 0})

        if _ == 0:
            for i in range(stabilizer_number // 2):
                circuit_z.append({"name": "DETECTOR", "target": [-stabilizer_number // 2 + i]})
                circuit_x.append({"name": "DETECTOR", "target": [-stabilizer_number + i]})
        else:
            for i in range(stabilizer_number):
                circuit_z.append({"name": "DETECTOR", "target": [-i - 1, -i - stabilizer_number - 1]})
                circuit_x.append({"name": "DETECTOR", "target": [-i - 1, -i - stabilizer_number - 1]})

    for i in range(code.physical_number):
        circuit_x.append({"name": "MPP", "target": PauliOperator([i], [], 1), 'p': 0})
        circuit_z.append({"name": "MPP", "target": PauliOperator([], [i], 1), 'p': 0})

    for i, stabilizer in enumerate(stabilizers_x):
        target = []
        for index in stabilizer.occupy_x:
            target.append(-code.physical_number + index)
        target.append(-code.physical_number - stabilizer_number + i)
        circuit_x.append({'name': 'DETECTOR', 'target': target})

    for i, stabilizer in enumerate(stabilizers_z):
        target = []
        for index in stabilizer.occupy_z:
            target.append(-code.physical_number + index)
        target.append(-code.physical_number - stabilizer_number // 2 + i)

        circuit_z.append({'name': 'DETECTOR', 'target': target})

    for i in range(len(logical_x)):
        circuit_z.append({"name": "MPP", "target": logical_z[i]})
        circuit_x.append({"name": "MPP", "target": logical_x[i]})
        circuit_z.append({"name": "OBSERVABLE_INCLUDE", "target": [len(circuit_z._measurements) - 1]})
        circuit_x.append({"name": "OBSERVABLE_INCLUDE", "target": [len(circuit_x._measurements) - 1]})

    ##  SECTION：----Details----
    return circuit_x, circuit_z


##  CHAPTER：====Details====
def _PauliCSSCode2CircuitLevelCircuit(code: PauliCSSCode, cycle_number: int) -> tuple[Circuit, Circuit]:
    """
    input.code: A PauliCSSCode.
    input.cycle_number: Number of syndrome measurements.
    output: The converted circuit (x and z).
    influence: None
    """

    ##  SECTION：----Details----
    stabilizers_x = code.generators_x
    stabilizers_z = code.generators_z

    logical_x = code.logical_operators_x
    logical_z = code.logical_operators_z

    stabilizer_number = len(stabilizers_x) + len(stabilizers_z)
    data_number = code.physical_number
    pauli_number = data_number + stabilizer_number

    ##  SECTION：----Details----
    circuit_x = Circuit()
    circuit_z = Circuit()
    circuit_x.append({'name': 'TRAP', 'majorana_number': 0, 'pauli_number': pauli_number})
    circuit_z.append({'name': 'TRAP', 'majorana_number': 0, 'pauli_number': pauli_number})
    circuit_z.append({'name': 'TICK'})
    circuit_x.append({'name': 'TICK'})
    circuit_x.append({'name': 'H', 'target': list(range(data_number))})

    for i in range(pauli_number):
        circuit_z.append({"name": "DEPOLARIZE1", "target": i, "p": 0})
        circuit_x.append({"name": "DEPOLARIZE1", "target": i, "p": 0})

    circuit_z.append({'name': 'TICK'})
    circuit_x.append({'name': 'TICK'})

    for _ in range(cycle_number):

        for i, stabilizer in enumerate(stabilizers_x):
            sequence_temp = _syndrome_pauli_css_measurement_circuit(stabilizer, i + data_number, 'x')
            for temp in sequence_temp:
                circuit_z.append(temp)
                circuit_x.append(temp)
            circuit_z.append({'name': 'TICK'})
            circuit_x.append({'name': 'TICK'})
        for i, stabilizer in enumerate(stabilizers_z):
            sequence_temp = _syndrome_pauli_css_measurement_circuit(stabilizer, i + data_number + len(stabilizers_x), 'z')
            for temp in sequence_temp:
                circuit_z.append(temp)
                circuit_x.append(temp)
            circuit_z.append({'name': 'TICK'})
            circuit_x.append({'name': 'TICK'})
        if _ == 0:
            for i in range(stabilizer_number // 2):
                circuit_z.append({"name": "DETECTOR", "target": [-stabilizer_number // 2 + i]})
                circuit_x.append({"name": "DETECTOR", "target": [-stabilizer_number + i]})
        else:
            for i in range(stabilizer_number):
                circuit_z.append({"name": "DETECTOR", "target": [-i - 1, -i - stabilizer_number - 1]})
                circuit_x.append({"name": "DETECTOR", "target": [-i - 1, -i - stabilizer_number - 1]})

    for i in range(code.physical_number):
        circuit_x.append({"name": "MPP", "target": PauliOperator([i], [], 1), 'p': 0})
        circuit_z.append({"name": "MPP", "target": PauliOperator([], [i], 1), 'p': 0})

    for i, stabilizer in enumerate(stabilizers_x):
        target = []
        for index in stabilizer.occupy_x:
            target.append(-code.physical_number + index)
        target.append(-code.physical_number - stabilizer_number + i)

        circuit_x.append({'name': 'DETECTOR', 'target': target})

    for i, stabilizer in enumerate(stabilizers_z):
        target = []
        for index in stabilizer.occupy_z:
            target.append(-code.physical_number + index)
        target.append(-code.physical_number - stabilizer_number // 2 + i)

        circuit_z.append({'name': 'DETECTOR', 'target': target})

    for i in range(len(logical_x)):
        circuit_z.append({"name": "MPP", "target": logical_z[i]})
        circuit_x.append({"name": "MPP", "target": logical_x[i]})
        circuit_z.append({"name": "OBSERVABLE_INCLUDE", "target": [len(circuit_z._measurements) - 1]})
        circuit_x.append({"name": "OBSERVABLE_INCLUDE", "target": [len(circuit_x._measurements) - 1]})

    ##  SECTION：----Details----
    return circuit_x, circuit_z


##  CHAPTER：====Details====
def _syndrome_majorana_css_measurement_circuit(stabilizer: MajoranaOperator, qubit_index: int, type: str) -> list[dict]:
    """
    input.stabilizer: A MajoranaOperator, represents stabilizer.
    input.qubit_index: An integer, represents the measured qubit index.
    input.type: A string, represents measurement type, can only be 'x', 'X', 'z', or 'Z'.
    output: Circuit sequence.
    influence: None
    """

    ##  SECTION：----Details----
    sequence = []  # Circuit sequence
    flag = True  # Gate type flag

    if type == 'x' or type == 'X':
        occupy = stabilizer.occupy_x
    elif type == 'z' or type == 'Z':
        occupy = stabilizer.occupy_z
    elif type == 'n' or type == 'N':
        assert np.all(stabilizer.occupy_z == stabilizer.occupy_x)
        occupy = stabilizer.occupy_x
    else:
        raise ValueError

    ##  SECTION：----Details----
    ##  N
    if type == 'n' or type == 'N':
        for j in range(len(occupy)):
            majorana_index_now = occupy[j]
            sequence.append({'name': 'CNX', 'target': [majorana_index_now, qubit_index], })
            sequence.append({'name': 'MIXTURE', 'target': [majorana_index_now, qubit_index], 'p': 0})

    ##  XZ
    elif type == 'x' or type == 'X' or type == 'z' or type == 'Z':

        ##  qubit
        sequence.append({'name': 'X', 'target': qubit_index})
        sequence.append({'name': 'DEPOLARIZE1', 'target': qubit_index, 'p': 0})

        for j in range(len(occupy)):
            majorana_index_now = occupy[j]

            ##  qubitCNX gate
            if j == len(occupy) - 1:
                sequence.append({'name': 'CNX', 'target': [majorana_index_now, qubit_index], })
                sequence.append({'name': 'MIXTURE', 'target': [majorana_index_now, qubit_index], 'p': 0})
                break

            majorana_index_down = occupy[j + 1]  # Next fermionic site

            ##  braid gate
            if flag:

                ##  braid
                if type == 'x' or type == 'X':
                    order_target = [majorana_index_down, majorana_index_now]
                elif type == 'z' or type == 'Z':
                    order_target = [majorana_index_now, majorana_index_down]
                else:
                    raise ValueError

                ##  braid gate
                sequence.append({"name": "BRAID", "target": order_target,})
                sequence.append({'name': 'FDEPOLARIZE2', 'target': order_target, 'p': 0})
                flag = False

            ##  CNN gate
            else:
                sequence.append({'name': 'CNN', 'target': [majorana_index_now, majorana_index_down], })
                sequence.append({'name': 'FDEPOLARIZE2', 'target': [majorana_index_now, majorana_index_down], 'p': 0})
                flag = True

        ##  syndrome extraction circuit
        flag = True
        for j in range(len(occupy) - 1):
            majorana_index_now = occupy[-1 - j]  # Current fermionic site
            majorana_index_up = occupy[-1 - j - 1]  # Previous fermionic site

            ##  braid gate
            if flag:
                if type == 'x' or type == 'X':
                    order_target = [majorana_index_now, majorana_index_up]
                elif type == 'z' or type == 'Z':
                    order_target = [majorana_index_up, majorana_index_now]
                else:
                    raise ValueError
                sequence.append({'name': 'N', 'target': [majorana_index_now]})
                sequence.append({'name': 'BRAID', 'target': order_target})
                sequence.append({'name': 'N', 'target': [majorana_index_now]})
                sequence.append({'name': 'FDEPOLARIZE2', 'target': order_target, 'p': 0})
                flag = False

            ##  CNN gate
            else:
                sequence.append({'name': 'CNN', 'target': [majorana_index_now, majorana_index_up]})
                sequence.append({'name': 'FDEPOLARIZE2', 'target': [majorana_index_now, majorana_index_up], 'p': 0})
                flag = True

    ##  qubit
    sequence.append({'name': 'MZ', 'target': qubit_index, 'p': 0})
    sequence.append({'name': 'R', 'target': qubit_index})
    sequence.append({'name': 'DEPOLARIZE1', 'target': qubit_index, 'p': 0})

    ##  SECTION：----Details----
    return sequence


##  CHAPTER：====Details====
def _syndrome_pauli_css_measurement_circuit(stabilizer: PauliOperator, qubit_index: int, type: str) -> list[dict]:
    """
    input.stabilizer: A PauliOperator, represents stabilizer.
    input.qubit_index: An integer, represents the measured qubit index.
    input.type: A string, represents measurement type, can only be 'x', 'X', 'z', or 'Z'.
    output: Circuit sequence.
    influence: None
    """

    ##  SECTION：----Details----
    sequence = []  # Circuit sequence

    if type == 'x' or type == 'X':
        occupy = stabilizer.occupy_x
    elif type == 'z' or type == 'Z':
        occupy = stabilizer.occupy_z
    else:
        raise ValueError

    ##  SECTION：----Details----
    if type == 'X' or type == 'x':
        sequence.append({'name': 'H', 'target': qubit_index})
        sequence.append({'name': 'DEPOLARIZE1', 'target': qubit_index, 'p': 0})
    for j in range(len(occupy)):
        if type == 'Z' or type == 'z':
            sequence.append({'name': 'CX', 'target': [occupy[j], qubit_index]})
            sequence.append({'name': 'DEPOLARIZE2', 'target': [occupy[j], qubit_index], 'p': 0})
        elif type == 'X' or type == 'x':
            sequence.append({'name': 'CX', 'target': [qubit_index, occupy[j]]})
            sequence.append({'name': 'DEPOLARIZE2', 'target': [occupy[j], qubit_index], 'p': 0})
        else:
            raise ValueError
    if type == 'X' or type == 'x':
        sequence.append({'name': 'H', 'target': qubit_index})
        sequence.append({'name': 'DEPOLARIZE1', 'target': qubit_index, 'p': 0})

    ##  qubit
    sequence.append({'name': 'MZ', 'target': qubit_index, 'p': 0})
    sequence.append({'name': 'R', 'target': qubit_index})
    sequence.append({'name': 'DEPOLARIZE1', 'target': qubit_index, 'p': 0})

    ##  SECTION：----Details----
    return sequence
