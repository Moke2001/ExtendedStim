import numpy as np
import stim
import stimbposd
from extendedstim.Physics.PauliOperator import PauliOperator
from extendedstim.Physics.MajoranaOperator import MajoranaOperator
from extendedstim.Circuit.Circuit import Circuit


def pauli_steane_capacity_circuit():
    """""
    生成Pauli Steane码的现象学电路。
    """""
    circuit = Circuit()
    for i in range(7):
        circuit.append("R", i)
    stabilizers_x = [PauliOperator([3, 4, 5, 6], [], 1), PauliOperator([1, 2, 5, 6], [], 1), PauliOperator([0, 2, 4, 6], [], 1)]
    stabilizers_z = [PauliOperator([], [3, 4, 5, 6], 1), PauliOperator([], [1, 2, 5, 6], 1), PauliOperator([], [0, 2, 4, 6], 1)]

    for stabilizer in stabilizers_x:
        circuit.append("MPP", stabilizer)
    for stabilizer in stabilizers_z:
        circuit.append("MPP", stabilizer)

    for i in range(7):
        circuit.append("DEPOLARIZE1", i, 0.1)

    for stabilizer in stabilizers_x:
        circuit.append("MPP", stabilizer)
    for stabilizer in stabilizers_z:
        circuit.append("MPP", stabilizer)
    for i in range(6):
        circuit.append("DETECTOR", [-i - 1, -i - 7])

    circuit.append("MPP", PauliOperator([], [0, 1, 2], 1))
    circuit.append("OBSERVABLE_INCLUDE", [-1])

    stim_circuit = stim.Circuit()
    stim_circuit.append('R', range(7))
    stabilizers_x = [stim.PauliString('___XXXX'), stim.PauliString('_XX__XX'), stim.PauliString('X_X_X_X')]
    stabilizers_z = [stim.PauliString('___ZZZZ'), stim.PauliString('_ZZ__ZZ'), stim.PauliString('Z_Z_Z_Z')]
    for stabilizer in stabilizers_x:
        stim_circuit.append("MPP", stabilizer)
    for stabilizer in stabilizers_z:
        stim_circuit.append("MPP", stabilizer)
    stim_circuit.append("DEPOLARIZE1",range(7),0.1)
    for stabilizer in stabilizers_x:
        stim_circuit.append("MPP", stabilizer)
    for stabilizer in stabilizers_z:
        stim_circuit.append("MPP", stabilizer)
    for i in range(6):
        stim_circuit.append("DETECTOR", [stim.target_rec(-i - 1), stim.target_rec(-i -7)])
    stim_circuit.append("MPP",stim.PauliString('ZZZ____'))
    stim_circuit.append('OBSERVABLE_INCLUDE', [stim.target_rec(-1)],0)

    number_shot = 10_0000
    detector_error_model = stim_circuit.detector_error_model(decompose_errors=True)
    sampler = stim_circuit.compile_detector_sampler()
    detection_events, observable_flips = sampler.sample(shots=number_shot, separate_observables=True)
    predictions = stimbposd.bp_osd.BPOSD(model=detector_error_model, bp_method='minimum_sum').decode_batch(detection_events)
    num_errors = 0
    for shot in range(number_shot):
        actual_for_shot = observable_flips[shot]
        predicted_for_shot = predictions[shot]
        if not np.array_equal(actual_for_shot, predicted_for_shot):
            num_errors += 1
    print(circuit.sample(10_0000), 'extendedstim')
    print(num_errors / number_shot, 'Stim')


def pauli_steane_phenomenological_circuit():
    """""
    生成Pauli Steane码的现象学电路。
    """""
    circuit = Circuit()
    for i in range(7):
        circuit.append("R", i)
    stabilizers_x = [PauliOperator([3, 4, 5, 6], [], 1), PauliOperator([1, 2, 5, 6], [], 1), PauliOperator([0, 2, 4, 6], [], 1)]
    stabilizers_z = [PauliOperator([], [3, 4, 5, 6], 1), PauliOperator([], [1, 2, 5, 6], 1), PauliOperator([], [0, 2, 4, 6], 1)]

    for stabilizer in stabilizers_x:
        circuit.append("MPP", stabilizer)
    for stabilizer in stabilizers_z:
        circuit.append("MPP", stabilizer)

    for roundr in range(3):
        for i in range(7):
            circuit.append("DEPOLARIZE1", i, 0.01)

        for stabilizer in stabilizers_x:
            circuit.append("MPP", stabilizer,0.001)
        for stabilizer in stabilizers_z:
            circuit.append("MPP", stabilizer,0.001)
        for i in range(6):
            circuit.append("DETECTOR", [-i - 1, -i - 7])

    for stabilizer in stabilizers_x:
        circuit.append("MPP", stabilizer)
    for stabilizer in stabilizers_z:
        circuit.append("MPP", stabilizer)
    for i in range(6):
        circuit.append("DETECTOR", [-i - 1, -i - 7])

    circuit.append("MPP", PauliOperator([], [0, 1, 2], 1))
    circuit.append("OBSERVABLE_INCLUDE", [-1])

    ##  下面是Stim的情况
    stim_circuit = stim.Circuit()
    stim_circuit.append('R', range(7))
    stabilizers_x = [stim.PauliString('___XXXX'), stim.PauliString('_XX__XX'), stim.PauliString('X_X_X_X')]
    stabilizers_z = [stim.PauliString('___ZZZZ'), stim.PauliString('_ZZ__ZZ'), stim.PauliString('Z_Z_Z_Z')]
    for stabilizer in stabilizers_x:
        stim_circuit.append("MPP", stabilizer)
    for stabilizer in stabilizers_z:
        stim_circuit.append("MPP", stabilizer)
    for roundr in range(3):
        stim_circuit.append("DEPOLARIZE1", range(7), 0.01)
        for stabilizer in stabilizers_x:
            stim_circuit.append("MPP", stabilizer,0.001)
        for stabilizer in stabilizers_z:
            stim_circuit.append("MPP", stabilizer,0.001)
        for i in range(6):
            stim_circuit.append("DETECTOR", [stim.target_rec(-i - 1), stim.target_rec(-i - 7)])

    for stabilizer in stabilizers_x:
        stim_circuit.append("MPP", stabilizer)
    for stabilizer in stabilizers_z:
        stim_circuit.append("MPP", stabilizer)
    for i in range(6):
        stim_circuit.append("DETECTOR", [stim.target_rec(-i - 1), stim.target_rec(-i - 7)])

    stim_circuit.append("MPP", stim.PauliString('ZZZ____'))
    stim_circuit.append('OBSERVABLE_INCLUDE', [stim.target_rec(-1)], 0)

    number_shot = 100_0000
    detector_error_model = stim_circuit.detector_error_model(decompose_errors=True)
    sampler = stim_circuit.compile_detector_sampler()
    detection_events, observable_flips = sampler.sample(shots=number_shot, separate_observables=True)
    predictions = stimbposd.bp_osd.BPOSD(model=detector_error_model, bp_method='minimum_sum').decode_batch(detection_events)
    num_errors = 0
    for shot in range(number_shot):
        actual_for_shot = observable_flips[shot]
        predicted_for_shot = predictions[shot]
        if not np.array_equal(actual_for_shot, predicted_for_shot):
            num_errors += 1

    print(circuit.sample(number_shot), 'extendedstim')
    print(num_errors / number_shot, 'Stim')


def majorana_steane_phenomenological_circuit():
    """""
    生成Majorana Steane码的现象学电路。
    """""
    circuit = Circuit()
    for i in range(7):
        circuit.append("FR", i)
    stabilizers_x = [MajoranaOperator([3, 4, 5, 6], [], 1), MajoranaOperator([1, 2, 5, 6], [], 1), MajoranaOperator([0, 2, 4, 6], [], 1)]
    stabilizers_z = [MajoranaOperator([], [3, 4, 5, 6], 1), MajoranaOperator([], [1, 2, 5, 6], 1), MajoranaOperator([], [0, 2, 4, 6], 1)]

    for stabilizer in stabilizers_x:
        circuit.append("MPP", stabilizer)
    for stabilizer in stabilizers_z:
        circuit.append("MPP", stabilizer)

    for roundr in range(1):
        for i in range(7):
            circuit.append("FDEPOLARIZE1", i, 0.01)

        for stabilizer in stabilizers_x:
            circuit.append("MPP", stabilizer,0.001)
        for stabilizer in stabilizers_z:
            circuit.append("MPP", stabilizer,0.001)
        for i in range(6):
            circuit.append("DETECTOR", [-i - 1, -i - 7])

    for stabilizer in stabilizers_x:
        circuit.append("MPP", stabilizer)
    for stabilizer in stabilizers_z:
        circuit.append("MPP", stabilizer)
    for i in range(6):
        circuit.append("DETECTOR", [-i - 1, -i - 7])

    circuit.append("MPP", MajoranaOperator([0,1,2], [0, 1, 2], 1j))
    circuit.append("OBSERVABLE_INCLUDE", [-1])
    number_shot = 10_0000
    circuit.draw('test.svg')
    print(circuit.sample(number_shot), 'extendedstim')


if __name__ == '__main__':
    majorana_steane_phenomenological_circuit()