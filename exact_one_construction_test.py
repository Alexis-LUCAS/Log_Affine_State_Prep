"""
The code below is based on construction test from https://github.com/BaptisteClaudon/Polylog_MCX-public.
"""

import collections
import itertools
import random
import numpy as np
import pytest
from qiskit import QuantumCircuit
from qiskit.compiler import transpile
from exact_one_x_cx_ccx import exact_one_gate

def bulk_simulate_result_of_applying_classical_circuit_to(circuit, input_states: np.ndarray) -> np.ndarray:
    states = np.copy(input_states)
    buffer = np.zeros_like(states[0])
    for instruction in circuit:
        if instruction.operation.name == 'x':
            assert len(instruction.qubits) == 1
            q = instruction.qubits[0]._index
            np.bitwise_not(states[q], out=states[q])
        elif instruction.operation.name == 'cx':
            assert len(instruction.qubits) == 2
            c = instruction.qubits[0]._index
            t = instruction.qubits[1]._index
            states[t] ^= states[c]
        elif instruction.operation.name == 'ccx':
            assert len(instruction.qubits) == 3
            a = instruction.qubits[0]._index
            b = instruction.qubits[1]._index
            t = instruction.qubits[2]._index
            np.bitwise_and(states[a], states[b], out=buffer)
            states[t] ^= buffer
        elif instruction.operation.name == 'barrier':
            pass
        elif instruction.operation.name == 'measure':
            pass
        else:
            raise NotImplementedError(f'{instruction=}')
    return states

def compute_depth(circuit) -> int:
    qubit_depth = collections.defaultdict(int)
    known = ['x', 'cx', 'ccx']
    for instruction in circuit:
        if instruction.operation.name in known:
            qs = [q._index for q in instruction.qubits]
            layer = max(qubit_depth[q] for q in qs) + 1
            for q in qs:
                qubit_depth[q] = layer
        elif instruction.operation.name == 'barrier':
            pass
        elif instruction.operation.name == 'measure':
            pass
        else:
            raise NotImplementedError(f'{instruction=}')
    return max(qubit_depth.values())

@pytest.mark.parametrize('bit_string_size', [10, 100, 200, 500, 1000])
def test_depth(bit_string_size: int):
    num_qubits = bit_string_size + 3
    gate = exact_one_gate(bit_string_size)
    circuit = QuantumCircuit(num_qubits)
    circuit.append(gate, list(range(num_qubits)))
    circuit = transpile(circuit, basis_gates=['x', 'cx', 'ccx'])
    depth = compute_depth(circuit)
    lg_n = bit_string_size.bit_length()
    expected = 634 * lg_n # Heuristic upper bound
    assert depth <= expected

@pytest.mark.parametrize('bit_string_size', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 32, 64, 100])
def test_fuzz_random_cases(bit_string_size: int):
    num_samples = 1024
    num_qubits = bit_string_size + 3
    target_index = bit_string_size + 2
    gate = exact_one_gate(bit_string_size)
    circuit = QuantumCircuit(num_qubits)
    circuit.append(gate, list(range(num_qubits)))
    circuit = transpile(circuit, basis_gates=['x', 'cx', 'ccx'])
    input_states = np.random.randint(
        low=0,
        high=(1 << 64) - 1,
        size=(num_qubits, num_samples // 64),
        dtype=np.uint64,
    )
    input_states[bit_string_size] = np.zeros_like(input_states[bit_string_size])
    input_states[bit_string_size+1] = np.zeros_like(input_states[bit_string_size+1])
    exact_one = np.zeros_like(input_states[0])
    for i in range(bit_string_size):
        term = input_states[i]
        for j in range(bit_string_size):
            if j != i:
                term &= ~input_states[j]
        exact_one |= term
    expected = np.copy(input_states)
    expected[target_index] ^= exact_one
    actual = bulk_simulate_result_of_applying_classical_circuit_to(
        circuit, input_states
    )
    assert np.array_equal(actual, expected)


@pytest.mark.parametrize('bit_string_size', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 32, 64, 100])
@pytest.mark.parametrize('num_ones_qubits', [0, 1, 2])
def test_exact_ones_cases(bit_string_size: int, num_ones_qubits: int):
    num_qubits = bit_string_size + 3
    target_index = bit_string_size + 2
    ancilla_index_0 = bit_string_size
    ancilla_index_1 = bit_string_size + 1
    gate = exact_one_gate(bit_string_size)
    circuit = QuantumCircuit(num_qubits)
    circuit.append(gate, list(range(num_qubits)))
    circuit = transpile(circuit, basis_gates=['x', 'cx', 'ccx'])
    cases = list(itertools.combinations(range(bit_string_size), num_ones_qubits))
    num_samples = len(cases)
    input_states = np.zeros((num_qubits, num_samples), dtype=np.bool_)
    for shot_index, on_bits in enumerate(cases):
        for q in on_bits:
            input_states[q, shot_index] = True
        input_states[target_index, shot_index] = random.random() < 0.5
        input_states[ancilla_index_0, shot_index] = False
        input_states[ancilla_index_1, shot_index] = False
    input_states_packed = np.packbits(input_states, axis=1)
    exact_one = np.zeros_like(input_states_packed[0])
    for i in range(bit_string_size):
        term = input_states_packed[i]
        for j in range(bit_string_size):
            if j != i:
                term &= ~input_states_packed[j]
        exact_one |= term
    expected = np.copy(input_states_packed)
    expected[target_index] ^= exact_one
    actual = bulk_simulate_result_of_applying_classical_circuit_to(
        circuit, input_states_packed
    )
    assert np.array_equal(actual, expected)