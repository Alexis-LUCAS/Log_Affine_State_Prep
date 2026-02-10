import numpy as np
from qiskit import QuantumCircuit
from exact_one_x_cx_ccx import exact_one_gate

### FAN-OUT GATE

def Fan_out(n):
    qc = QuantumCircuit(n+1)
    for i in range(np.floor(np.log2(n+1)).astype(int)+1):
        for j in range(2**i):
            if 2**i + j < n+1:
                qc.cx(j, 2**i + j)
    gate = qc.to_gate(label = f"Fan_out_{n}")
    return gate

### SPUE CIRCUIT CONSTRUCTION

# A_n PREPARE
def A_n(n):
    qc = QuantumCircuit(n, name=f"A_{n}")
    for k in range(1, n + 1):
        th = np.arcsin(np.sqrt(1 / (2**k + 1)))
        qc.ry(2 * th, n - k)
    return qc.to_gate(label=f"A_{n}")


# PREPARE + EXACT-1
def build_circuit_from_A(n):
    """Qubit order: [q0, q1, ..., q_{n-1}, anc1, anc2, flag]"""
    qc = QuantumCircuit(n + 3, name=f"PREPARE_{n}")
    qc.append(A_n(n), list(range(n)))
    qc.append(exact_one_gate(n), qc.qubits)
    return qc

# PREPARE – SELECT – PREPARE†
def build_prepare_select_prepare_dag(n):
    """
    Qubit layout (Qiskit bitstring, q0 on the right):
    [anc_(n-1) ... anc_0 | b_(n-1) ... b_0 | flag | q_(n-1) ... q_0 ]
    """
    if n == 1: # 0 ancilla
         qc = QuantumCircuit(3, name="PREPARE_SELECT_PREPAREdag")
         qc.h(2)
         qc.append(A_n(1), qargs = [0])
         qc.cx(0,1)
         qc.ccz(1,0,2)
         qc.cx(0,1)
         qc.append(A_n(1).inverse(), qargs = [0])
         return qc
    
    qc = QuantumCircuit(3*n + 1, name="PREPARE_SELECT_PREPAREdag")

    q    = list(range(n))
    flag = n
    b    = list(range(n + 1, 2 * n + 1))
    anc = list(range(2 * n + 1, 3 * n + 1))

    # Initialize b ancillas in |+>
    for bi in b:
        qc.h(bi)

    prepare = build_circuit_from_A(n).to_gate(label=f"PREPARE_{n}")
    qc.append(prepare, q + [anc[0], anc[1], flag])

    # SELECT
    fan_out_gate = Fan_out(n)
    qc.append(fan_out_gate,[flag]+anc)
    for i in range(n):
        qc.ccz(anc[i], q[i], b[i])
    qc.append(fan_out_gate.inverse(),[flag]+anc)

    qc.append(prepare.inverse(), q + [anc[0], anc[1], flag])

    return qc
