import numpy as np
from qiskit import QuantumCircuit
from exact_one_x_cx_ccx import exact_one_gate

# Implementation of the Fan_out_1 gate is based on Algorithm 1 from paper : Ancilla-free Quantum Adder with Sublinear Depth (https://arxiv.org/pdf/2501.16802)
# Authors : Maxime Remaud and Vivien Vandaele

def L1(n):
    qc =  QuantumCircuit(n)
    if n == 1:
        return qc
    if n == 2:
        qc.cx(0,1)
        return qc
    else :
        X = [1]
        CR = [[0,1]]
        CL = [[n-2,n-1]]
        for i in range(1,int(np.ceil(n/2))-1):
            CL.append([2*i-1,2*i])
            CR.append([2*i,2*i+1])
            X.append(2*i+1)
        if n%2 == 0:
            X.append(n-2)
        for k in CL:
                qc.cx(k[0],k[1])
        qc_mid = L1(len(X))
        qc.compose(qc_mid, qubits=X, inplace=True)
        for m in CR:
                qc.cx(m[0],m[1])
        return qc
      
def F1(n):
     qc = QuantumCircuit(n+1)
     qc1 = L1(n+1)
     qc2 = L1(n)
     qc.compose(qc1, qubits = qc.qubits, inplace=True)
     qc.compose(qc2.inverse(), qubits = list(range(1,n+1)), inplace=True)
     gate = qc.to_gate(label = f"F1_{n}")
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
    """Qubit order: [q0, q1, ..., q_{n-1}, anc, flag]"""
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
    fan_out_gate = F1(n)
    qc.append(fan_out_gate,[flag]+anc)
    for i in range(n):
        qc.ccz(anc[i], q[i], b[i])
    qc.append(fan_out_gate.inverse(),[flag]+anc)

    qc.append(prepare.inverse(), q + [anc[0], anc[1], flag])

    return qc
