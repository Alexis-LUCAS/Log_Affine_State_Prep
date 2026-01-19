from qiskit import QuantumCircuit
from qiskit.compiler import transpile
import numpy as np

exact_one_ladder_memory = {}
exact_one_ladder_depth = {}
exact_one_ladder_size = {}

print('The circuits are compiled in the [x, cx, ccx] basis.')

n = 1
qc = QuantumCircuit(2)
exact_one_ladder_memory[n] = qc.to_gate()
exact_one_ladder_depth[n] = qc.depth()
exact_one_ladder_size[n] = qc.size()

def exact_one_ladder_power_2(n):
    qc = QuantumCircuit(3*n-2)
    if n == 2:
        qc.x(1)
        qc.ccx(0,1,2)
        qc.x(1)
        qc.x(0)
        qc.ccx(0,1,2)
        qc.x(0)
        qc.ccx(0,1,3)
        gate = qc.to_gate()
        return gate
    else:
        qc1 = exact_one_ladder_power_2(n//2)
        reg1 = [k for k in range(0,n//2)] + [k for k in range(n,n +2*(n//2)-2)]
        reg2 = [k for k in range(n//2, n)] + [k for k in range(n +2*(n//2)-2, 3*n -4)]
        qc.compose(qc1, qubits = reg1, inplace=True)
        qc.compose(qc1, qubits = reg2, inplace=True)
        qc.x(reg1[-1])
        qc.x(reg2[-1])
        qc.x(3*n-3)
        qc.ccx(reg1[-2],reg2[-1],3*n-4)
        qc.ccx(reg1[-1],reg2[-2],3*n-4)
        qc.ccx(reg1[-1],reg2[-1],3*n-3)
        qc.ccx(reg1[-2],reg2[-2],3*n-3)
        qc.x(reg1[-1])
        qc.x(reg2[-1])
        gate = qc.to_gate()
        return gate
    
    
def exact_one_ladder(n, trace_depth_and_size=False):
    if n in exact_one_ladder_memory:
        return exact_one_ladder_memory[n]
    k = int(np.floor(np.log2(n)))
    if n == 2**k:
        qc = exact_one_ladder_power_2(2**k)
        exact_one_ladder_memory[n] = qc
        if trace_depth_and_size:
            exact_one_ladder_depth[n] = qc.depth()
            exact_one_ladder_size[n] = qc.size()
        return qc
    else:
        qc1 = exact_one_ladder_power_2(2**k)
        qc2 = exact_one_ladder(n-2**k)
        nb_ancilla_1 = qc1.num_qubits - 2**k
        nb_ancilla_2 = qc2.num_qubits - (n - 2**k)
        reg1 = [k for k in range(0,2**k)] + [k for k in range(n,n +nb_ancilla_1)]
        reg2 = [k for k in range(2**k, n)] + [k for k in range(n +nb_ancilla_1, n + nb_ancilla_1 + nb_ancilla_2)]
        qc_size = n + nb_ancilla_1 + nb_ancilla_2 + 2
        qc = QuantumCircuit(qc_size)
        qc.compose(qc1, qubits = reg1, inplace=True)
        qc.compose(qc2, qubits = reg2, inplace=True)
        qc.x(reg1[-1])
        qc.x(reg2[-1])
        qc.x(qc_size-1)
        qc.ccx(reg1[-2],reg2[-1],qc_size-2)
        qc.ccx(reg1[-1],reg2[-2],qc_size-2)
        qc.ccx(reg1[-1],reg2[-1],qc_size-1)
        qc.ccx(reg1[-2],reg2[-2],qc_size-1)
        qc.x(reg1[-1])
        qc.x(reg2[-1])
        exact_one_ladder_memory[n] = qc
        if trace_depth_and_size:
            exact_one_ladder_depth[n] = qc.depth()
            exact_one_ladder_size[n] = qc.size()
        return qc