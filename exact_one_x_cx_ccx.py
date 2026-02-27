from qiskit import QuantumCircuit, QuantumRegister
from qiskit.compiler import transpile
from log_mcx_x_cx_ccx import mcx_log_gate
import numpy as np

print('The circuits are compiled in the [x, cx, ccx] basis.')

###### Excat-One Ladder ######

exact_one_ladder_memory = {}
exact_one_ladder_depth = {}
exact_one_ladder_size = {}

n = 1
qc = QuantumCircuit(2)
exact_one_ladder_memory[n] = qc.to_gate()
exact_one_ladder_depth[n] = qc.depth()
exact_one_ladder_size[n] = qc.size()

def exact_one_ladder_power_2(n):
    qc = QuantumCircuit(3*n-2)
    if n == 2:
        qc.cx(0,2)
        qc.cx(1,2)
        qc.ccx(0,1,3)
        return qc
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
        return qc
    
    
def exact_one_ladder(n, trace_depth_and_size=False):
    if n in exact_one_ladder_memory.keys() and not(trace_depth_and_size and n not in exact_one_ladder_depth.keys()):
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
            qc = transpile(qc, basis_gates=['x','cx','ccx'])
            exact_one_ladder_depth[n] = qc.depth()
            exact_one_ladder_size[n] = qc.size()
        return qc
    

###### Exact-One Oracle ######

exact_one_memory = {}
exact_one_depth = {}
exact_one_size = {}


exact_one_memory = {}
exact_one_depth = {}
exact_one_size = {}

n = 1
r = QuantumRegister(n, "r")
anc = QuantumRegister(3,"anc")
qc = QuantumCircuit(r, anc)
qc.cx(r[0], anc[2])
exact_one_memory[n] = qc.to_gate()
exact_one_depth[n] = qc.depth()
exact_one_size[n] = qc.size()


n = 2
r = QuantumRegister(n, "r")
anc = QuantumRegister(3,"anc")
qc = QuantumCircuit(r, anc)
qc.x(r[1])
qc.ccx(r[0], r[1], anc[2])
qc.x(r[1])
qc.x(r[0])
qc.ccx(r[0], r[1], anc[2])
qc.x(r[0])
exact_one_memory[n] = qc.to_gate()
exact_one_depth[n] = qc.depth()
exact_one_size[n] = qc.size()

def exact_one_gate(n, trace_depth_and_size=False):
    if n in exact_one_memory.keys() and not(trace_depth_and_size and n not in exact_one_depth.keys()):
        return exact_one_memory[n]

    k = n//3
    r = n%3
    r1 = QuantumRegister(k + 1*(r != 0), "r1")
    r2 = QuantumRegister (k + r//2, "r2")
    r3 = QuantumRegister (k , "r3")
    anc = QuantumRegister(3,"anc")
    qregs = [r1,
            r2,
            r3,
            anc,
    ]
    qc = QuantumCircuit(*qregs)

    regs_permut = [r1, r2, r3]
    
    for i in range(3):
        r1_bis, r2_bis, r3_bis = regs_permut
        
        reg_anc = [i for i in r2_bis] + [i for i in r3_bis]
        mcx_gate = mcx_log_gate(len(reg_anc))
        for i in reg_anc:
            qc.x(i)
        qc.append(mcx_gate, qargs = reg_anc + [anc[0], anc[2], anc[1]])
        for i in reg_anc:
            qc.x(i)

        reg_ctrl = [i for i in r1_bis]
        qc_r1 = exact_one_ladder(len(reg_ctrl))
        size_anc = qc_r1.num_qubits - len(reg_ctrl)
        if size_anc == 1 : #Can be True only for 3 <= n < 6
            reg_hamming = reg_ctrl + [anc[0]]
        else:
            reg_hamming = reg_ctrl + reg_anc[:size_anc-2] + [anc[0]] + [reg_anc[size_anc-2]]
        qc.compose(qc_r1, qubits = reg_hamming, inplace=True)
        qc.ccx(reg_hamming[-2],anc[1], anc[2]) #If size_anc == 1, reg_hamming[-2] == reg_ctrl[0], else reg_hamming[-2] == anc[0] 
        qc.compose(qc_r1.inverse(), qubits = reg_hamming, inplace=True)
        
        for i in reg_anc:
            qc.x(i)
        qc.append(mcx_gate, qargs = reg_anc + [anc[0], anc[2], anc[1]])
        for i in reg_anc:
            qc.x(i)

        regs_permut = regs_permut[1:] + [regs_permut[0]]
    
    gate = qc.to_gate()
    exact_one_memory[n] = gate
    if trace_depth_and_size:
        qc = transpile(qc, basis_gates=['x','cx','ccx'])
        exact_one_depth[n] = qc.depth()
        exact_one_size[n] = qc.size()

    return gate

def access_exact_one_depth(nc):
    if nc in exact_one_depth.keys():
        return exact_one_depth[nc]
    else:
        print('Such gate has not been compiled yet.')

def access_exact_one_size(nc):
    if nc in exact_one_size.keys():
        return exact_one_size[nc]
    else:
        print('Such gate has not been compiled yet.')

def exact_one_circuit(n):
    nqubits = n+3
    exact_one = exact_one_gate(n, trace_depth_and_size=True)
    log_qc = QuantumCircuit(nqubits)
    log_qc.append(exact_one, list(range(nqubits)))
    log_qc = transpile(log_qc, basis_gates=['x', 'cx', 'ccx'])
    return log_qc