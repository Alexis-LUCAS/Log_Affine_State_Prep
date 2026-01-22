from qiskit import QuantumCircuit, QuantumRegister
from qiskit.compiler import transpile
from log_mcx_x_cx_ccx import mcx_log_gate
from exact_one_ladder import exact_one_ladder

exact_one_memory = {}
exact_one_depth = {}
exact_one_size = {}


exact_one_memory = {}
exact_one_depth = {}
exact_one_size = {}

print('The circuits are compiled in the [x, cx, ccx] basis.')

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