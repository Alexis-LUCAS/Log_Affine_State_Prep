"""
The following code implementing multi-control-X gates (MCX) with logarithmic depth is based on paper : 
Rise of conditionally clean ancillae for efficient quantum circuit constructions (https://quantum-journal.org/papers/q-2025-05-21-1752/pdf/)
Authors : Tanuj Khattar and Craig Gidney
"""

from qiskit import QuantumCircuit
from qiskit.compiler import transpile

linear_mcx_memory = {}
linear_mcx_depth = {}
linear_mcx_size = {}


log_mcx_memory = {}
log_mcx_depth = {}
log_mcx_size = {}

print('The circuits are compiled in the [x, cx, ccx] basis.')

nc = 1
qc = QuantumCircuit(nc+2)
qc.cx(control_qubit=0, target_qubit=nc+1)
linear_mcx_memory[nc] = qc.to_gate()
linear_mcx_depth[nc] = qc.depth()
linear_mcx_size[nc] = qc.size()

nc = 1
qc = QuantumCircuit(nc+3)
qc.cx(control_qubit=0, target_qubit=nc+2)
log_mcx_memory[nc] = qc.to_gate()
log_mcx_depth[nc] = qc.depth()
log_mcx_size[nc] = qc.size()



nc = 2
qc = QuantumCircuit(nc+2)
qc.ccx(control_qubit1=0, control_qubit2=1, target_qubit=nc+1)
linear_mcx_memory[nc] = qc.to_gate()
linear_mcx_depth[nc] = qc.depth()
linear_mcx_size[nc] = qc.size()

nc = 2
qc = QuantumCircuit(nc+3)
qc.ccx(control_qubit1=0, control_qubit2=1, target_qubit=nc+2)
log_mcx_memory[nc] = qc.to_gate()
log_mcx_depth[nc] = qc.depth()
log_mcx_size[nc] = qc.size()



def get_linear_depth_ladder_ops(ncontrol):
    #Initialize ladder ops
    ladder_ops = []
    ncontrol +=1 # Ancilla included

    #Up-ladder
    for i in range(0, ncontrol - 2, 2):
        x, y, t = i + 1, i + 2, i
        ladder_ops.append((x, y, t))

    #Down-ladder

    if ncontrol & 1:
        x = ncontrol - 3
        y = ncontrol - 5
        t = ncontrol - 6
    else:
        x = ncontrol - 1
        y = ncontrol - 4
        t = ncontrol - 5
    if t > 0:
        ladder_ops.append((x, y, t))


    for i in range(t, 2, -2):
        ladder_ops.append((i, i - 1, i - 2))

    #Final controls

    mid_second_ctrl = 1 + max(0, 6 - ncontrol)

    return ladder_ops, mid_second_ctrl



def mcx_linear_gate(ncontrol, clean=False, trace_depth_and_size=False):
    """Multi-control-single-target $C^{n}X$ using 1 dirty ancilla (or clean when clean=True) and depth n. Register is of form [controls...,ancilla, target]"""

    if ncontrol in linear_mcx_memory.keys() and not(trace_depth_and_size and ncontrol not in linear_mcx_depth.keys()):
        return linear_mcx_memory[ncontrol]

    nqubits = ncontrol + 2
    controls = list(range(ncontrol))
    anc = ncontrol
    target = ncontrol + 1

    ladder_ops, mid_second_ctrl = get_linear_depth_ladder_ops(ncontrol)

    qc = QuantumCircuit(nqubits)
    
    #Ladder_ops
    for idx, (x, y, t) in enumerate(ladder_ops):
        if t == 0:
            qc.ccx(x, y, t)
        else:
            if idx <= (len(ladder_ops)//2):
                qc.ccx(x, y, t)
                qc.x(t)
            else :
                qc.x(t)
                qc.ccx(x, y, t)

    #Final controls
    qc.ccx(0, mid_second_ctrl, target)

    #Reverse Ladder_ops
    for idx, (x, y, t) in enumerate(reversed(ladder_ops)):
        if t == 0:
            qc.ccx(x, y, t)
        else:
            if idx <= (len(ladder_ops)//2) -1 : # Trick to reduce depth manually (qiskit does not optimize well here)
                # Toffoli + X
                qc.ccx(x, y, t)
                qc.x(t)
            else :
                qc.x(t)
                qc.ccx(x, y, t)
    
    if not clean:
        #Ladder_ops
        for idx, (x, y, t) in enumerate(ladder_ops):
            if t != 0:
                if idx <= (len(ladder_ops)//2):
                    # Toffoli + X
                    qc.ccx(x, y, t)
                    qc.x(t)
                else :
                    qc.x(t)
                    qc.ccx(x, y, t)
        #Final controls
        qc.ccx(0, mid_second_ctrl, target)

        #Reverse_Ladder_ops
        for idx, (x, y, t) in enumerate(reversed(ladder_ops)):
            if t != 0:
                if idx <= (len(ladder_ops)//2) -1 :
                    qc.ccx(x, y, t)
                    qc.x(t)
                else :
                    qc.x(t)
                    qc.ccx(x, y, t)
            
        
    qc_reordered = QuantumCircuit(nqubits)
    qc_reordered.compose(qc, qubits = [anc] + controls + [target], inplace=True)
    gate = qc_reordered.to_gate()

    if not clean : #Storing gates only if ancilla is dirty (what we use in final mcx decomposition)
        linear_mcx_memory[ncontrol] = gate
        if trace_depth_and_size:
            linear_mcx_depth[ncontrol] = qc_reordered.depth()
            linear_mcx_size[ncontrol] = qc_reordered.size()

    return gate


def get_log_depth_ladder_ops(ncontrol):

    ancilla = 0
    ctrls = list(range(1, ncontrol + 1))

    anc = [ancilla]         
    ladder_ops = []
    final_ctrls = []

    while len(ctrls) > 1:
        next_batch_len = min(len(anc) + 1, len(ctrls))
        batch = ctrls[:next_batch_len]
        ctrls = ctrls[next_batch_len:]

        new_anc = []

        while len(batch) > 1:
            ccx_n = len(batch) // 2
            st = len(batch) & 1  

            xs = batch[st : st + ccx_n]
            ys = batch[st + ccx_n : st + 2 * ccx_n]
            ts = anc[-ccx_n:]

            assert len(xs) == len(ys) == len(ts) == ccx_n >= 1

            ladder_ops.append((xs, ys, ts))

            new_anc += batch[st:]
            batch = ts + batch[:st]
            anc = anc[:-ccx_n]

        anc = sorted(anc + new_anc)
        final_ctrls += batch

    final_ctrls += ctrls
    return ladder_ops, sorted(final_ctrls)


def mcx_log_gate(ncontrol, clean=False, trace_depth_and_size=False):
    """Multi-control-single-target $C^{n}X$ using 2 dirty ancillae (or clean when clean=True) and depth log(n). Register is of form [controls..., ancilla1, ancilla2, target]"""

    if ncontrol in log_mcx_memory.keys() and not(trace_depth_and_size and ncontrol not in log_mcx_depth.keys()):
        return log_mcx_memory[ncontrol]

    nqubits = ncontrol + 3
    controls = list(range(ncontrol))
    anc1 = ncontrol
    anc2 = ncontrol + 1
    target = ncontrol + 2

    ladder_ops, final_ctrls = get_log_depth_ladder_ops(ncontrol)

    qc = QuantumCircuit(nqubits)
    
    #Ladder_ops
    for (x, y, t) in ladder_ops:
        if t[0] == 0:
            qc.ccx(x, y, t)
        else:
            qc.ccx(x, y, t)
            qc.x(t)

    #Central linear MCX
    linear_mcx = mcx_linear_gate(len(final_ctrls))
    qc.append(linear_mcx, final_ctrls + [anc2, target])

    #Reverse_Ladder_ops
    for (x, y, t) in reversed(ladder_ops):
        if t[0] == 0:
            qc.ccx(x, y, t)
        else:
            qc.x(t)
            qc.ccx(x, y, t)

    if not clean:
        #Ladder_ops
        for (x, y, t) in ladder_ops[1:]:
                qc.ccx(x, y, t)
                qc.x(t)

        #Central linear MCX
        linear_mcx = mcx_linear_gate(len(final_ctrls))
        qc.append(linear_mcx, final_ctrls + [anc2, target])
        #Reverse_Ladder_ops
        for (x, y, t) in reversed(ladder_ops[1:]):
                qc.x(t)
                qc.ccx(x, y, t)
        
    qc_reordered = QuantumCircuit(nqubits)
    qc_reordered.compose(qc, qubits = [anc1] + controls + [anc2, target], inplace=True)
    gate = qc_reordered.to_gate()

    if not clean : #Storing gates only if ancilla is dirty
        log_mcx_memory[ncontrol] = gate
        if trace_depth_and_size:
            qc_reordered = transpile(qc_reordered, basis_gates=['x','cx','ccx'])
            log_mcx_depth[ncontrol] = qc_reordered.depth()
            log_mcx_size[ncontrol] = qc_reordered.size()

    return gate

def access_mcx_log_depth(nc):
    if nc in log_mcx_depth.keys():
        return log_mcx_depth[nc]
    else:
        print('Such gate has not been compiled yet.')

def access_mcx_size(nc):
    if nc in log_mcx_size.keys():
        return log_mcx_size[nc]
    else:
        print('Such gate has not been compiled yet.')

def mcx_log_circuit(ncontrol):
    nqubits = ncontrol+3
    log_gate = mcx_log_gate(ncontrol, trace_depth_and_size=True)
    log_qc = QuantumCircuit(nqubits)
    log_qc.append(log_gate, list(range(nqubits)))
    log_qc = transpile(log_qc, basis_gates=['x', 'cx', 'ccx'])
    return log_qc