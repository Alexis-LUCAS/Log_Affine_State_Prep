# Log_Affine_State_Prep


This repository contains the code used in the paper "Logarithmic-depth quantum state preparation of polynomials" (Lien Arxiv) to prepare the operator $1-2p_n L_n$.

In particular, it contains:

* a python file implementing MCX gates (Khattar and Gidney's method with clean zeroed ancillae and dirty ancillae) with Qiskit, entitled log_mcx_x_cx_ccx.py, and a test file log_mcx_construction_test.py (to run from the terminal with the command 'pytest log_mcx_construction_test.py'),
* a python file implementing the Exact-one oracle, entitled exact_one_x_cx_ccx.py, and a test file exact_one_construction_test.py (to run from the terminal with the command 'pytest exact_one_construction_test.py'),
* a python file to prepare the operator $1 - 2 p_n L_n$ , entitled spue_circuit.py,
* a Notebook simulating the state preparation of  $(1 - 2 p_n L_n)H^{\otimes n} \ket 0^{\otimes n} = \frac{1}{\sqrt{2^{n}}} \sum\limits_{x \in B_n}{\left( 1-\frac{2 p_n x}{1- \frac{1}{2^n}}\right) \ket x}$  , entitled spue_simulation.ipynb,
* a Notebook giving depth and size plots of the MCX gate with first ancilla zeroed and second dirty, the Exact-one oracle and the block-encoding circuit of $1 - 2 p_n L_n$ , entitled depth_size_plots.ipynb.
