## Q Chemistry 

One of the most promising near term applications of quantum computing in the NISQ era is the simulation of quantum chemical dynamics in increasing complex chemical systems. This directory is for storing code related to quantum chemistry projects and experiments

### Toeplitz Decomposition

How might the state dynamics of a Hamiltonian $H$ with a Toeplitz form be efficiently simulated with trapped ions? In part we are answering this question because we know that with efficiently simulated state dynamics we can reconstruct the transition energies between eigenstates in quantum chemical systems. It might also be used in a quantum algorithm, such as QPE, that measures exact eigenvalues of $H$. It would also be useful in real-space dynamics simulations or even something akin to a finite difference algorithm. We are also interested in answering this question because it shows there are some problems where global MS gates are more desirable than pair-wise MS gates. The main advantage being that global MS gates can generate multi-qubit gates with fewer laser drives and fewer calibrations likely leading to higher fidelity and a less engineering intensive experimental setup.

For the quantum chemical systems, the Hamiltonian $H$ is the sum of the potential and kinetic energies: $H = V + K$. The potential energy has a diagonal matrix form which can be decomposed into a sum of appropriately tensored and weighted $Z$ and $I$ operators. Since $V$ is already diagonal, we will focus our attention on the kinetic energy term $K$ which has a symmetric Toeplitz matrix form

![snip](https://github.com/wburkle11/trapped-ions/assets/92954143/48c6dd1a-157c-42cc-b854-52935191a234)

Where we decompose this Toeplitz matrix $K$ into a sum of generalized Pauli operators. This is a hardware inspired (and much more natural) decomposition to choose instead of Clifford gates because quantum platforms such as trapped-ions and superconducting circuits---at the hardware level---drive quantum gates with generalized Pauli operator interactions, not Clifford gates.

![decomp](https://github.com/wburkle11/trapped-ions/assets/92954143/e175e2b2-1124-49f9-a20c-dab0a80d2762)


