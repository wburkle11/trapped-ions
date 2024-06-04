## **Q Chemistry**

One of the most promising near term applications of quantum computing in the NISQ era is the simulation of quantum chemical dynamics in increasing complex chemical systems. This directory is for storing code related to quantum chemistry projects and experiments

### Toeplitz Decomposition

For quantum chemical systems, the Hamiltonian $H$ is the sum of the potential and kinetic energies: $H = V + K$. The potential energy has a diagonal matrix form which can be decomposed into a sum of appropriately tensored and weighted $Z$ and $I$ operators. Since $V$ is already diagonal, we will focus our attention on the kinetic energy term $K$ which has a symmetric Toeplitz matrix form

![snip](https://github.com/wburkle11/trapped-ions/assets/92954143/48c6dd1a-157c-42cc-b854-52935191a234)

Where we decompose this Toeplitz matrix $K$ into a sum of generalized Pauli operators. This is a hardware inspired (and much more natural) decomposition to choose instead of Clifford gates because quantum platforms such as trapped-ions and superconducting circuits---at the hardware level---drive quantum gates with generalized Pauli operator interactions, not Clifford gates.

![decomp](https://github.com/wburkle11/trapped-ions/assets/92954143/e175e2b2-1124-49f9-a20c-dab0a80d2762)


