# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 10:50:52 2024

@author: iontrap
"""

#Computing Eigenvectors and Eigenvalues from a given hamiltonian: 
import numpy as np
import sympy as sp 

#Example

H = np.array([[1, 0, 0, 0],
              [0, 2, 0, 0],
              [0, 0, 3, 0],
              [0, 0, 0, 4]])

#eigenvalues, eigenvectors = np.linalg.eigh(H)

#print(eigenvalues)
#print(eigenvectors)

#Now we are going to try with our real hamiltonian

#Define Rabi Freq. 
omega_sigm = 2*(np.pi) * 17e6
omega_sigp = omega_sigm
omega_pi = 2*(np.pi) * 4e6

#Define Detunings
delta = 2*(np.pi) * 50e6
zeeman = 2*(np.pi) * 5e6

H1 = np.array([[0, omega_sigm/2, -omega_pi/2, omega_sigp/2],
               [omega_sigm/2, delta, 0, 0],
               [-omega_pi/2, 0, delta, 0],
               [omega_sigp/2, 0, 0, delta - 2*zeeman]])

eigenvalues1, eigenvectors1 = np.linalg.eigh(H1)

#print(eigenvalues1)
#print(eigenvectors1)


#Lets see if we can do this symbolically: 
    
#define symbolic variables

a, b, c, d, e = sp.symbols('a b c d e')

#define hamiltonian matrix

H2 = sp.Matrix([[0, a/2, -b/2, c/2],
               [a/2, d, 0, 0],
               [-b/2, 0, d, 0],
               [c/2, 0, 0, d-2*e]])

eigenvalues2 = H2.eigenvals()

eigenvectors2 = H2.eigenvects()

#print(eigenvalues2)

#DID NOT WORK 