# -*- coding: utf-8 -*-
"""
This is a Decomposition of 4-Qubit unitary toeplitz matrices into a pauli basis
"""

# Importing standard Qiskit libraries
from qiskit import QuantumCircuit, transpile, Aer, IBMQ
from qiskit.tools.jupyter import *
from qiskit.visualization import *
from qiskit.providers.aer import QasmSimulator



#Function used for plotting Toeplitz sideband matrices

def plot_matrix(m):
    fig = plt.figure(figsize = (12, 14))
    ax1 = plt.subplot(1, 2, 1)
    ax2 = plt.subplot(1, 2, 2)
    im1 = ax1.imshow(np.imag(m))
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im1, cax=cax)
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    im2 = ax2.imshow(np.real(m))
    fig.colorbar(im2, cax=cax)
    ax1.set_xticks(np.arange(1, 17, dtype=np.int))
    ax2.set_xticks(np.arange(1, 17, dtype=np.int))
    ax1.set_yticks(np.arange(1, 17, dtype=np.int))
    ax2.set_yticks(np.arange(1, 17, dtype=np.int))
    plt.show()
    
from qiskit.opflow import X, Y, Z, I, MatrixOp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

xxxx = (X^X^X^X)

ixxx = (I^X^X^X)
xixx = (X^I^X^X)
xxix = (X^X^I^X)
xxxi = (X^X^X^I)

iixx = (I^I^X^X)
ixix = (I^X^I^X)
xiix = (X^I^I^X)
xxii = (X^X^I^I)
ixxi = (I^X^X^I)
xixi = (X^I^X^I)

iiix = (I^I^I^X)
xiii = (X^I^I^I)
ixii = (I^X^I^I)
iixi = (I^I^X^I)

yyyy = (Y^Y^Y^Y)

#iyyy = (I^Y^Y^Y)
#yiyy = (Y^I^Y^Y)
#yyiy = (Y^Y^I^Y) Any odd number of 'y' lead to an imaginary component. 
#yyyi = (Y^Y^Y^I)

iiyy = (I^I^Y^Y)
iyiy = (I^Y^I^Y)
yiiy = (Y^I^I^Y)
yyii = (Y^Y^I^I)
iyyi = (I^Y^Y^I)
yiyi = (Y^I^Y^I)

#iiiy = (I^I^I^Y)
#yiii = (Y^I^I^I)

ixxx = (I^X^X^X)
ixyy = (I^X^Y^Y)
iyyx = (I^Y^Y^X)
iyxy = (I^Y^X^Y)
xiyy = (X^I^Y^Y)

xyyi = (X^Y^Y^I)
yyxi = (Y^Y^X^I)
yxyi = (Y^X^Y^I)

xyiy = (X^Y^I^Y)
yyix = (Y^Y^I^X)
yiyx = (Y^I^Y^X)
yixy = (Y^I^X^Y)
yxiy = (Y^X^I^Y)

xyyx = (X^Y^Y^X)
yxyx = (Y^X^Y^X)
yyxx = (Y^Y^X^X)
xyxy = (X^Y^X^Y)
xxyy = (X^X^Y^Y)
yxxy = (Y^X^X^Y)


n1_est = iiix + 1/2 * (iixx + iiyy) + 1/4 * (ixxx -ixyy + iyyx + iyxy) + 1/8 * (xxxx - yyyy - xyyx + yxyx + yyxx - xyxy - xxyy + yxxy)
plot_matrix(n1_est.to_matrix())

n2_est = iixi + 1/2*(ixxi + iyyi) + 1/4 * (yxyi + xxxi - xyyi + yyxi)
plot_matrix(n2_est.to_matrix())

n3_est = 1/2 * (iixx - iiyy) + 1/2 * (ixix + iyiy)+ 1/4 * (xxix - xyiy + yyix + yxiy) + 1/4 * (ixxx +ixyy + iyyx - iyxy) + 1/8 * (xxxx + yyyy - xyyx + yxyx + yyxx + xyxy + xxyy - yxxy)
plot_matrix(n3_est.to_matrix())

n4_est = ixii + 1/2 * (xxii + yyii)
plot_matrix(n4_est.to_matrix())

n5_est = 1/2*(ixix - iyiy) + 1/4*(xxix + yyix + xixx + ixxx + yiyx - iyyx -yxiy + xyiy + yixy + iyxy - xiyy + ixyy) +1/8 * (xxxx + yyxx -yxyx + xyyx + yxxy -xyxy + xxyy + yyyy) 
plot_matrix(n5_est.to_matrix())

#n6_est = 1/2 * (ixxi - iyyi) + 1/2 * (xixi + yiyi) + 1/4 * (xxxi + (xyyi + yyxi - yxyi) 
#plot_matrix(n6_est.to_matrix())

n7_est = 1/2* (xiix + yiiy) + 1/4*(xixx + ixxx + yiyx - iyyx - yixy - iyxy + xiyy - ixyy) + 1/8*(xxxx + yyxx - yxyx + xyyx - yxxy + xyxy - xxyy - yyyy) 
plot_matrix(n7_est.to_matrix())

n8_est = xiii
plot_matrix(n8_est.to_matrix())

n9_est = 1/2* (xiix - yiiy) + 1/4*(xixx - yiyx + yixy + xiyy) + 1/8*(xxxx - yyxx + yxyx + xyyx + yxxy + xyxy - xxyy + yyyy)
plot_matrix(n9_est.to_matrix())

n10_est = 1/2 * (xixi - yiyi) + 1/4 * (xxxi + (xyyi - yyxi + yxyi))
plot_matrix(n10_est.to_matrix())

n11_est = 1/4*(xxix - yyix + xixx - yiyx + yxiy + xyiy - yixy - xiyy) + 1/8*(xxxx - yyxx + yxyx + xyyx - yxxy -xyxy + xxyy -yyyy)
plot_matrix(n11_est.to_matrix())

n12_est = xxii - yyii 
plot_matrix(n12_est.to_matrix())

n13_est = 1/4*(xxix - yyix - yxiy - xyiy) + 1/8*(xxxx - yyxx - yxyx - xyyx + yxxy + xyxy + xxyy -yyyy)
plot_matrix(n13_est.to_matrix())

n14_est = 1/4 * (xxxi - (xyyi + yyxi + yxyi))
plot_matrix(n14_est.to_matrix())

n15_est = 1/8 * (xxxx + yyyy - xyyx - yxyx - yyxx - xyxy - xxyy - yxxy)
plot_matrix(n15_est.to_matrix())