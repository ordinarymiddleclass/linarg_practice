"""
    cohen_15_2.py
    ~~~~~~~~~~~~~
    Code challenge from Cohen's linear algebra book, chapter 15, problem 2.
    This problem demonstrates eigendecomposition of diagonal matrices.
"""

from scipy.linalg import eig
from numpy import array
import numpy as np
from random import randint

# Create random 5x5 diagonal matrix A
A = array([[randint(1, 10) if i == j else 0 for i in range(5)] for j in range(5)])

# Pretty print the matrix
print("5x5 Diagonal Matrix A:")
print(A)

# Perform eigendecomposition on A
eigvals, eigvecs = eig(A)
# Pretty print the eigenvalues and eigenvectors
print("Eigenvalues by eigendecomposition:")
print(eigvals)
print("Eigenvectors by eigendecomposition:")
print(eigvecs)
