"""
    cohen_15_1.py
    ~~~~~~~~~~~~~
    Code challenge from Cohen's linear algebra book, chapter 15, problem 1.
    This problem demonstrates generalized eigendecomposition.
"""
from scipy.linalg import eig
from numpy import array
import numpy as np
from random import randint

# Create random 2x2 matrix A
A = array([[randint(1, 10) for i in range(2)] for j in range(2)])
# Create random 2x2 matrix B
B = array([[randint(1, 10) for i in range(2)] for j in range(2)])

#Pretty print the matrices
print("2x2 Matrix A:")
print(A)
print("2x2 Matrix B:")
print(B)

# Perform generalized eigendecomposition on A and B
eigvals, eigvecs = eig(A, B)
# Pretty print the eigenvalues and eigenvectors
print("Eigenvalues by generalized eigendecomposition:")
print(eigvals)
print("Eigenvectors by generalized eigendecomposition:")
print(eigvecs)

# Perform eigendecomposition on inverse of B times A
eigvals, eigvecs = eig(np.linalg.inv(B).dot(A))
# Pretty print the eigenvalues and eigenvectors
print("Eigenvalues by eigendecomposition of inverse of B times A:")
print(eigvals)
print("Eigenvectors by eigendecomposition of inverse of B times A:")
print(eigvecs)

#Perform the procedure above for a 10x10 matrix
A = array([[randint(1, 10) for i in range(10)] for j in range(10)])
B = array([[randint(1, 10) for i in range(10)] for j in range(10)])

print("10x10 Matrix A:")
print(A)
print("10x10 Matrix B:")
print(B)

eigvals, eigvecs = eig(A, B)
print("Eigenvalues by generalized eigendecomposition:")
print(eigvals)
#print("Eigenvectors by generalized eigendecomposition:")
#print(eigvecs)

eigvals, eigvecs = eig(np.linalg.inv(B).dot(A))
print("Eigenvalues by eigendecomposition of inverse of B times A:")
print(eigvals)
#print("Eigenvectors by eigendecomposition of inverse of B times A:")
#print(eigvecs)

#Perform the procedure above for a 50x50 matrix
A = array([[randint(1, 10) for i in range(50)] for j in range(50)])
B = array([[randint(1, 10) for i in range(50)] for j in range(50)])

print("50x50 Matrix A:")
print(A)
print("50x50 Matrix B:")
print(B)

eigvals, eigvecs = eig(A, B)
print("Eigenvalues by generalized eigendecomposition:")
print(eigvals)
#print("Eigenvectors by generalized eigendecomposition:")
#print(eigvecs)

eigvals, eigvecs = eig(np.linalg.inv(B).dot(A))
print("Eigenvalues by eigendecomposition of inverse of B times A:")
print(eigvals.sort())
#print("Eigenvectors by eigendecomposition of inverse of B times A:")
#print(eigvecs)



