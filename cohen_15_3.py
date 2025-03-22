"""
    cohen_15_3.py
    ~~~~~~~~~~~~~
    Code challenge from Cohen's linear algebra book, chapter 15, problem 3.
    This problem demonstrates eigendecomposition of Hankel matrices.
"""

from scipy.linalg import eig
from scipy.linalg import hankel
from numpy import array
import numpy as np
from random import randint
from matplotlib import pyplot as plt

# Create random 50x50 Hankel matrix A
A = hankel([randint(1, 10) for i in range(50)], [randint(1, 10) for i in range(50)])

# Pretty print the matrix
print("50x50 Hankel Matrix A:")
print(A)

# Perform eigendecomposition on A
eigvals, eigvecs = eig(A)
# Sort the eigenvalues and eigenvectors by descending order of eigenvalues
idx = eigvals.argsort()[::-1]
eigvals = eigvals[idx]
eigvecs = eigvecs[:, idx]
# Pretty print the eigenvalues and eigenvectors
print("Eigenvalues by eigendecomposition:")
print(eigvals)
print("Eigenvectors by eigendecomposition:")
print(eigvecs)

#Produce a figure that shows Hankel matrix A and its eigenvectors
plt.figure()
plt.imshow(A)
plt.title("Hankel Matrix A")
plt.show()

plt.figure()
plt.imshow(eigvecs)
plt.title("Eigenvectors of Hankel Matrix A")
plt.show()

#Plot the first 10 eigenvectors
plt.figure()
plt.plot(eigvecs[:, 0:3])
plt.title("First 3 Eigenvectors of Hankel Matrix A")
plt.show()
