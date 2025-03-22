"""
    cohen_13_3.py
    ~~~~~~~~~~~~~~~~
    Code challenge from Mike X Cohen's linear algebra book chapter 13, problem 3.
    This is a code that use Gram-Schmidt orthogonalization in cohen_13_2.py to compute the QR decomposition of a predefined wide matrix.

    The matrix is:
    M = [[1,1,-2],
         [3,-1,1]]
"""
import numpy as np
import numpy.linalg as la
# Import the gram_schmidt_qr function from cohen_13_2.py
from gscohen import gram_schmidt_qr
import sys

sys.stdout = open("cohen_13_3_output.txt", "w")

# Define the wide matrix M
M = np.array([[1, 1, -2], [3, -1, 1]])
# Print the matrix
print("Wide matrix M: ")
print(M)

# Perform the QR decomposition of the matrix using the Gram-Schmidt orthogonalization
Q, R = gram_schmidt_qr(M)

# Print the Q and R matrices
print("\nQ matrix:")
print(Q)
print("\nR matrix:")
print(R)

# Check the product of Q and R to see if it equals the original matrix M
print("\nQ * R:")
print(np.dot(Q, R))
print("\nOriginal matrix M:")
print(M)
print("\nAre Q * R and M equal?")
print(np.allclose(np.dot(Q, R), M))