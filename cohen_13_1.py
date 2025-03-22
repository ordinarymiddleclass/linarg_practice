"""
    cohen_13_1.py
    ~~~~~~~~~~~~~~~~
    Code challenge from Mike X Cohen's linear algebra book chapter 13, problem 1.
    This is a code that demonstrates sizes of Q and R matrices in the QR decomposition of a matrix.
"""

import numpy as np
import numpy.linalg as la
import random
import sys

#sys.stdout = open("cohen_13_1_output.txt", "w")

def generate_square_matrix(n):
    return np.array([[random.randint(1, 20) for _ in range(n)] for _ in range(n)])

def generate_rectangular_matrix(m, n):
    return np.array([[random.randint(1, 20) for _ in range(n)] for _ in range(m)])

if __name__ == "__main__":
    # Generate a random 3x3 square matrix with integer values between 1 and 20
    M = generate_square_matrix(3)
    # Print the matrix
    print("Generated square matrix M: ")
    print(M)

    # Perform the QR decomposition of the matrix
    Q, R = la.qr(M)

    # Print the sizes of the Q and R matrices
    print(f"Size of the Q matrix for square matrix M: {Q.shape}")
    print(f"Size of the R matrix for square matrix M: {R.shape}")

    # Print the ranks of the Q and R matrices
    print(f"Rank of the Q matrix for square matrix M: {la.matrix_rank(Q)}")
    print(f"Rank of the R matrix for square matrix M: {la.matrix_rank(R)}")

    # Generate a random 3x2 wide matrix with integer values between 1 and 20
    M_wide = generate_rectangular_matrix(3, 2)
    # Print the rectangular matrix
    print("\nGenerated wide matrix M_wide: ")
    print(M_wide)

    # Perform the QR decomposition of the rectangular matrix
    Q_wide, R_wide = la.qr(M_wide)

    # Print the sizes of the Q and R matrices for the rectangular matrix
    print(f"Size of the Q matrix for the wide matrix: {Q_wide.shape}")
    print(f"Size of the R matrix for the wide matrix: {R_wide.shape}")

    # Print the ranks of the Q and R matrices for the rectangular matrix
    print(f"Rank of the Q matrix for the wide matrix: {la.matrix_rank(Q_wide)}")
    print(f"Rank of the R matrix for the wide matrix: {la.matrix_rank(R_wide)}")

    # Generate a random 2x3 tall matrix with integer values between 1 and 20
    M_tall = generate_rectangular_matrix(2, 3)
    # Print the rectangular matrix
    print("\nGenerated tall matrix M_tall: ")
    print(M_tall)

    # Perform the QR decomposition of the rectangular matrix
    Q_tall, R_tall = la.qr(M_tall)

    # Print the sizes of the Q and R matrices for the rectangular matrix
    print(f"Size of the Q matrix for the tall matrix: {Q_tall.shape}")
    print(f"Size of the R matrix for the tall matrix: {R_tall.shape}")

    # Print the ranks of the Q and R matrices for the rectangular matrix
    print(f"Rank of the Q matrix for the wide matrix: {la.matrix_rank(Q_wide)}")
    print(f"Rank of the R matrix for the wide matrix: {la.matrix_rank(R_wide)}")
