"""
    cohen_11_2.py
    ~~~~~~~~~~~~~~~~
    Code challenge from Mike X Cohen's linear algebra book chapter 11, problem 2.
    This is a code that demonstrates numerical instability of the determinant calculation.
"""
import numpy as np
import numpy.linalg as la
import random
import sys
sys.stdout = open("cohen_11_2_output.txt", "w")

def generate_matrix(n):
    """
    Generate a random n x n square matrix with integer values between 1 and 10, and ensure that it is reduced rank (singular).
    """
    # Generate a random n-1 x n square matrix with integer values between 1 and 10
    A = np.array([[random.randint(1,10) for _ in range(n)] for _ in range(n-1)])
    # Generate the last row of the matrix to make it singular
    A = np.vstack((A, A[0] + A[1]))
    return A

def calculate_determinant(A):
    """
    Calculate the determinant of the matrix A.
    """
    return la.det(A)

def determinant_singular(n):
    """
    Calculate the determinant of a singular matrix 
    """
    A = generate_matrix(n)
    det_A = calculate_determinant(A)
    return det_A

if __name__ == "__main__":    
    for m in range(3,31):
        for i in range(10):
            det_A = determinant_singular(m)
            print(f"Determinant of the {m}x{m} singular matrix A, repeat {i+1}: {det_A}")