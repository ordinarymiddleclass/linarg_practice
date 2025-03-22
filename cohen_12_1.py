"""
    cohen_12_1.py
    ~~~~~~~~~~~~~~~~
    Code challenge from Mike X Cohen's linear algebra book chapter 12, problem 1.
    This is a code that implements the matrix inverse using MCA (matrix of cofactors and adjugate matrix) algorithm.
"""
import numpy as np
import numpy.linalg as la
import random
import sys
sys.stdout = open("cohen_12_1_output.txt", "w")

def cofactor_matrix(A):
    """
    Calculate the cofactor matrix of the matrix A.
    """
    n = A.shape[0]
    C = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            C[i,j] = (-1)**(i+j) * la.det(np.delete(np.delete(A,i,0),j,1))
    return C

def adjugate_matrix(A):
    """
    Calculate the adjugate matrix of the matrix A.
    """
    return cofactor_matrix(A).T

def inverse_matrix(A):
    """
    Calculate the inverse of the matrix A.
    """
    return adjugate_matrix(A) / la.det(A)

def calculate_determinant(A):
    """
    Calculate the determinant of the matrix A.
    """
    return la.det(A)

def generate_matrix(m):
    """
    Generate a random m x m matrix, with integer values between -10 and 10.
    """
    return np.random.randint(-10,10,(m,m))

if __name__ == "__main__":    
    m = 3
    A = generate_matrix(m)
    print("Matrix A:")
    print(A)
    print("\nInverse of A using MCA algorithm:")
    print(inverse_matrix(A))
    print("\nDeterminant of A:")
    print(calculate_determinant(A))
    print("\nInverse of A using numpy:")
    print(la.inv(A))