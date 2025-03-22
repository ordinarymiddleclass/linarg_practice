"""
    cohen_10_2.py
    ~~~~~~~~~~~~~~~~
    Code challenge from Mike X Cohen's linear algebra book chapter 10, problem 2.
"""
from sympy import Matrix
import random
#import a module to screen-print to a file
import sys
sys.stdout = open("cohen_10_2_output.txt", "w")

# Generate a random 3x3 square matrix with integer values between 1 and 10
A = Matrix([[random.randint(1,10) for _ in range(3)] for _ in range(3)])
# Print the matrix
print("Generated square matrix A: ")
print(A)

# Reduce the matrix to reduced row echelon form
rref, pivots = A.rref()
# Print the reduced row echelon form of the matrix
print("Reduced row echelon form of the square matrix A: ")
print(rref)

# Generate a wide matrix with 3 rows and 4 columns
B = Matrix([[random.randint(1,10) for _ in range(4)] for _ in range(3)])
# Print the matrix
print("Generated wide matrix B: ")
print(B)

# Reduce the wide matrix to reduced row echelon form
rref, pivots = B.rref()
# Print the reduced row echelon form of the matrix
print("Reduced row echelon form of the wide matrix B: ")
print(rref)

# Generate a tall matrix with 4 rows and 3 columns
C = Matrix([[random.randint(1,10) for _ in range(3)] for _ in range(4)])
# Print the matrix
print("Generated tall matrix C: ")
print(C)

# Reduce the tall matrix to reduced row echelon form
rref, pivots = C.rref()
# Print the reduced row echelon form of the matrix
print("Reduced row echelon form of the tall matrix C: ")
print(rref)
