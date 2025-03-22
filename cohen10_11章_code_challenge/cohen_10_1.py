"""
    cohen_10_1.py
    ~~~~~~~~~~~~~~~~
    Code challenge from Mike X Cohen's linear algebra book chapter 10, problem 1.
"""
import numpy as np
import numpy.linalg as la
import random
import sys
sys.stdout = open("cohen_10_1_output.txt", "w")

# Generate a random vector of length 3 with integer values between 1 and 10
x = np.array([random.randint(1,10) for _ in range(3)])
# Print the vector
print("Generated vector [v1, v2, v3]: ")
print(f'v1: {x[0]}', f'v2: {x[1]}', f'v3: {x[2]}')

# Generate a random matrix of size 3x3 with integer values between 1 and 10
A = np.array([[random.randint(1,10) for _ in range(3)] for _ in range(3)])
# System of equations
# a11x1 + a12x2 + a13x3 = y1
# a21x1 + a22x2 + a23x3 = y2
# a31x1 + a32x2 + a33x3 = y3
# Print the system of equations
print("Generated matrix [a11, a12, a13; a21, a22, a23; a31, a32, a33]: ")
print(f'a11: {A[0,0]}', f'a12: {A[0,1]}', f'a13: {A[0,2]}')
print(f'a21: {A[1,0]}', f'a22: {A[1,1]}', f'a23: {A[1,2]}')
print(f'a31: {A[2,0]}', f'a32: {A[2,1]}', f'a33: {A[2,2]}')

print("System of equations: ")
print(f'{A[0,0]}*v1 + {A[0,1]}*v2 + {A[0,2]}*v3 = y1')
print(f'{A[1,0]}*v1 + {A[1,1]}*v2 + {A[1,2]}*v3 = y2')
print(f'{A[2,0]}*v1 + {A[2,1]}*v2 + {A[2,2]}*v3 = y3')

# Calculate the dot product of the vector and the matrix
y = np.dot(A,x)
# Print the result
print("Result of the matrix - vector multiplication: ")
print(f'y1: {y[0]}', f'y2: {y[1]}', f'y3: {y[2]}')

# Print the matrix - vector multiplication process as a system of equations
print("Matrix - vector multiplication as a system of equations: ")
print(f'{A[0,0]}*v1 + {A[0,1]}*v2 + {A[0,2]}*v3 = {y[0]}')
print(f'{A[1,0]}*v1 + {A[1,1]}*v2 + {A[1,2]}*v3 = {y[1]}')
print(f'{A[2,0]}*v1 + {A[2,1]}*v2 + {A[2,2]}*v3 = {y[2]}')