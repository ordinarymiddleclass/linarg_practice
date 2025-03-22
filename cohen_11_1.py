"""
    cohen_11_1.py
    ~~~~~~~~~~~~~~~~
    Code challenge from Mike X Cohen's linear algebra book chapter 11, problem 1.
    This is a code that demonstrates that det(beta*M) = beta^n * det(M) for an n x n matrix M and a scalar beta.
"""
import numpy as np
import numpy.linalg as la
import random
import sys
sys.stdout = open("cohen_11_1_output.txt", "w")

# Generate a random 3x3 square matrix with integer values between 1 and 10
M = np.array([[random.randint(1,10) for _ in range(3)] for _ in range(3)])
# Print the matrix
print("Generated square matrix M: ")
print(M)

# Generate a random scalar beta
beta = random.randint(1,10)
# Print the scalar
print(f"Generated scalar beta: {beta}")

# Calculate the determinant of the matrix
det_M = la.det(M)
# Print the determinant of the matrix
print(f"Determinant of the square matrix M: {det_M}")

# Calculate the determinant of the beta*M matrix
det_beta_M = la.det(beta*M)
# Print the determinant of the beta*M matrix
print(f"Determinant of the beta*M matrix: {det_beta_M}")

# Calculate the determinant of the beta*M matrix using the property det(beta*M) = beta^n * det(M)
det_beta_M_property = beta**3 * det_M

print(f"Determinant of the beta*M matrix using the property det(beta*M) = beta^n * det(M): {det_beta_M_property}")

