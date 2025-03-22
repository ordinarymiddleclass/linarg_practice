"""
    demonstrate.py
    ~~~~~~~~~~~~~~
    code for demonstrating lstsq in scipy.  
"""
import numpy as np
from scipy.linalg import lstsq

#create a matrix A and a vector b
A = np.array([[1, 2], [3, 4], [5, 6]])
b = np.array([[3], [7], [11]])

#solve the system of equations
x, residuals, rank, s = lstsq(A, b)

#pretty print the results
print("A = ")
print(A)
print("b = ")
print(b)
print("x = ")
print(x)
print("residuals = ")
print(residuals)
print("rank = ")
print(rank)
print("s = ")
print(s)
