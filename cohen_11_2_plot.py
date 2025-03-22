"""
    cohen_11_2.py
    ~~~~~~~~~~~~~~~~
    Code challenge from Mike X Cohen's linear algebra book chapter 11, problem 2.
    This is a code that demonstrates numerical instability of the determinant calculation, with plot from size 3 to size 30.
"""
import numpy as np
import numpy.linalg as la
import random
import matplotlib.pyplot as plt
from cohen_11_2 import generate_matrix, calculate_determinant
import sys
sys.stdout = open("cohen_11_2_plot_output.txt", "w")

def determinant_singular(m):
    # calculate the determinant of a singular matrix of size m x m 100 times and return the average
    det_sum = 0
    for _ in range(100):
        A = generate_matrix(m)
        det_sum += calculate_determinant(A)
    return det_sum/100

if __name__ == "__main__":
    x = list(range(3, 31))    
    y = [(np.log(np.abs(determinant_singular(m)))) for m in x]
    for m in range(3, 31):
        det = determinant_singular(m)
        print(f"Average determinant of singular matrix of size {m} x {m}: {det}")
    # Plot the average determinant of singular matrices of size 3 x 3 to 30 x 30
    plt.figure()
    plt.xlabel('Matrix size')
    plt.ylabel('Log average determinant')
    plt.title('Average determinant of singular matrices')    
    plt.plot(x, y)
    plt.show()