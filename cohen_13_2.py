"""
    cohen_13_2.py
    ~~~~~~~~~~~~~~~~
    Code challenge from Mike X Cohen's linear algebra book chapter 13, problem 2.
    This is a code that implements the Gram-Schmidt orthogonalization process.
"""
import numpy as np
import numpy.linalg as la
import random
import sys

sys.stdout = open("cohen_13_2_output.txt", "w")

def generate_square_matrix(n):
    return np.array([[random.randint(1, 20) for _ in range(n)] for _ in range(n)])

def generate_rectangular_matrix(m, n):
    return np.array([[random.randint(1, 20) for _ in range(n)] for _ in range(m)])

def normalize(v):
    """
    Normalize a vector v.
    """
    return v / la.norm(v)

def row_separator(M):
    """
    Separate the rows of a matrix M into vectors.
    """
    return [M[i, :] for i in range(M.shape[0])]

def column_separator(M):
    """
    Separate the columns of a matrix M into vectors.
    """
    return [M[:, i] for i in range(M.shape[1])]

def gram_schmidt(M):
    """
    Perform the Gram-Schmidt orthogonalization process on a matrix M.
    """
    # Get the rows of the matrix as vectors
    columns = column_separator(M)    
    # Initialize the list of orthogonal vectors
    ortho_vectors = []    
    # Iterate over the columns of the matrix
    for i in range(len(columns)):
        # Get the current column vector and convert to float 
        v = columns[i].astype(float)        
        # Subtract the projection of v onto the previous orthogonal vectors
        for j in range(i):
            v -= np.dot(v, ortho_vectors[j]) * ortho_vectors[j]
        
        # Normalize the orthogonal vector
        v = normalize(v)
        
        # Add the orthogonal vector to the list
        ortho_vectors.append(v)       
    
    #combine the orthogonal vectors to form the Q matrix
    Q = np.array(ortho_vectors).T
    # Return the orthogonalized matrix as a numpy array
    return Q

def gram_schmidt_qr(M):
    """
    calculate Q and R using the gram-schmidt process
    """
    Q = gram_schmidt(M)
    # calculate R using the property that R = Q^T * M
    R = Q.T @ M
    return Q, R

if __name__ == "__main__":
    # Generate a random 3x3 square matrix with integer values between 1 and 20
    M = generate_square_matrix(3)
    # Print the matrix
    print("Generated square matrix M: ")
    print(M)

    # Perform the QR decomposition of the matrix using the Gram-Schmidt process
    Q_gs, R_gs = gram_schmidt_qr(M)

    # Print the Q and R matrices obtained from the Gram-Schmidt process
    print("\nQ matrix obtained from the Gram-Schmidt process: ")
    print(Q_gs)
    print("\nR matrix obtained from the Gram-Schmidt process: ")
    print(R_gs)

    # Print Q^T*Q
    print("\nQ^T*Q: ")
    print(Q_gs.T @ Q_gs)
    
    # Print Q*R 
    print("\nQ*R: ")
    print(Q_gs @ R_gs)

    # Perform the QR decomposition of the matrix using the numpy library
    Q_np, R_np = la.qr(M)

    # Print the Q and R matrices obtained from the numpy library
    print("\nQ matrix obtained from the numpy library: ")
    print(Q_np)
    print("\nR matrix obtained from the numpy library: ")
    print(R_np)

    # Print Q*R
    print("\nQ*R obtained from the numpy library: ")
    print(Q_np @ R_np)
