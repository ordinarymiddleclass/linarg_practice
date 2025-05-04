"""
    cohen_16_2.py
    ~~~~~~~~~~~~~
    Code challenge from Cohen's linear algebra book, chapter 16, problem 2.
    This problem compares svd obtained from eigendeposition of A^T A and the SVD of A.
    The SVD of A is obtained using the numpy function scipy.linalg.svd.
"""
import numpy as np
from scipy.linalg import svd
from scipy.linalg import eig

# Create reproducible 10x3 matrix
np.random.seed(0)
A = np.random.rand(10, 3)
print("Original Matrix A:")
print(np.round(A, 3))
# -----------------------
# Full SVD
# -----------------------
U_full, S_full, Vt_full = svd(A, full_matrices=True)
# Build full S matrix (10x3) for multiplication
S_full_matrix = np.zeros((10, 3))
S_full_matrix[:3, :3] = np.diag(S_full)
print("\n=== FULL SVD ===")
print("U_full (10x10):")
print(np.round(U_full, 3))
print("\nS_full (10x3):")
print(np.round(S_full_matrix, 3))
print("\nVt_full (3x3):")
print(np.round(Vt_full, 3))
# Reconstruct matrix
A_full_reconstructed = U_full @ S_full_matrix @ Vt_full
print("\nReconstructed A (Full SVD):")
print(np.round(A_full_reconstructed, 3))
# -----------------------
# SVD from eigen decomposition of A^T A
# -----------------------
# Compute A^T A
AtA = A.T @ A
print("\nA^T A:")
print(np.round(AtA, 3))
# Compute eigenvalues and eigenvectors
eigenvalues, eigenvectors = eig(AtA)
print("\nEigenvalues of A^T A:")
print(np.round(eigenvalues, 3))
print("\nEigenvectors of A^T A:")
print(np.round(eigenvectors, 3))
# Sort eigenvalues and eigenvectors
sorted_indices = np.argsort(eigenvalues)[::-1]
eigenvalues_sorted = eigenvalues[sorted_indices]
eigenvectors_sorted = eigenvectors[:, sorted_indices]
print("\nSorted Eigenvalues of A^T A:")
print(np.round(eigenvalues_sorted, 3))
print("\nSorted Eigenvectors of A^T A:")
print(np.round(eigenvectors_sorted, 3))
# Compute singular values
singular_values = np.sqrt(eigenvalues_sorted)
print("\nSingular Values from Eigenvalues of A^T A:")
print(np.round(singular_values, 3))
# Compute U from eigenvectors
U_eigen = A @ eigenvectors_sorted
# Normalize U
U_eigen_norm = np.linalg.norm(U_eigen, axis=0)
U_eigen_normalized = U_eigen / U_eigen_norm
print("\nU from Eigenvectors of A^T A (normalized):")
print(np.round(U_eigen_normalized, 3))
# Compute S from singular values
S_eigen = np.zeros((10, 3))
S_eigen[:3, :3] = np.diag(singular_values)
print("\nS from Singular Values of A^T A:")
print(np.round(S_eigen, 3))
# Compute V from eigenvectors
V_eigen = eigenvectors_sorted
# Normalize V
V_eigen_norm = np.linalg.norm(V_eigen, axis=0)
V_eigen_normalized = V_eigen / V_eigen_norm
print("\nV from Eigenvectors of A^T A (normalized):")
print(np.round(V_eigen_normalized, 3))
# Reconstruct matrix
A_eigen_reconstructed = U_eigen_normalized @ S_eigen @ V_eigen_normalized.T
print("\nReconstructed A (Eigen SVD):")
print(np.round(A_eigen_reconstructed, 3))
# -----------------------
# Compare reconstruction errors
# -----------------------
diff_full = np.linalg.norm(A - A_full_reconstructed)
diff_eigen = np.linalg.norm(A - A_eigen_reconstructed)
print("\nReconstruction error (Full SVD):", diff_full)
print("Reconstruction error (Eigen SVD):", diff_eigen)
print("Difference between Full SVD and Eigen SVD reconstruction:", np.linalg.norm(A_full_reconstructed - A_eigen_reconstructed))
print("Difference between Full SVD and Eigen SVD reconstruction (normalized):", np.linalg.norm(A_full_reconstructed - A_eigen_reconstructed) / np.linalg.norm(A))
print("Difference between Full SVD and Eigen SVD reconstruction (normalized):", np.linalg.norm(A_full_reconstructed - A_eigen_reconstructed) / np.linalg.norm(A))