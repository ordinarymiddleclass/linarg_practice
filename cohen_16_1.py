"""
    cohen_16_1.py
    ~~~~~~~~~~~~~
    Code challenge from Cohen's linear algebra book, chapter 16, problem 1.
    This problem compares the full and economy SVDs of a 10x3 matrix A.
"""
import numpy as np

# Create reproducible 10x3 matrix
np.random.seed(0)
A = np.random.rand(10, 3)
print("Original Matrix A:")
print(np.round(A, 3))

# -----------------------
# Full SVD
# -----------------------
U_full, S_full, Vt_full = np.linalg.svd(A, full_matrices=True)

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
# Economy SVD
# -----------------------
U_econ, S_econ, Vt_econ = np.linalg.svd(A, full_matrices=False)
S_econ_matrix = np.diag(S_econ)

print("\n=== ECONOMY SVD ===")
print("U_econ (10x3):")
print(np.round(U_econ, 3))

print("\nS_econ (3x3):")
print(np.round(S_econ_matrix, 3))

print("\nVt_econ (3x3):")
print(np.round(Vt_econ, 3))

# Reconstruct matrix
A_econ_reconstructed = U_econ @ S_econ_matrix @ Vt_econ
print("\nReconstructed A (Economy SVD):")
print(np.round(A_econ_reconstructed, 3))

# -----------------------
# Compare reconstruction errors
# -----------------------
diff_full = np.linalg.norm(A - A_full_reconstructed)
diff_econ = np.linalg.norm(A - A_econ_reconstructed)

print(f"\nReconstruction error (Full SVD): {diff_full:.2e}")
print(f"Reconstruction error (Economy SVD): {diff_econ:.2e}")

# print the difference between economic U and first three columns of full U
print(U_full[:, :3] - U_econ)