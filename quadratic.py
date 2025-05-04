import numpy as np 
import matplotlib.pyplot as plt
import random 

def compute_Z(X, Y, A, eps=1e-2):
    denom = X**2 + Y**2
    mask = denom > eps  # avoid division by near-zero
    Z = np.zeros_like(X)
    
    # Compute where denominator is safe
    Z[mask] = (A[0, 0]*X[mask]**2 + A[1, 1]*Y[mask]**2 + 2*A[0, 1]*X[mask]*Y[mask]) / denom[mask]
    
    # Handle origin (0,0) if needed — define a value or keep it 0
    Z[~mask] = np.nan  # or set to some known limit if direction-independent
    
    return Z

# Create a random 2x2 matrix
random.seed(0)
A = np.random.rand(2, 2)
print("Original Matrix A:")
print(A)
"""
# Plot the quadratic form of matrix A
x = np.linspace(-2, 2, 100)
y = np.linspace(-2, 2, 100)
X, Y = np.meshgrid(x, y)
Z = A[0, 0]*X**2 + A[1, 1]*Y**2 + 2*A[0, 1]*X*Y
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
plt.title('Quadratic Form of Matrix A')
plt.show()
"""
# Plot the normalized quadratic form
x = np.linspace(-2, 2, 100)
y = np.linspace(-2, 2, 100)
X, Y = np.meshgrid(x, y)
Z = compute_Z(X, Y, A)
# Plot the normalized quadratic form
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
plt.title('Normalized Quadratic Form of Matrix A')
plt.show()
