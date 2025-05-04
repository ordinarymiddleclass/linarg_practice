# Re-import necessary modules after code execution environment reset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
import os 

# Set random seed for reproducibility
np.random.seed(0)
# switch to the working directory to the directory of the script
os.chdir(os.path.dirname(os.path.abspath(__file__)))
# Re-define the symmetric matrix A and the vectors
A = np.array([[3, 1],
              [1, 2]])

theta = np.linspace(0, 2 * np.pi, 300)
vectors = np.array([np.cos(theta), np.sin(theta)])  # shape (2, N)

# Compute the normalized quadratic form
quad_form = np.einsum('ij,ji->i', vectors.T @ A, vectors)

# Re-create animation components
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlim(-4, 4)
ax.set_ylim(-4, 4)
ax.set_aspect('equal')
ax.grid(True)
ax.set_title("Dynamic Mapping Over Magnitude")

vec_line, = ax.plot([], [], 'r-', lw=2, label='Direction v')
mapped_line, = ax.plot([], [], 'b--', lw=2, label='Mapped length')
circle = plt.Circle((0, 0), 1, color='gray', fill=False, linestyle='--')
ax.add_patch(circle)
ax.legend()

def init():
    vec_line.set_data([], [])
    mapped_line.set_data([], [])
    return vec_line, mapped_line

def animate(i):
    v = vectors[:, i]
    z = quad_form[i]
    v_scaled = v * z
    vec_line.set_data([0, v[0]], [0, v[1]])
    mapped_line.set_data([0, v_scaled[0]], [0, v_scaled[1]])
    return vec_line, mapped_line

ani = animation.FuncAnimation(fig, animate, frames=len(theta),
                              init_func=init, blit=True, interval=30)

# Save as GIF
gif_path = "mapping_animation.gif"
ani.save(gif_path, writer=PillowWriter(fps=30))

# Display the GIF

