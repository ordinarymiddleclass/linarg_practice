"""
    cohen_16_5.py
    ~~~~~~~~~~~~~
    Code challenge from Cohen's linear algebra book, chapter 16, problem 5.
    This problem download a picture from wikimedia common, convert it to grayscale,
    and perform SVD on the grayscale image.
    The SVD is performed using the numpy function scipy.linalg.svd.
    The image is then reconstructed using the first k singular values and vectors.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import svd
from PIL import Image
import requests
from io import BytesIO
import os

def download_image(url):
    """
    Download an image from a URL and return it as a PIL Image object.
    """
    response = requests.get(url)
    img = Image.open(BytesIO(response.content))
    return img

def image_k (image, k):
    """
    Reconstruct the image using the first k singular values and vectors.
    """
    # Convert the image to grayscale
    img_gray = image.convert("L")
    
    # Convert the grayscale image to a numpy array
    A = np.array(img_gray)
    
    # Perform SVD
    U, S, Vt = svd(A, full_matrices=False)
    
    # Reconstruct the image using the first k singular values and vectors
    S_k = np.zeros((k, k))
    np.fill_diagonal(S_k, S[:k])
    
    A_reconstructed = U[:, :k] @ S_k @ Vt[:k, :]
    m  = U.shape[0]
    k = S_k.shape[0]
    n = Vt.shape[1]
    A_size = k*(m + n + k)

    return A_reconstructed, A_size

if __name__ == "__main__":
    #switch to the working directory to the directory of the script
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    #open "mercedes.jpg" from the working directory
    img = Image.open("mercede.jpg")    
    # Convert the image to grayscale
    img = img.convert("L")

    # Convert the grayscale image to a numpy array
    A = np.array(img)
    # print the array 
    #print(A)
    m = A.shape[0]
    n = A.shape[1]
    size_original = m * n
    print(f"Original image size: {size_original} pixels")

    # Perform SVD with k=10, k=20, k=50, k=100, k=500
    ks = [10, 20, 50, 100, 500]
    images = []
    size_reconstructed = []

    for k in ks:
        A_reconstructed, A_size = image_k(img, k)
        images.append(A_reconstructed)
        size_reconstructed.append(A_size/size_original)
    
    # Plot the original and reconstructed images on a 2x3 grid
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes[0, 0].imshow(A, cmap='gray')
    axes[0, 0].set_title("Original Image")
    axes[0, 0].axis('off')
    for i, k in enumerate(ks):
        j = i + 1 
        x = j // 3
        y = j % 3
        axes[x, y].imshow(images[i], cmap='gray')
        axes[x, y].set_title(f"Reconstructed Image (k={k})")
        #display the size of the reconstructed image at the bottom of the image
        axes[x, y].text(0, -10, f"Size: {size_reconstructed[i]:.2f}x", fontsize=12, ha='center', va='top', color='white')
        axes[x, y].set_xticks([])
        axes[x, y].set_yticks([])
        axes[x, y].axis('off')
    plt.tight_layout()
    plt.savefig("mercedes_reconstructed.png", dpi=300, bbox_inches='tight')
    plt.show()        
    plt.close()
    
    