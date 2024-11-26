import numpy as np
import numba

#@numba.njit
def generate_unique_triplets(N, L):
    """
    Generate N unique triplets of integers (x, y, z) where 0 <= x, y, z <= L.

    Args:
        N (int): Number of unique triplets to generate.
        L (int): Maximum value for each element in the triplet (0 to L).

    Returns:
        numpy.ndarray: Array of shape (N, 3) containing N unique triplets.
    """
    # Total possible unique triplets
    total_triplets = (L + 1) ** 3
    
    # Check if N is larger than the number of possible unique triplets
    if N > total_triplets:
        raise ValueError("N is larger than the total number of unique triplets available.")
    
    # Generate N unique indices within the range of possible triplets
    selected_indices = np.random.choice(total_triplets, size=N, replace=False)
    
    # Convert the indices back to triplets (x, y, z)
    z = selected_indices % (L + 1)
    y = (selected_indices // (L + 1)) % (L + 1)
    x = (selected_indices // ((L + 1) * (L + 1))) % (L + 1)
    
    
    # Stack them into an N x 3 array
    return np.stack((x, y, z), axis=-1)