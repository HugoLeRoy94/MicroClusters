import numpy as np
import numba
from numba.experimental import jitclass
from numba import int32, float64, boolean
from numba.typed import List  # Import Numba's List
from typing import List as PyList, Tuple
from numpy.typing import NDArray

@numba.njit
def generate_unique_triplets(N: int, L: int) -> np.ndarray:
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

@numba.njit
def to_single_index(x: int, y: int, z: int, L: int) -> int:
    """Convert 3D indices (x, y, z) to a single index."""
    return x * L * L + y * L + z

@numba.njit
def to_xyz(index: int, L: int) -> Tuple[int, int, int]:
    """Convert a single index back to 3D indices (x, y, z)."""
    z = index % L
    y = (index // L) % L
    x = index // (L * L)
    return x, y, z

spec = [
    ('lattice', boolean[:,:,:]),
    ('L', int32),
    ('E0', float64),
    ('E1', float64),
    ('T', float64),
    ('true_sites', int32[:,:])
]


@jitclass(spec)
class BOX:
    def __init__(self, size: int, nparticles: int, npolymers: int, lpolymer: int, Interactions: List[float], temperature: float) -> None:
        self.lattice = np.zeros((size, size, size), dtype=np.bool_)

        self.generate_polymers(npolymers,lpolymer)

        particles = generate_unique_triplets(nparticles, size - 1)
        self.L = size
        self.true_sites = np.empty((nparticles, 3), dtype=np.int32)
        for idx,particle in enumerate(particles):
            self.lattice[particle[0], particle[1], particle[2]] = True
            self.true_sites[idx] = particle
        self.E0 = Interactions[0]
        self.E1 = Interactions[1]
        self.T = temperature
 
    
    def add_random_poly(self,lpolymer:int)-> bool:
        # place the first seed in a random lattice
        return
             #while True:
             #    x,y,z = np.random.randint(self.lattice.shape[0]),np.random.randint(self.lattice.shape[1]),np.random.randint(self.lattice.shape[2])
             #    if self.lattice[x,y,z]==0 and self.has_free_neighbor(x,y,z):
             #        break
             #lpoly = 0
             #while lpoly < lpolymer:
             #    while True:
             #        nxyz = self.get_neighbors(x,y,z)
             #        nx,ny,nz = nxyz[0],nxyz[1],nxyz[2]
             #        if self.lattice[nx,ny,nz]==0 and  
 #
             #    x,y,z=nx,ny,nz
             #    lpoly+=1

    def generate_polymers(self, npolymers: int,lpolymer:int)->None:
        """
        generates n polymers made of l monomers.
        Each polymer is a self avoiding random walk.
        """
        npoly = 0
        while npoly < npolymers:
            if self.add_random_poly(lpolymer):
                npoly+=1


        return

    def has_free_neighbor(self, x:int,y:int,z:int)-> bool:
        """
        Check if a site has a surounding free site
        """
        neighbors = self.get_neighbors(x, y, z)
        for nxyz in neighbors:
            nx, ny, nz = nxyz[0], nxyz[1], nxyz[2]
            if self.lattice[nx,ny,nz]==0:
                return True
        return False
    
    def get_neighbors(self, x: int, y: int, z: int) -> List[np.ndarray]:
        """ Returns the coordinates of the 26 nearest neighbors of a given lattice site (x, y, z). """
        neighbors = []
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                for dz in [-1, 0, 1]:
                    if dx == 0 and dy == 0 and dz == 0:
                        continue
                    neighbors.append(np.array([(x + dx) % self.L, 
                                               (y + dy) % self.L, 
                                               (z + dz) % self.L]))
        return neighbors
    
    def compute_local_energy(self, x: int, y: int, z: int) -> float:
        """ Computes the local energy for a given lattice site. """
        if not self.lattice[x, y, z]:
            return 0.0
        neighbors = self.get_neighbors(x, y, z)
        local_energy = 0.0
        
        neigh_count = 0.0
        for nxyz in neighbors:
            nx, ny, nz = nxyz[0], nxyz[1], nxyz[2]
            if self.lattice[nx, ny, nz]:  # Interaction with other particles
                local_energy -= self.E0
                if neigh_count < 1:
                    local_energy -= self.E1
                    neigh_count += 1

        return local_energy
    
    def total_energy(self) -> float:
        """ Computes the total energy of the system. 
            The function is really not efficient, but only used once at the beginning.
        """
        energy = 0.0
        for x in range(self.L):
            for y in range(self.L):
                for z in range(self.L):
                    if self.lattice[x, y, z]:
                        energy += self.compute_local_energy(x, y, z)
        return energy / 2  # To correct for double counting of pairs
    
    def monte_carlo_step(self) -> bool:
        """Performs a single Monte Carlo step using the Metropolis algorithm with site exchange."""
        # Select a random True site
        idx = np.random.randint(0, len(self.true_sites))
        x, y, z = self.true_sites[idx]

        # Select a random site to exchange with
        x2 = np.random.randint(0, self.L)
        y2 = np.random.randint(0, self.L)
        z2 = np.random.randint(0, self.L)
        
        # Compute the energy before the move
        initial_energy = self.compute_local_energy(x, y, z) + self.compute_local_energy(x2, y2, z2)
        
        # Swap the particles (True <-> False)
        self.lattice[x, y, z], self.lattice[x2, y2, z2] = self.lattice[x2, y2, z2], self.lattice[x, y, z]
        
        # Compute the energy after the move
        final_energy = self.compute_local_energy(x, y, z) + self.compute_local_energy(x2, y2, z2)
        
        # Calculate the energy difference
        delta_e = final_energy - initial_energy
        
        # Decide whether to accept the move
        if delta_e > 0 and np.random.rand() >= np.exp(-delta_e / self.T):
            # Reject the move (revert)
            self.lattice[x, y, z], self.lattice[x2, y2, z2] = self.lattice[x2, y2, z2], self.lattice[x, y, z]
            return False  # Move was rejected
        # Move was accepted, update the list of True sites
        if self.lattice[x, y, z]:  # (x, y, z) is now True
            self.true_sites[idx] = np.array([x, y, z], dtype=np.int32)
        else:  # (x2, y2, z2) is now True
            self.true_sites[idx] = np.array([x2, y2, z2], dtype=np.int32)

        return True  # Move was accepted
    def monte_carlo_steps(self, steps: int) -> NDArray[np.bool_]:
        success = np.zeros(steps,dtype = np.bool_)
        for step in range(steps):
            success[step] = self.monte_carlo_step()
        return success
    def build_clusters(self) -> Tuple[NDArray[np.int_], NDArray[np.int_]]:
        """
        Computes the distribution of cluster sizes in the lattice.
        A cluster is defined as a set of connected indices where lattice[index] = True.
        The connected condition is defined by the get_neighbors function.

        Returns:
            Tuple[NDArray[np.int_], NDArray[np.int_]]: A tuple containing an array of cluster indices and an array of cluster start indices.
        """
        visited = np.zeros(self.L ** 3, dtype=np.bool_)  # Flattened 3D lattice
        cluster_indices = np.empty(self.L ** 3, dtype=np.int32)  # Allocate maximum size, we'll trim later
        cluster_starts = np.empty(self.L ** 3, dtype=np.int32)  # Same here
        cluster_count = 0
        current_index = 0

        def flood_fill(start_index: int) -> int:
            """Performs flood fill from the start index to find the entire cluster."""
            nonlocal current_index,cluster_indices,visited,cluster_starts,cluster_count
            stack = List([start_index])
            cluster_start = current_index

            while len(stack) > 0:
                current_idx = stack.pop()
                if not visited[current_idx]:
                    visited[current_idx] = True
                    cluster_indices[current_index] = current_idx
                    current_index += 1
                    x, y, z = to_xyz(current_idx, self.L)
                    for neighbor in self.get_neighbors(x, y, z):
                        nx, ny, nz = neighbor[0], neighbor[1], neighbor[2]
                        neighbor_idx = to_single_index(nx, ny, nz, self.L)
                        if self.lattice[nx, ny, nz] and not visited[neighbor_idx]:
                            stack.append(neighbor_idx)
            if current_index > cluster_start:
                cluster_starts[cluster_count] = cluster_start
                cluster_count += 1

        # Iterate over all lattice sites
        for x in range(self.L):
            for y in range(self.L):
                for z in range(self.L):
                    idx = to_single_index(x, y, z, self.L)
                    if self.lattice[x, y, z] and not visited[idx]:
                        # Start a new cluster
                        flood_fill(idx)

        return cluster_indices[:current_index], cluster_starts[:cluster_count]
            

    def cluster_size(self) -> NDArray[np.int32]:
        cluster_indices, cluster_starts = self.build_clusters()
        cluster_starts = np.append(cluster_starts, len(cluster_indices))
        sizes = np.diff(cluster_starts)
        return sizes
    def compute_av_Nneigh(self) -> int:
        av_Nneigh= 0
        for x in range(self.L):
            for y in range(self.L):
                for z in range(self.L):
                    if self.lattice[x,y,z]:
                        counter = 0
                        for neighbor in self.get_neighbors(x, y, z):
                            nx, ny, nz = neighbor[0], neighbor[1], neighbor[2]
                            if self.lattice[nx,ny,nz]:
                                counter+=1
                        av_Nneigh+=counter/self.true_sites.shape[0]
        return av_Nneigh