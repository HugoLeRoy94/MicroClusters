import numpy as np
import Objects
from numpy.typing import NDArray
from typing import List as List, Tuple

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

def to_single_index(x: int, y: int, z: int, L: int) -> int:
    """Convert 3D indices (x, y, z) to a single index."""
    return x * L * L + y * L + z

def to_xyz(index: int, L: int) -> Tuple[int, int, int]:
    """Convert a single index back to 3D indices (x, y, z)."""
    z = index % L
    y = (index // L) % L
    x = index // (L * L)
    return x, y, z

class box:
    def __init__(self, 
                 size: int,
                 nobjects:int,
                 Interactions: List[List[float]]) -> None:
        self.size = size
        self.lattice = np.zeros((self.size, self.size, self.size), dtype=Objects.Object)
        for x in range(self.size):
            for y in range(self.size):
                for z in range(self.size):
                    self.lattice[x,y,z] = Objects.Empty((x,y,z))
        
        self.objects = np.empty((nobjects), dtype=Objects.Object)

        self.E = Interactions # symetric matrix of unspecific interactions

    def set_lattice(self,site:Tuple[int,int,int],object:Objects)->None:
        self.lattice[site]=object
        object.position = site

    def get_lattice(self,site:Tuple[int,int,int])->None:
        return self.lattice[site]
        

    def compute_local_energy(self, xyz:Tuple) -> float:
        """ Computes the local energy for a given lattice site. """
        if self.lattice[xyz].isempty():
            return 0.0
        neighbors = self.get_neighbors(xyz)
        local_energy = 0.0


        neigh_count = 0.0 # for later
        for nxyz in neighbors:            
            local_energy -= self.E[self.lattice[xyz].Index()][self.lattice[nxyz].Index()]

        return local_energy
    def total_energy(self) -> float:
        """ Computes the total energy of the system. 
            The function is really not efficient, but only used once at the beginning.
        """
        energy = 0.0
        for x in range(self.size):
            for y in range(self.size):
                for z in range(self.size):
                    if not self.lattice[x,y,z].isempty():
                        xyz=(x,y,z)
                        energy += self.compute_local_energy(xyz)
        return energy / 2  # To correct for double counting of pairs


    def get_neighbors(self, xyz:Tuple[int,int,int]) -> List[np.ndarray]:
        """ Returns the coordinates of the 26 nearest neighbors of a given lattice site (x, y, z). """
        neighbors = []
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                for dz in [-1, 0, 1]:
                    if dx == 0 and dy == 0 and dz == 0:
                        continue
                    neighbors.append(((xyz[0] + dx) % self.size, 
                                                (xyz[1] + dy) % self.size, 
                                                (xyz[2] + dz) % self.size))
        return neighbors

    def has_free_neighbor(self, xyz:Tuple[int,int,int])-> bool:
        """
        Check if a site has a surounding free site
        """
        neighbors = self.get_neighbors(xyz)
        for nxyz in neighbors:
            if self.lattice[nxyz].isempty():
                return True
        return False

        
    def build_clusters(self) -> Tuple[NDArray[np.int_], NDArray[np.int_]]:
        """
        Computes the distribution of cluster sizes in the lattice.
        A cluster is defined as a set of connected indices where lattice[index] = True.
        The connected condition is defined by the get_neighbors function.

        Returns:
            Tuple[NDArray[np.int_], NDArray[np.int_]]: A tuple containing an array of cluster indices and an array of cluster start indices.
        """
        visited = np.zeros(self.size ** 3, dtype=np.bool_)  # Flattened 3D lattice
        cluster_indices = np.empty(self.size ** 3, dtype=np.int32)  # Allocate maximum size, we'll trim later
        cluster_starts = np.empty(self.size ** 3, dtype=np.int32)  # Same here
        cluster_count = 0
        current_index = 0

        def flood_fill(start_index: int) -> int:
            """Performs flood fill from the start index to find the entire cluster."""
            nonlocal current_index,cluster_indices,visited,cluster_starts,cluster_count
            stack = list([start_index])
            cluster_start = current_index

            while len(stack) > 0:
                current_idx = stack.pop()
                if not visited[current_idx]:
                    visited[current_idx] = True
                    cluster_indices[current_index] = current_idx
                    current_index += 1
                    x, y, z = to_xyz(current_idx, self.size)
                    for neighbor in self.get_neighbors((x, y, z)):
                        nx, ny, nz = neighbor[0], neighbor[1], neighbor[2]
                        neighbor_idx = to_single_index(nx, ny, nz, self.size)
                        if self.lattice[nx, ny, nz] and not visited[neighbor_idx]:
                            stack.append(neighbor_idx)
            if current_index > cluster_start:
                cluster_starts[cluster_count] = cluster_start
                cluster_count += 1

        # Iterate over all lattice sites
        for x in range(self.size):
            for y in range(self.size):
                for z in range(self.size):
                    idx = to_single_index(x, y, z, self.size)
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
        for x in range(self.size):
            for y in range(self.size):
                for z in range(self.size):
                    if self.lattice[x,y,z]:
                        counter = 0
                        for neighbor in self.get_neighbors((x, y, z)):
                            nx, ny, nz = neighbor[0], neighbor[1], neighbor[2]
                            if self.lattice[nx,ny,nz]:
                                counter+=1
                        av_Nneigh+=counter/self.true_sites.shape[0]
        return av_Nneigh