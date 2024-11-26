import numpy as np
import BOX
import Objects
from numpy.typing import NDArray
from typing import List as List, Tuple

def to_single_index(x: int, y: int, z: int, L: int) -> int:
    """Convert 3D indices (x, y, z) to a single index."""
    return x * L * L + y * L + z

def to_xyz(index: int, L: int) -> Tuple[int, int, int]:
    """Convert a single index back to 3D indices (x, y, z)."""
    z = index % L
    y = (index // L) % L
    x = index // (L * L)
    return x, y, z
class MC:
    def __init__(self, 
                 size: int, 
                 nparticles: int, 
                 npolymers: int, 
                 lpolymer: int, 
                 interactions: List[List[float]], 
                 temperature: float) -> None:
        self.npolymers=npolymers
        self.lpolymer=lpolymer
        self.nparticles=nparticles
        self.E = interactions
        self.T = temperature
        
        self.box = BOX.box(size,nparticles+npolymers,interactions)

        self.generate_polymers(npolymers,lpolymer)
        self.generate_particles(nparticles)

    def generate_polymers(self,npolymers:int,lpolymer:int)->None:            
        """
        generates n polymers made of l monomers.
        Each polymer is a self avoiding random walk.
        """
        npoly = 0
        while npoly < npolymers:
            if self.add_random_poly(lpolymer):
                npoly+=1
        return True
    def add_random_poly(self,lpolymer:int):
        return
    
    def generate_particles(self,nparticles:int)->None:
        particles = BOX.generate_unique_triplets(nparticles, self.box.size - 1)
        for idx,particle in enumerate(particles):
            newdhh1 = Objects.DHH1((particle[0], particle[1], particle[2]))
            self.box.set_lattice((particle[0],particle[1],particle[2]),newdhh1)
            self.box.objects[self.npolymers+idx]=newdhh1 # polymer must be generated first

    def add_random_poly(self,lpolymer:int)-> bool:
        # place the first seed in a random lattice
            while True:
                x,y,z = np.random.randint(self.lattice.shape[0]),np.random.randint(self.lattice.shape[1]),np.random.randint(self.lattice.shape[2])
                if self.lattice[x,y,z]==0 and self.has_free_neighbor(x,y,z):
                    break
            lpoly = 0
            while lpoly < lpolymer:
                while True:
                    nxyz = self.box.get_neighbors(x,y,z)
                    nx,ny,nz = nxyz[0],nxyz[1],nxyz[2]
                    if self.lattice[nx,ny,nz]==0:
                        break
                x,y,z=nx,ny,nz
                lpoly+=1
    def monte_carlo_step(self) -> bool:
        """Performs a single Monte Carlo step using the Metropolis algorithm with site exchange."""
        # Select a random True site
        idx = np.random.randint(0, len(self.box.objects))
        object1 = self.box.objects[idx]
        
        site1 = object1.position
        site2 = object1.get_site_to_exchange(self.box)
        object2 = self.box.get_lattice(site2)

        # Compute the energy before the move
        initial_energy = self.box.compute_local_energy(site1) + self.box.compute_local_energy(site2)
        
        # Swap the particles (True <-> False)
        self.box.set_lattice(site1,object2), self.box.set_lattice(site2,object1)
        
        # Compute the energy after the move
        final_energy = self.box.compute_local_energy(site1) + self.box.compute_local_energy(site2)
        
        # Calculate the energy difference
        delta_e = final_energy - initial_energy
        
        # Decide whether to accept the move
        if delta_e > 0 and np.random.rand() >= np.exp(-delta_e / self.T):
            # Reject the move (revert)
            self.box.set_lattice(site1,object1), self.box.set_lattice(site2,object2)
            return False  # Move was rejected

        return True  # Move was accepted
    def monte_carlo_steps(self, steps: int) -> NDArray[np.bool_]:
        success = np.zeros(steps,dtype = np.bool_)
        for step in range(steps):
            success[step] = self.monte_carlo_step()
        return success