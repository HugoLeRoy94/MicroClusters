import numpy as np
import numba

# Parameters
L = 10  # Size of the 3D box (LxLxL)
n_steps = 10000  # Number of Monte Carlo steps
temperature = 1.0  # Temperature of the system
interaction_strength = 1.0  # Interaction strength scale

@numba.njit
def get_neighbors(x, y, z, L):
    """ Returns the coordinates of the 8 nearest neighbors of a given lattice site (x, y, z). """
    neighbors = []
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dz in [-1, 0, 1]:
                if dx == 0 and dy == 0 and dz == 0:
                    continue
                neighbors.append(((x + dx) % L, (y + dy) % L, (z + dz) % L))
    return neighbors

@numba.njit
def compute_local_energy(lattice, x, y, z, L, interaction_strength):
    """ Computes the local energy for a given lattice site. """
    neighbors = get_neighbors(x, y, z, L)
    local_energy = 0
    
    for nx, ny, nz in neighbors:
        if lattice[nx, ny, nz] == 1:  # Interaction with other particles
            local_energy -= interaction_strength
    
    # Specific interaction with at most 2 peers
    if lattice[x, y, z] == 1:
        specific_neighbors = np.random.choice(len(neighbors), 2, replace=False)
        for idx in specific_neighbors:
            nx, ny, nz = neighbors[idx]
            if lattice[nx, ny, nz] == 1:
                local_energy -= interaction_strength

    return local_energy

@numba.njit
def total_energy(lattice, L, interaction_strength):
    """ Computes the total energy of the system. """
    energy = 0
    for x in range(L):
        for y in range(L):
            for z in range(L):
                if lattice[x, y, z] == 1:
                    energy += compute_local_energy(lattice, x, y, z, L, interaction_strength)
    return energy / 2  # To correct for double counting of pairs

@numba.njit
def monte_carlo_step(lattice, L, interaction_strength, temperature):
    """ Performs a single Monte Carlo step using the Metropolis algorithm. """
    # Select a random site
    x, y, z = np.random.randint(0, L, size=3)
    
    # Compute the energy before the move
    initial_energy = compute_local_energy(lattice, x, y, z, L, interaction_strength)
    
    # Flip the particle (0 -> 1 or 1 -> 0)
    lattice[x, y, z] = 1 - lattice[x, y, z]
    
    # Compute the energy after the move
    final_energy = compute_local_energy(lattice, x, y, z, L, interaction_strength)
    
    # Calculate the energy difference
    delta_e = final_energy - initial_energy
    
    # Decide whether to accept the move
    if delta_e > 0 and np.random.rand() >= np.exp(-delta_e / temperature):
        lattice[x, y, z] = 1 - lattice[x, y, z]  # Reject the move (revert)

@numba.njit
def run_simulation(L, n_steps, temperature, interaction_strength):
    """ Runs the Monte Carlo simulation. """
    # Initialize lattice
    lattice = np.random.choice([0, 1], size=(L, L, L))
    
    # Compute initial energy
    energy = total_energy(lattice, L, interaction_strength)
    print(f"Initial Energy: {energy}")
    
    # Perform Monte Carlo steps
    for step in range(n_steps):
        monte_carlo_step(lattice, L, interaction_strength, temperature)
        
        # Optionally, track energy (e.g., print every 1000 steps)
        if step % 1000 == 0:
            energy = total_energy(lattice, L, interaction_strength)
            print(f"Step {step}, Energy: {energy}")
    
    # Final energy after simulation
    final_energy = total_energy(lattice, L, interaction_strength)
    print(f"Final Energy: {final_energy}")

# Run the simulation
run_simulation(L, n_steps, temperature, interaction_strength)
