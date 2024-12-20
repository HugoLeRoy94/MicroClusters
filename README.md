# Monte Carlo Simulation of Polymers and Particles

This project provides a Monte Carlo simulation of polymers and particles in a 3D lattice, written in C++ with a Python interface. The simulation allows you to model interactions between different objects (e.g., polymers and particles) and compute properties like total energy and cluster sizes.

## Table of Contents

- Overview
- Features
- Requirements
- Compilation
- Python Interface
  - Installation
  - Usage
- Examples
  - Running a Simulation
  - Analyzing Results
- API Reference
  - MC Class
- Contributing
- License

## Overview

This simulation models the behavior of polymers and particles in a 3D lattice using Monte Carlo methods. The primary components are:

- **Polymers**: Represented as chains of connected monomers.
- **Particles**: Individual entities that can interact with polymers and other particles.
- **Lattice**: A 3D grid where polymers and particles reside.

Interactions between objects are defined via an interaction matrix, allowing for customizable simulations.

## Features

- **Monte Carlo Simulation**: Perform Monte Carlo steps to simulate the movement and interaction of objects.
- **Energy Calculation**: Compute the total energy of the system based on interactions.
- **Cluster Analysis**: Identify and analyze clusters formed by objects in the lattice.
- **Python Interface**: Interact with the simulation using a Python class for ease of use.

## Requirements

- **C++ Compiler**: Supports C++11 or later.
- **Python**: Version 3.x.
- **NumPy**: For handling numerical arrays in Python.
- **CMake** (optional): For building the project.

## Compilation

To use the simulation from Python, you need to compile the C++ code into a shared library that can be loaded by Python.

### Steps:

1. **Clone the Repository**:

   ```bash
   git clone https://github.com/HugoLeRoy94/MicroClusters.git
   cd MicroClusters
   ```

2. **Organize the Files**:

   Ensure the following files are in the same directory:

   - `MC.h`
   - `MC.cpp`
   - `BOX.h`
   - `BOX.cpp`
   - `Objects.h`
   - `Objects.cpp`
   - `Move.h`
   - `Move.cpp`
   - `MC_front.cpp`
   - `MC.py`

3. **Compile the Shared Library**:

   Depending on your operating system, run the appropriate command.

   - **Linux**:

     ```bash
     make
     ```

   - **macOS**:

     ```bash
     g++ -shared -dynamiclib -std=c++11 MC.cpp BOX.cpp Objects.cpp Move.cpp MC_front.cpp -o libmc.dylib
     ```

   - **Windows**:

     ```bash
     g++ -shared -std=c++11 MC.cpp BOX.cpp Objects.cpp Move.cpp MC_front.cpp -o mc.dll
     ```

   Ensure that you have `g++` installed and accessible from your command line. You might need to adjust the compiler flags based on your environment.

4. **Verify the Library**:

   After compilation, you should have a shared library file (`libmc.so`, `libmc.dylib`, or `mc.dll`) in your directory.

## Python Interface

The Python interface allows you to use the simulation conveniently within Python scripts or Jupyter notebooks.

### Installation

1. **Place the Shared Library**:

   Ensure that the compiled shared library is in the same directory as `MC.py`.

2. **Install Dependencies**:

   make sure that the modules in requirements.txt are all installed

   ```bash
   pip install -r requirements.txt
   ```

3. **Set Up Environment Variables** (Optional):

   If the shared library is not in the same directory, you may need to set the `LD_LIBRARY_PATH` (Linux), `DYLD_LIBRARY_PATH` (macOS), or add the directory to your system's PATH (Windows).

### Usage

Import the `MC` class from `MC.py` in your Python script:

```python
import sys
sys.path.append('src/') # the path to the library
from MC_wrapper import MC
```

## Examples

### Running a Simulation

Here's a step-by-step example of how to set up and run a simulation.

#### 1. **Import Necessary Modules**:

```python
import numpy as np
import sys
sys.path.append('src/')
import MC_wrapper
```

#### 2. **Define Simulation Parameters**:

```python
# Lattice size must be a power of 2
size = 32

# Number of particles and polymers
nparticles = 100
npolymers = 50
lpolymer = 10  # Length of each polymer

# Interaction matrix (example with 3 object types: Empty, DHH1, RNA)
# Positive interaction value corresponds to attractive force
# The matrix should be symmetric with zero interaction for Empty (index 0)
interactions = [
    [0.0, 0.0, 0.0],    # Interactions with Empty
    [0.0, 1.0, 0.],  # Interactions for DHH1 (index 1)
    [0.0, 0, 1.0]   # Interactions for RNA (index 2)
]

# Valence energy (additional energy term)
Evalence = -0.3

# Temperature set the energy unit
temperature = 1.0

# seed to make reproduceable simulations
seed = 9876

# ratio of diffusive vs conformational moves for the RNA
# Careful, the diffusives moves slow-down the simulation by a lot
diff_moves_ratio = 0.1
```

#### 3. **Initialize the Simulation**:

```python
mc = MC(size, nparticles, npolymers, lpolymer, interactions, Evalence, temperature, seed,diff_moves_ratio)
```

#### 4. **Perform Monte Carlo Steps**:

```python
# Number of Monte Carlo steps
steps = 1000

# Run the simulation
mc.monte_carlo_steps(steps)
```

#### 5. **Retrieve Results**:

```python
# Get total energy
energy = mc.get_energy()
print(f"Total Energy: {energy}")

# Get average cluster size
avg_cluster_size = mc.average_cluster_size()
print(f"Average Cluster Size: {avg_cluster_size}")

# Get cluster details
indices_array, starts_array = mc.get_clusters()
# instead of making a list of list, clusters are represented as a single list,
# and a list of start indices.
# indices_array gives the position of the particles
print(f"Cluster Indices: {indices_array}")
print(f"Cluster Starts: {starts_array}")

# Get sizes of clusters
cluster_sizes = mc.get_cluster_size()
print(f"Cluster Sizes: {cluster_sizes}")
```

### Analyzing Results

You can use NumPy and Matplotlib to analyze and visualize the simulation results.

```python
import matplotlib.pyplot as plt

# Histogram of cluster sizes
plt.hist(cluster_sizes, bins=range(1, max(cluster_sizes) + 1))
plt.xlabel('Cluster Size')
plt.ylabel('Frequency')
plt.title('Distribution of Cluster Sizes')
plt.show()
```

### Plot the results

```python
# create a plot
plotter = MC_wrapper.create_plotter(mc.size)
# add a system to plot in
MC_wrapper.add_system(mc,plotter)
# export to a file
MC_wrapper.export_to_html(plotter, 'simulation.html')
```

This is equivalent to run:
```python
Mc_wrapper.plot_simulation(mc)
```

However, it is also possible to plot several systems:

```python
plotter = MC_wrapper.create_plotter(mc.size)
MC_wrapper.add_system(mc,plotter)
for step in range(10):
    mc.monte_carlo_steps(100)
    # plot the system every 100 steps with a gradient color-coding
    MC_wrapper.add_system(mc,plotter,dhh1_color=viridis(step/10.),rna_color=plasma(step/10.))
# export the whole plot
MC_wrapper.export_to_html(plotter, 'simulation.html')
```

## API Reference

### MC Class

The `MC` class is the primary interface to the simulation.

#### Constructor

```python
MC(size, nparticles, npolymers, lpolymer, interactions, Evalence, temperature,seed, diff_moves_ratio)
```

- **Parameters**:
  - `size` (int): Size of the lattice (the lattice is size x size x size).
  - `nparticles` (int): Number of particles to place in the lattice.
  - `npolymers` (int): Number of polymers to place in the lattice.
  - `lpolymer` (int): Length of each polymer.
  - `interactions` (list of lists): Interaction matrix between different object types.
  - `Evalence` (float): Valence energy term.
  - `temperature` (float): Temperature of the system.
  - `seed` (int): seed to make reproducible simulations
  - `diff_moves_ratio` (float) : percentage of diffusive to conformational moves for RNA

#### Methods

- `monte_carlo_step()`: Performs a single Monte Carlo step. Returns `True` if the move was accepted.
- `monte_carlo_steps(steps)`: Performs multiple Monte Carlo steps. Returns a list of booleans indicating whether each move was accepted.
- `get_energy()`: Returns the total energy of the system.
- `average_cluster_size()`: Returns the average cluster size.
- `get_clusters()`: Returns two NumPy arrays, `indices_array` and `starts_array`, representing clusters.
- `get_cluster_size()`: Returns a NumPy array of cluster sizes.

#### Properties

- **Address**: Memory address of the C++ MC instance (used internally).

## Contributing

Contributions are welcome! Please submit a pull request or open an issue to discuss your ideas.

## License

This project is licensed under the MIT License.

---

**Note**: Ensure that the compiled shared library (`libmc.so`, `libmc.dylib`, or `mc.dll`) is accessible to your Python script, either by placing it in the same directory or by adjusting your system's library path.