# MC.py
import ctypes
import numpy as np
import os
import sys
import pyvista as pv
def is_power_of_two(n):
    return n > 0 and (n & (n - 1)) == 0

class MC:
    def __init__(self, size, nparticles, npolymers, lpolymer, interactions,Evalence, temperature=1.,seed=98765):
        if not is_power_of_two(size):
            raise ValueError("Size must be a power of 2.")
        # Load the shared library
        if sys.platform.startswith('win'):
            lib_name = 'mc.dll'
        elif sys.platform.startswith('darwin'):
            lib_name = 'libmc.dylib'
        else:
            lib_name = 'libmc.so'

        # Get the directory where the current script is located
        lib_dir = os.path.dirname(os.path.abspath(__file__))

        # Build the full path to the library
        lib_path = os.path.join(lib_dir, lib_name)

        # Load the shared library
        self.lib = ctypes.CDLL(lib_path)

        
        # Define the argument and return types of the functions
        self.lib.MC_new.argtypes = [
            ctypes.c_int,  # size
            ctypes.c_int,  # nparticles
            ctypes.c_int,  # npolymers
            ctypes.c_int,  # lpolymer
            ctypes.POINTER(ctypes.c_float),  # interactions_flat
            ctypes.c_int,  # interactions_size
            ctypes.c_double, # limited valence interaction
            ctypes.c_float,  # temperature
            ctypes.c_int
        ]
        self.lib.MC_new.restype = ctypes.c_void_p  # Return a pointer to MC

        self.lib.MC_delete.argtypes = [ctypes.c_void_p]
        self.lib.MC_delete.restype = None

        self.lib.MC_monte_carlo_steps.argtypes = [
            ctypes.c_void_p,  # MC* mc
            ctypes.c_int,     # int steps
            ctypes.POINTER(ctypes.c_bool)  # bool* success
        ]
        self.lib.MC_monte_carlo_steps.restype = None

        self.lib.MC_monte_carlo_step.argtypes = [ctypes.c_void_p]
        self.lib.MC_monte_carlo_step.restype = ctypes.c_bool

        self.lib.MC_total_energy.argtypes = [ctypes.c_void_p]
        self.lib.MC_total_energy.restype = ctypes.c_float

        self.lib.MC_average_cluster_size.argtypes = [ctypes.c_void_p]
        self.lib.MC_average_cluster_size.restype = ctypes.c_double

        self.lib.MC_get_cluster_indices_size.argtypes = [ctypes.c_void_p]
        self.lib.MC_get_cluster_indices_size.restype = ctypes.c_int

        self.lib.MC_get_cluster_starts_size.argtypes = [ctypes.c_void_p]
        self.lib.MC_get_cluster_starts_size.restype = ctypes.c_int

        self.lib.MC_fill_cluster_indices.argtypes = [
            ctypes.c_void_p,                # MC* mc
            ctypes.POINTER(ctypes.c_int),   # int* indices_array
            ctypes.c_int                    # int array_length
        ]
        self.lib.MC_fill_cluster_indices.restype = ctypes.c_int

        self.lib.MC_fill_cluster_starts.argtypes = [
            ctypes.c_void_p,                # MC* mc
            ctypes.POINTER(ctypes.c_int),   # int* starts_array
            ctypes.c_int                    # int array_length
        ]
        self.lib.MC_fill_cluster_starts.restype = ctypes.c_int

        self.lib.MC_average_cluster_size.argtypes = [ctypes.c_void_p]
        self.lib.MC_average_cluster_size.restype=ctypes.c_double

        self.lib.MC_get_energy.argtypes = [ctypes.c_void_p]
        self.lib.MC_get_energy.restype = ctypes.c_double

        # Flatten the interactions matrix
        interactions = np.array(interactions, dtype=np.float32)
        n = interactions.shape[0]
        interactions_flat = interactions.flatten()
        interactions_flat_p = interactions_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        self.size = size  # Store the size for later use
        
        self.nparticles = nparticles
        self.npolymers = npolymers
        self.lpolymer = lpolymer

        self.lib.MC_fill_DHH1_positions.argtypes = [ctypes.c_void_p,
                                                    ctypes.POINTER(ctypes.c_int),
                                                    ctypes.c_int]
        self.lib.MC_fill_DHH1_positions.restype = ctypes.c_int

        self.lib.MC_fill_RNA_positions.argtypes = [ctypes.c_void_p,
                                                   ctypes.POINTER(ctypes.c_int), ctypes.c_int,
                                                   ctypes.POINTER(ctypes.c_int), ctypes.c_int]
        self.lib.MC_fill_RNA_positions.restype = ctypes.c_int

        # Call MC_new to create a new MC instance
        self.Address = self.lib.MC_new(
            size,
            nparticles,
            npolymers,
            lpolymer,
            interactions_flat_p,
            n,
            Evalence,
            temperature,
            seed
        )

    def __del__(self):
        # Call MC_delete to free the MC instance
        if hasattr(self, 'Address') and self.Address:
            self.lib.MC_delete(self.Address)
            self.Address = None

    def monte_carlo_steps(self, steps):
        # Prepare the success array
        success_array = (ctypes.c_bool * steps)()
        self.lib.MC_monte_carlo_steps(self.Address, steps, success_array)
        # Convert to Python list
        success_list = [bool(success_array[i]) for i in range(steps)]
        return success_list

    def monte_carlo_step(self):
        result = self.lib.MC_monte_carlo_step(self.Address)
        return result

    def total_energy(self):
        energy = self.lib.MC_total_energy(self.Address)
        return energy    
    
    def average_cluster_size(self):
        """Returns the average cluster size in the box."""
        return self.lib.MC_average_cluster_size(self.Address)  
    
    def get_cluster_indices_size(self):
        size = self.lib.MC_get_cluster_indices_size(self.Address)
        if size < 0:
            raise ValueError("Error getting cluster indices size")
        return size
    
    def get_cluster_size(self):
        indices_array,starts_array = self.get_clusters()
        starts_array = np.append(starts_array,indices_array.shape[0])
        return np.diff(starts_array)
    
    def get_cluster_starts_size(self):
        size = self.lib.MC_get_cluster_starts_size(self.Address)
        if size < 0:
            raise ValueError("Error getting cluster starts size")
        return size

    def get_clusters(self):
        """Returns cluster_indices and cluster_starts as NumPy arrays."""
        indices_size = self.get_cluster_indices_size()
        starts_size = self.get_cluster_starts_size()

        if indices_size == 0 or starts_size == 0:
            # No clusters
            return np.array([], dtype=np.int32), np.array([], dtype=np.int32)

        # Allocate NumPy arrays
        indices_array = np.zeros(indices_size, dtype=np.int32)
        starts_array = np.zeros(starts_size, dtype=np.int32)

        # Get pointers to the arrays
        indices_ptr = indices_array.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        starts_ptr = starts_array.ctypes.data_as(ctypes.POINTER(ctypes.c_int))

        # Fill the arrays
        res_indices = self.lib.MC_fill_cluster_indices(self.Address, indices_ptr, indices_size)
        if res_indices != indices_size:
            raise ValueError("Error filling cluster indices")

        res_starts = self.lib.MC_fill_cluster_starts(self.Address, starts_ptr, starts_size)
        if res_starts != starts_size:
            raise ValueError("Error filling cluster starts")

        return indices_array, starts_array
    
    def get_energy(self):
        return self.lib.MC_get_energy(self.Address)
    
    def get_DHH1_positions(self):
        positions_array = np.zeros(self.nparticles, dtype=np.int32)
        positions_ptr = positions_array.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        res = self.lib.MC_fill_DHH1_positions(self.Address, positions_ptr, self.nparticles)
        if res!=self.nparticles:
            raise ValueError("Error filling DHH1 positions")
        return positions_array

    def get_RNA_positions(self):
        if self.npolymers*self.lpolymer <= 0 or self.npolymers <= 0:
            return []

        positions_array = np.zeros(self.npolymers*self.lpolymer, dtype=np.int32)
        lengths_array = np.zeros(self.npolymers, dtype=np.int32)

        positions_ptr = positions_array.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        lengths_ptr = lengths_array.ctypes.data_as(ctypes.POINTER(ctypes.c_int))

        res = self.lib.MC_fill_RNA_positions(self.Address,
                                             positions_ptr, self.npolymers*self.lpolymer,
                                             lengths_ptr, self.npolymers)
        if res != self.npolymers*self.lpolymer:
            print(res)
            print(self.npolymers)
            print(self.lpolymer)
            raise ValueError("Error filling RNA positions")

        # Split positions_array into individual RNA polymers using lengths_array
        rna_positions = []
        pos_index = 0
        for length in lengths_array:
            positions = positions_array[pos_index: pos_index + length]
            rna_positions.append(positions)
            pos_index += length
        
        return rna_positions

            

def to_xyz(index, size):
    z = index % size
    y = (index // size) % size
    x = index // (size * size)
    return np.array([x, y, z])

def minimal_distance(p1, p2, box_size):
    delta = (p2 - p1) % box_size
    delta = np.where(delta > box_size / 2, delta - box_size, delta)
    return np.linalg.norm(delta)

def chebyshev_distance(p1,p2,size):
    x1,y1,z1 = to_xyz(p1,size)
    x2,y2,z2 = to_xyz(p2,size)
    return np.max(np.abs([x2-x1,y2-y1,z2-z1]))

def clip_segment_to_box(p1, p2, size):
    x_min,y_min,z_min=-0.5,-0.5,-0.5
    x_max,y_max,z_max=size-0.5,size-0.5,size-0.5
    x1,y1,z1 = p1[0],p1[1],p1[2]
    x2,y2,z2 = p2[0],p2[1],p2[2]
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    p = [-dx, dx, -dy, dy, -dz, dz]
    q = [x1 - x_min, x_max - x1, y1 - y_min, y_max - y1, z1 - z_min, z_max - z1]
    t_enter = 0.0
    t_exit = 1.0

    for i in range(6):
        if p[i] == 0:  # Check if line is parallel to the clipping boundary
            if q[i] < 0:
                return None  # Line is outside and parallel, so completely discarded
        else:
            t = q[i] / p[i]
            if p[i] < 0:
                if t > t_enter:
                    t_enter = t
            else:
                if t < t_exit:
                    t_exit = t

    if t_enter > t_exit:
        return None  # Line is completely outside

    x1_clip = x1 + t_enter * dx
    y1_clip = y1 + t_enter * dy
    x2_clip = x1 + t_exit * dx
    y2_clip = y1 + t_exit * dy
    z1_clip = z1 + t_enter * dz
    z2_clip = z2 + t_exit * dz


    return np.array([[x1_clip, y1_clip, z1_clip],[x2_clip, y2_clip, z2_clip]])

def generate_lattice_vertices(size):
    # Generate all combinations of x, y, z in 0..size-1
    x = np.arange(size)
    y = np.arange(size)
    z = np.arange(size)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    grid_points = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
    return grid_points

def plot_simulation(mc, output_filename='simulation.html'):
    size = mc.size

    # Get DHH1 positions and convert to coordinates
    dhh1_positions = mc.get_DHH1_positions()
    dhh1_coords = np.array([to_xyz(pos, size) for pos in dhh1_positions])

    # Get RNA positions and convert to coordinates
    rna_positions_list = mc.get_RNA_positions()
    rna_coords_list = []
    for rna_positions in rna_positions_list:
        coords = np.array([to_xyz(pos, size) for pos in rna_positions])
        rna_coords_list.append(coords)

    # Generate lattice vertices
    lattice_vertices = generate_lattice_vertices(size)

    # Create a PyVista plotter
    plotter = pv.Plotter()

    # Add the box contour as a light black wireframe cube
    box_center = (size / 2 - 0.5, size / 2 - 0.5, size / 2 - 0.5)
    box = pv.Cube(center=box_center, x_length=size, y_length=size, z_length=size)
    plotter.add_mesh(box, color='black', opacity=0.2, style='wireframe', line_width=1)

    # Add lattice vertices as small black dots
    #lattice_points = pv.PolyData(lattice_vertices)
    #plotter.add_mesh(lattice_points, color='black', point_size=2.5, render_points_as_spheres=True)

    # Add DHH1 particles as points
    if dhh1_coords.size > 0:
        dhh1_points = pv.PolyData(dhh1_coords)
        plotter.add_mesh(dhh1_points, color='red', point_size=10.0, render_points_as_spheres=True)

    # Plot RNA polymers as segments
    for coords in rna_coords_list:
        if coords.size > 0:
            # Initialize lists to store line segments
            segments = []
            num_monomers = coords.shape[0]
            # plot segment from last to first and first to last, for correct boundary conditions reprsentation
            for i in range(-num_monomers+1,num_monomers - 1):
                if i==-1: #avoid looping segment
                    continue
                p1 = coords[i]
                p2 = coords[i + 1]

                # Unwrap positions to account for periodic boundary conditions
                delta = p2 - p1
                for dim in range(3):
                    if delta[dim] > size / 2:
                        delta[dim] -= size
                    elif delta[dim] < -size / 2:
                        delta[dim] += size
                p2_unwrapped = p1 + delta

                # Check if the segment crosses the boundary
                if np.any(p2_unwrapped < 0) or np.any(p2_unwrapped > size - 1):
                    # Clip the segment at the box boundaries
                    clipped_segment = clip_segment_to_box(p1, p2_unwrapped, size)
                    if clipped_segment is not None:
                        segments.append(clipped_segment)
                else:
                    # Segment is inside the box, add it directly
                    segments.append(np.array([p1, p2_unwrapped]))
            # Plot the segments
            for segment in segments:
                line = pv.Line(segment[0], segment[1])
                plotter.add_mesh(line, color='blue', line_width=3)

            # Also add monomer points
            monomer_points = pv.PolyData(coords)
            plotter.add_mesh(monomer_points, color='blue', point_size=10.0, render_points_as_spheres=True)

    # Set plotter options
    plotter.set_background('white')
    plotter.show_axes()
    plotter.show_bounds(
    grid='front',          # Show grid on the front face
    location='outer',      # Place labels outside the bounding box
    all_edges=True,        # Show all edges of the grid
    xtitle='X Axis',       # Label for the X-axis
    ytitle='Y Axis',       # Label for the Y-axis
    ztitle='Z Axis'        # Label for the Z-axis
    )
    # Optional: Set the camera position
    plotter.view_isometric()

    # Save the plot as an HTML file
    plotter.export_html(output_filename)