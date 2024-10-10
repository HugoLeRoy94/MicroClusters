# MC.py
import ctypes
import numpy as np
import os
import sys
import pyvista as pv

class MC:
    def __init__(self, size, nparticles, npolymers, lpolymer, interactions,Evalence, temperature):
        # Load the shared library
        if sys.platform.startswith('win'):
            lib_name = 'mc.dll'
        elif sys.platform.startswith('darwin'):
            lib_name = 'libmc.dylib'
        else:
            lib_name = 'libmc.so'

        lib_path = os.path.abspath(lib_name)
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
            ctypes.c_float  # temperature
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
            temperature
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

    # Create a PyVista plotter
    plotter = pv.Plotter()

    # Add the box contour as a light black wireframe cube
    box_center = (size / 2 - 0.5, size / 2 - 0.5, size / 2 - 0.5)
    box = pv.Cube(center=box_center, x_length=size, y_length=size, z_length=size)
    plotter.add_mesh(box, color='black', opacity=0.2, style='wireframe', line_width=1)

    # Add DHH1 particles as points
    if dhh1_coords.size > 0:
        dhh1_points = pv.PolyData(dhh1_coords)
        plotter.add_mesh(dhh1_points, color='red', point_size=10.0, render_points_as_spheres=True)

    # Set threshold slightly above 1 to account for numerical errors
    threshold = 1.1

    # Add RNA polymers
    #for coords in rna_coords_list:
    #    if len(coords) >= 2:
    #        segments = []
    #        current_segment = [coords[0]]
    #        for i in range(1, len(coords)):
    #            prev = coords[i - 1]
    #            curr = coords[i]
    #            #dist = minimal_distance(prev, curr, size)
    #            dist = chebyshev_distance(prev,curr,size)
    #            if dist < threshold:
    #                current_segment.append(curr)
    #            else:
    #                # Finalize current segment
    #                if len(current_segment) >= 1:
    #                    segments.append(np.array(current_segment))
    #                current_segment = [curr]
    #        # Add the last segment
    #        if len(current_segment) >= 1:
    #            segments.append(np.array(current_segment))
    #        # Plot all segments
    #        for segment in segments:
    #            if len(segment) >= 2:
    #                lines = pv.lines_from_points(segment)
    #                plotter.add_mesh(lines, color='blue', line_width=2.0)
    #    elif len(coords) == 1:
    #        # Single monomer (will plot as point below)
    #        pass  # We will plot all monomers as points next

    # Plot all RNA monomers as points on top of the segments
    for coords in rna_coords_list:
        if coords.size > 0:
            monomer_points = pv.PolyData(coords)
            plotter.add_mesh(monomer_points, color='blue', point_size=5.0, render_points_as_spheres=True)

    # Set plotter options
    plotter.set_background('white')
    plotter.show_axes()

    # Optional: Set the camera position
    plotter.view_isometric()

    # Save the plot as an HTML file
    plotter.export_html(output_filename)