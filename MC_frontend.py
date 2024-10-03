# MC.py
import ctypes
import numpy as np
import os
import sys

class MC:
    def __init__(self, size, nparticles, npolymers, lpolymer, interactions, temperature):
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

        # Flatten the interactions matrix
        interactions = np.array(interactions, dtype=np.float32)
        n = interactions.shape[0]
        interactions_flat = interactions.flatten()
        interactions_flat_p = interactions_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

        # Call MC_new to create a new MC instance
        self.Address = self.lib.MC_new(
            size,
            nparticles,
            npolymers,
            lpolymer,
            interactions_flat_p,
            n,
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
    def average_cluster_size(self):
        return self.lib.MC_average_cluster_size(self.Address)
    
