// MC_wrapper.cpp
#include "MC.h"
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>

extern "C" {

// Function to create a new MC instance.
MC* MC_new(int size, int nparticles, int npolymers, int lpolymer, const float* interactions_flat, int interactions_size, float temperature) {
    // Convert interactions_flat to std::vector<std::vector<float>>
    int n = interactions_size;
    std::vector<std::vector<float>> interactions(n, std::vector<float>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            interactions[i][j] = interactions_flat[i * n + j];
        }
    }
    // Create the MC instance.
    MC* mc = new MC(size, nparticles, npolymers, lpolymer, interactions, temperature);
    return mc;
}

// Function to delete an MC instance.
void MC_delete(MC* mc) {
    delete mc;
}

// Function to perform multiple Monte Carlo steps.
void MC_monte_carlo_steps(MC* mc, int steps, bool* success) {
    std::vector<bool> result = mc->monte_carlo_steps(steps);
    // Copy the result into the success array element by element.
    for (int i = 0; i < steps; ++i) {
        success[i] = result[i];
    }
}

// Function to perform a single Monte Carlo step.
bool MC_monte_carlo_step(MC* mc) {
    return mc->monte_carlo_step();
}

// Function to get the total energy.
float MC_total_energy(MC* mc) {
    return mc->box.total_energy();
}

double MC_average_cluster_size(MC* mc) {
    if (!mc) {
        return -1.0; // Return an error value or handle appropriately
    }
    return mc->box.average_cluster_size();
}

int MC_get_cluster_indices_size(MC* mc) {
    if (!mc) {
        return -1; // Error
    }
    return static_cast<int>(mc->box.get_cluster_indices_size());
}

int MC_get_cluster_starts_size(MC* mc) {
    if (!mc) {
        return -1; // Error
    }
    return static_cast<int>(mc->box.get_cluster_starts_size());
}

int MC_fill_cluster_indices(MC* mc, int* indices_array, int array_length) {
    if (!mc || !indices_array || array_length <= 0) {
        return -1; // Error
    }
    const std::vector<int>& indices = mc->box.get_cluster_indices();
    int size = static_cast<int>(indices.size());
    if (array_length < size) {
        return -1; // Not enough space
    }
    std::copy(indices.begin(), indices.end(), indices_array);
    return size; // Number of elements filled
}

int MC_fill_cluster_starts(MC* mc, int* starts_array, int array_length) {
    if (!mc || !starts_array || array_length <= 0) {
        return -1; // Error
    }
    const std::vector<int>& starts = mc->box.get_cluster_starts();
    int size = static_cast<int>(starts.size());
    if (array_length < size) {
        return -1; // Not enough space
    }
    std::copy(starts.begin(), starts.end(), starts_array);
    return size; // Number of elements filled
}
}