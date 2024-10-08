// BOX.cpp
#include "BOX.h"
#include "Objects.h"
#include <stdexcept>
#include <random>
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <cmath>
#include <functional>
#include <stack>
#include <numeric>
#include <vector>

int to_single_index(int x, int y, int z, int L) {
    return x * L * L + y * L + z;
}

std::tuple<int, int, int> to_xyz(int index, int L) {
    int z = index % L;
    int y = (index / L) % L;
    int x = index / (L * L);
    return std::make_tuple(x, y, z);
}

BOX::BOX(int size_, int nobjects, const std::vector<std::vector<float>>& Interactions,double Evalence_)
    : size(size_), E(Interactions), Evalence(Evalence_){
    int total_sites = size * size * size;
    lattice.resize(total_sites);

    // Initialize lattice with Empty objects
    for (int idx = 0; idx < total_sites; ++idx) {
        lattice[idx] = Empty::make_shared_ptr(idx);
    }

    // Initialize objects array
    objects.reserve(nobjects);
}

void BOX::create_new_DHH1(int index){
    auto new_dhh1 = DHH1::make_shared_ptr(index);
    lattice[index] = new_dhh1;
    objects.push_back(new_dhh1);
}

std::shared_ptr<RNA> BOX::add_RNA(int length) {
    int start_index = random_free_site();
    // Initialize data structures
    std::vector<int> monomers;
    monomers.reserve(length);
    monomers.push_back(start_index);

    std::unordered_set<int> monomer_set;
    monomer_set.insert(start_index);

    // Check if the starting position is empty
    if (!get_lattice(start_index)->isempty()) {
        throw std::runtime_error("Starting position is already occupied");
    }
    
    // Initialize RNG   
    static std::mt19937 rng(std::random_device{}());

    // Build the polymer
    for (int i = 1; i < length; ++i) {
        std::vector<int> neighbors = get_neighbors(monomers[i - 1]);
        std::shuffle(neighbors.begin(), neighbors.end(), rng);

        bool found = false;
        for (const auto& neighbor_idx : neighbors) {
            if (get_lattice(neighbor_idx)->isempty() && monomer_set.find(neighbor_idx) == monomer_set.end()) {
                monomers.push_back(neighbor_idx);
                monomer_set.insert(neighbor_idx);
                found = true;
                break;
            }
        }

        if (!found) {
            // Dead end encountered; clean up and throw an exception
            for (const auto& monomer_idx : monomers) {
                lattice[monomer_idx] = Empty::make_shared_ptr(monomer_idx);
            }
            throw std::runtime_error("Unable to place RNA polymer due to dead end");
        }
    }

    // At this point, the monomer positions are determined
    // Now, create the RNA object
    auto rna = RNA::make_shared_ptr(monomers);
    objects.push_back(rna);

    // Update the lattice with the RNA object
    for (const auto& monomer_idx : monomers) {
        lattice[monomer_idx] = rna;
    }

    return rna;
}

void BOX::swap(int idx1, int idx2) {
    std::swap(lattice[idx1], lattice[idx2]);

    // Update the positions of the objects
    lattice[idx1]->setPosition(idx1);
    lattice[idx2]->setPosition(idx2);

    clusters_valid = false; // Invalidate cluster data
}

std::shared_ptr<Object> BOX::get_lattice(int index) const {
    return lattice[index];
}

void BOX::set_lattice(int index, std::shared_ptr<Object> obj) {
    lattice[index] = obj;
}

float BOX::compute_local_energy(int index) const {
    auto obj = get_lattice(index);
    if (obj->isempty()) {
        return 0.0f;
    }

    return obj->compute_local_energy(*this);
}

float BOX::total_energy() const {
    float energy = 0.0f;
    for (int idx = 0; idx < size * size * size; ++idx) {
        auto obj = get_lattice(idx);
        if (!obj->isempty()) {
            energy += obj->compute_local_energy(*this);
        }
    }
    return energy / 2.0f;  // Correct for double counting
}

std::vector<int> BOX::get_neighbors(int index) const {
    int x, y, z;
    std::tie(x, y, z) = to_xyz(index, size);

    std::vector<int> neighbors;
    neighbors.reserve(26);
    for (int dx = -1; dx <=1; ++dx) {
        for (int dy = -1; dy <=1; ++dy) {
            for (int dz = -1; dz <=1; ++dz) {
                if (dx == 0 && dy == 0 && dz == 0) {
                    continue;
                }
                int nx = (x + dx + size) % size;
                int ny = (y + dy + size) % size;
                int nz = (z + dz + size) % size;
                int neighbor_idx = to_single_index(nx, ny, nz, size);
                neighbors.push_back(neighbor_idx);
            }
        }
    }
    return neighbors;
}

bool BOX::has_free_neighbor(int index) const {
    auto neighbors = get_neighbors(index);
    for (const auto& nidx : neighbors) {
        auto neighbor_obj = get_lattice(nidx);
        if (neighbor_obj->isempty()) {
            return true;
        }
    }
    return false;
}

// Implement cluster-related methods similarly, adjusting for index-based positions

std::vector<int> BOX::generate_unique_indices(int N) {
    int total_sites = size * size * size;
    if (N > total_sites) {
        throw std::invalid_argument("N is larger than the total number of sites available.");
    }

    std::vector<int> indices(total_sites);
    for (int i = 0; i < total_sites; ++i) {
        indices[i] = i;
    }

    static std::mt19937 rng(std::random_device{}());
    std::shuffle(indices.begin(), indices.end(), rng);

    indices.resize(N);
    return indices;
}

int BOX::random_free_site(){
    static std::mt19937 rng(std::random_device{}());
    while(true){
        std::uniform_int_distribution<int> dist(0, lattice.size() - 1);
        int idx = dist(rng);
        if(lattice[idx]->isempty()){return idx;}
    }
}
// Placeholder implementations for build_clusters, cluster_size, and compute_av_Nneigh
// Clustering Functions

// Build clusters using Depth-First Search (DFS)
void BOX::build_clusters() {
    int total_sites = size * size * size;
    std::vector<bool> visited(total_sites, false);

    cluster_indices.clear();
    cluster_indices.reserve(total_sites);

    cluster_starts.clear();
    cluster_starts.reserve(total_sites);

    int current_cluster_start = 0;

    // Define the flood_fill function using a lambda
    std::function<void(int)> flood_fill = [&](int start_index) {
        std::stack<int> stack;
        stack.push(start_index);
        cluster_starts.push_back(static_cast<int>(cluster_indices.size()));

        while (!stack.empty()) {
            int current_idx = stack.top();
            stack.pop();

            if (!visited[current_idx]) {
                visited[current_idx] = true;
                cluster_indices.push_back(current_idx);

                auto neighbors = get_neighbors(current_idx);
                for (const auto& neighbor_idx : neighbors) {
                    auto neighbor_obj = get_lattice(neighbor_idx);
                    if (!neighbor_obj->isempty() && !visited[neighbor_idx]) {
                        stack.push(neighbor_idx);
                    }
                }
            }
        }
    };

    // Iterate over all lattice sites
    for (int idx = 0; idx < total_sites; ++idx) {
        if (!visited[idx]) {
            auto obj = get_lattice(idx);
            if (!obj->isempty()) {
                // Start a new cluster
                flood_fill(idx);
            } else {
                visited[idx] = true;
            }
        }
    }

    clusters_valid = true; // Mark clusters as valid
}

// Compute sizes of all clusters
std::vector<int> BOX::cluster_size() {
    if (!clusters_valid) {
        build_clusters();
    }

    std::vector<int> sizes;
    sizes.reserve(cluster_starts.size());

    for (size_t i = 0; i < cluster_starts.size(); ++i) {
        if (i < cluster_starts.size() - 1) {
            sizes.push_back(cluster_starts[i + 1] - cluster_starts[i]);
        } else {
            sizes.push_back(static_cast<int>(cluster_indices.size()) - cluster_starts[i]);
        }
    }

    return sizes;
}

// Compute the average cluster size
double BOX::average_cluster_size() {
    std::vector<int> sizes = cluster_size();
    if (sizes.empty()) {
        return 0.0;
    }
    double total_size = std::accumulate(sizes.begin(), sizes.end(), 0.0);
    return total_size / static_cast<double>(sizes.size());
}

// Compute the average number of neighbors for occupied sites
double BOX::compute_av_Nneigh() const {
    double total_neighbors = 0.0;
    int occupied_sites = 0;

    for (int idx = 0; idx < size * size * size; ++idx) {
        auto obj = get_lattice(idx);
        if (!obj->isempty()) {
            int neighbor_count = 0;
            auto neighbors = get_neighbors(idx);
            for (const auto& neighbor_idx : neighbors) {
                auto neighbor_obj = get_lattice(neighbor_idx);
                if (!neighbor_obj->isempty()) {
                    ++neighbor_count;
                }
            }
            total_neighbors += neighbor_count;
            ++occupied_sites;
        }
    }

    return occupied_sites > 0 ? (total_neighbors / static_cast<double>(occupied_sites)) : 0.0;
}

// Get cluster indices
const std::vector<int>& BOX::get_cluster_indices() {
    if (!clusters_valid) {
        build_clusters();
    }
    return cluster_indices;
}

// Get cluster starts
const std::vector<int>& BOX::get_cluster_starts() {
    if (!clusters_valid) {
        build_clusters();
    }
    return cluster_starts;
}

// Get the size of cluster_indices
size_t BOX::get_cluster_indices_size() {
    if (!clusters_valid) {
        build_clusters();
    }
    return cluster_indices.size();
}

// Get the size of cluster_starts
size_t BOX::get_cluster_starts_size() {
    if (!clusters_valid) {
        build_clusters();
    }
    return cluster_starts.size();
}