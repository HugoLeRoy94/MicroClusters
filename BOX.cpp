// BOX.cpp
#include "BOX.h"
#include "Objects.h"
#include <stdexcept>
#include <random>
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <cmath>

int to_single_index(int x, int y, int z, int L) {
    return x * L * L + y * L + z;
}

std::tuple<int, int, int> to_xyz(int index, int L) {
    int z = index % L;
    int y = (index / L) % L;
    int x = index / (L * L);
    return std::make_tuple(x, y, z);
}

std::vector<std::array<int, 3>> generate_unique_triplets(int N, int L) {
    int total_triplets = (L + 1) * (L + 1) * (L + 1);

    if (N > total_triplets) {
        throw std::invalid_argument("N is larger than the total number of unique triplets available.");
    }

    std::vector<int> indices(total_triplets);
    for (int i = 0; i < total_triplets; ++i) {
        indices[i] = i;
    }

    static std::mt19937 rng(std::random_device{}());
    std::shuffle(indices.begin(), indices.end(), rng);

    std::vector<std::array<int, 3>> triplets(N);
    for (int i = 0; i < N; ++i) {
        int selected_index = indices[i];

        int z = selected_index % (L + 1);
        int y = (selected_index / (L + 1)) % (L + 1);
        int x = (selected_index / ((L + 1) * (L + 1))) % (L + 1);

        triplets[i] = { x, y, z };
    }
    return triplets;
}

BOX::BOX(int size_, int nobjects, const std::vector<std::vector<float>>& Interactions)
    : size(size_), E(Interactions){
    int total_sites = size * size * size;
    lattice.resize(total_sites);

    // Initialize lattice with Empty objects
    for (int x = 0; x < size; ++x) {
        for (int y = 0; y < size; ++y) {
            for (int z = 0; z < size; ++z) {
                int idx = to_single_index(x, y, z, size);
                lattice[idx] = new Empty(std::make_tuple(x, y, z));
            }
        }
    }

    // Initialize objects array
    objects.resize(nobjects, nullptr);
}

BOX::~BOX() {
    // Delete all unique Objects in lattice
    std::unordered_set<Object*> unique_objects;
    for (auto obj : lattice) {
        unique_objects.insert(obj);
    }
    for (auto obj : unique_objects) {
        delete obj;
    }
}

void BOX::create_new_DHH1(const std::tuple<int,int,int>& site, int object_idx){
    DHH1* new_dhh1 = new DHH1(site);
    int idx = to_single_index(std::get<0>(site), std::get<1>(site), std::get<2>(site), size);
    delete lattice[idx];
    lattice[idx] = new_dhh1;
    objects[object_idx] = new_dhh1;
    new_dhh1->setPosition(site);
}

void BOX::swap(const std::tuple<int, int, int>& site1, const std::tuple<int, int, int>& site2) {
    int idx1 = to_single_index(std::get<0>(site1), std::get<1>(site1), std::get<2>(site1), size);
    int idx2 = to_single_index(std::get<0>(site2), std::get<1>(site2), std::get<2>(site2), size);

    std::swap(lattice[idx1], lattice[idx2]);

    // Update the positions of the objects
    lattice[idx1]->setPosition(site1);
    lattice[idx2]->setPosition(site2);
}

Object* BOX::get_lattice(const std::tuple<int, int, int>& site) const {
    int x = std::get<0>(site);
    int y = std::get<1>(site);
    int z = std::get<2>(site);
    int idx = to_single_index(x, y, z, size);
    return lattice[idx];
}

float BOX::compute_local_energy(const std::tuple<int, int, int>& xyz) const {
    Object* obj = get_lattice(xyz);
    if (obj->isempty()) {
        return 0.0f;
    }
    auto neighbors = get_neighbors(xyz);
    float local_energy = 0.0f;

    for (const auto& nxyz : neighbors) {
        Object* neighbor_obj = get_lattice(nxyz);
        local_energy -= E[obj->Index()][neighbor_obj->Index()];
    }

    return local_energy;
}

float BOX::total_energy() const {
    float energy = 0.0f;
    for (int x = 0; x < size; ++x) {
        for (int y = 0; y < size; ++y) {
            for (int z = 0; z < size; ++z) {
                Object* obj = get_lattice(std::make_tuple(x, y, z));
                if (!obj->isempty()) {
                    energy += compute_local_energy(std::make_tuple(x, y, z));
                }
            }
        }
    }
    return energy / 2.0f;  // Correct for double counting
}

std::vector<std::tuple<int, int, int>> BOX::get_neighbors(const std::tuple<int, int, int>& xyz) const {
    int x = std::get<0>(xyz);
    int y = std::get<1>(xyz);
    int z = std::get<2>(xyz);

    std::vector<std::tuple<int, int, int>> neighbors;
    for (int dx = -1; dx <=1; ++dx) {
        for (int dy = -1; dy <=1; ++dy) {
            for (int dz = -1; dz <=1; ++dz) {
                if (dx == 0 && dy == 0 && dz == 0) {
                    continue;
                }
                int nx = (x + dx + size) % size;
                int ny = (y + dy + size) % size;
                int nz = (z + dz + size) % size;
                neighbors.push_back(std::make_tuple(nx, ny, nz));
            }
        }
    }
    return neighbors;
}

bool BOX::has_free_neighbor(const std::tuple<int, int, int>& xyz) const {
    auto neighbors = get_neighbors(xyz);
    for (const auto& nxyz : neighbors) {
        Object* neighbor_obj = get_lattice(nxyz);
        if (neighbor_obj->isempty()) {
            return true;
        }
    }
    return false;
}

// Placeholder implementations for build_clusters, cluster_size, and compute_av_Nneigh

std::tuple<std::vector<int>, std::vector<int>> BOX::build_clusters() {
    // Implement cluster building algorithm
    return std::make_tuple(std::vector<int>(), std::vector<int>());
}

std::vector<int> BOX::cluster_size() {
    auto [cluster_indices, cluster_starts] = build_clusters();
    std::vector<int> sizes;
    // Compute sizes based on cluster_starts and cluster_indices
    return sizes;
}

int BOX::compute_av_Nneigh() const {
    int av_Nneigh = 0;
    int true_sites = 0;
    for (int x = 0; x < size; ++x) {
        for (int y = 0; y < size; ++y) {
            for (int z = 0; z < size; ++z) {
                Object* obj = get_lattice(std::make_tuple(x, y, z));
                if (!obj->isempty()) {
                    int counter = 0;
                    auto neighbors = get_neighbors(std::make_tuple(x, y, z));
                    for (const auto& neighbor : neighbors) {
                        Object* neighbor_obj = get_lattice(neighbor);
                        if (!neighbor_obj->isempty()) {
                            ++counter;
                        }
                    }
                    av_Nneigh += counter;
                    ++true_sites;
                }
            }
        }
    }
    return true_sites > 0 ? av_Nneigh / true_sites : 0;
}
