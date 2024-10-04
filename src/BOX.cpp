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

struct TupleHash {
    std::size_t operator()(const std::tuple<int, int, int>& t) const {
        auto h1 = std::hash<int>{}(std::get<0>(t));
        auto h2 = std::hash<int>{}(std::get<1>(t));
        auto h3 = std::hash<int>{}(std::get<2>(t));
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

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

BOX::BOX(int size_, int nobjects, const std::vector<std::vector<float>>& Interactions,double Evalence_)
    : size(size_), E(Interactions), Evalence(Evalence_){
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

void BOX::create_new_DHH1(const std::tuple<int,int,int>& site){
    DHH1* new_dhh1 = new DHH1(site);
    int idx = to_single_index(std::get<0>(site), std::get<1>(site), std::get<2>(site), size);
    delete lattice[idx];
    lattice[idx] = new_dhh1;
    objects.push_back(new_dhh1);
    new_dhh1->setPosition(site);
}
RNA* BOX::add_RNA(int length) {
    std::tuple<int,int,int> start_position(random_free_site());
    // Initialize data structures
    std::vector<std::tuple<int, int, int>> monomers;
    monomers.reserve(length);
    monomers.push_back(start_position);

    std::unordered_set<std::tuple<int, int, int>, TupleHash> monomer_set;
    monomer_set.insert(start_position);

    // Check if the starting position is empty
    if (!get_lattice(start_position)->isempty()) {
        throw std::runtime_error("Starting position is already occupied");
    }
    

    // Initialize RNG   
    static std::mt19937 rng(std::random_device{}());

    // Build the polymer
    for (int i = 1; i < length; ++i) {
        std::vector<std::tuple<int, int, int>> neighbors = get_neighbors(monomers[i - 1]);
        std::shuffle(neighbors.begin(), neighbors.end(), rng);

        bool found = false;
        for (const auto& neighbor : neighbors) {
            if (get_lattice(neighbor)->isempty() && monomer_set.find(neighbor) == monomer_set.end()) {
                monomers.push_back(neighbor);
                monomer_set.insert(neighbor);
                found = true;
                break;
            }
        }

        if (!found) {
            // Dead end encountered; clean up and throw an exception or return nullptr
            for (const auto& monomer_pos : monomers) {
                lattice[to_single_index(std::get<0>(monomer_pos),std::get<1>(monomer_pos),std::get<2>(monomer_pos),size)]= new Empty(monomer_pos);
            }
            throw std::runtime_error("Unable to place RNA polymer due to dead end");
        }
    }

    // At this point, the monomer positions are determined
    // Now, create the RNA object
    RNA* rna = new RNA(monomers);
    objects.push_back(rna);

    // Update the lattice with the RNA object
    for (const auto& monomer_pos : monomers) {
        delete lattice[to_single_index(std::get<0>(monomer_pos),std::get<1>(monomer_pos),std::get<2>(monomer_pos),size)];
        lattice[to_single_index(std::get<0>(monomer_pos),std::get<1>(monomer_pos),std::get<2>(monomer_pos),size)]=rna;
    }

    return rna;
}
void BOX::swap(const std::tuple<int, int, int>& site1, const std::tuple<int, int, int>& site2) {
    int idx1 = to_single_index(std::get<0>(site1), std::get<1>(site1), std::get<2>(site1), size);
    int idx2 = to_single_index(std::get<0>(site2), std::get<1>(site2), std::get<2>(site2), size);

    std::swap(lattice[idx1], lattice[idx2]);

    // Update the positions of the objects
    lattice[idx1]->setPosition(site1);
    lattice[idx2]->setPosition(site2);

    clusters_valid = false; // Invalidate cluster data
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

    return obj->compute_local_energy(*this);
}

float BOX::total_energy() const {
    float energy = 0.0f;
    for (int x = 0; x < size; ++x) {
        for (int y = 0; y < size; ++y) {
            for (int z = 0; z < size; ++z) {
                Object* obj = get_lattice(std::make_tuple(x, y, z));
                if (!obj->isempty()) {
                    energy += obj->compute_local_energy(*this);
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

void BOX::build_clusters() {
    int total_sites = size * size * size;
    std::vector<bool> visited(total_sites, false);

    cluster_indices.clear();
    cluster_indices.reserve(total_sites);

    cluster_starts.clear();
    cluster_starts.reserve(total_sites);

    int current_index = 0;

    // Define the flood_fill function
    std::function<void(int)> flood_fill = [&](int start_index) {
        std::stack<int> stack;
        stack.push(start_index);
        int cluster_start = current_index;

        while (!stack.empty()) {
            int current_idx = stack.top();
            stack.pop();
            if (!visited[current_idx]) {
                visited[current_idx] = true;
                cluster_indices.push_back(current_idx);
                current_index++;
                int x, y, z;
                std::tie(x, y, z) = to_xyz(current_idx, size);
                auto neighbors = get_neighbors(std::make_tuple(x, y, z));
                for (const auto& neighbor : neighbors) {
                    int nx = std::get<0>(neighbor);
                    int ny = std::get<1>(neighbor);
                    int nz = std::get<2>(neighbor);
                    int neighbor_idx = to_single_index(nx, ny, nz, size);
                    if (!visited[neighbor_idx]) {
                        Object* neighbor_obj = get_lattice(neighbor);
                        if (!neighbor_obj->isempty()) {
                            stack.push(neighbor_idx);
                        }
                    }
                }
            }
        }
        if (current_index > cluster_start) {
            cluster_starts.push_back(cluster_start);
        }
    };

    // Iterate over all lattice sites
    for (int x = 0; x < size; ++x) {
        for (int y = 0; y < size; ++y) {
            for (int z = 0; z < size; ++z) {
                int idx = to_single_index(x, y, z, size);
                if (!visited[idx]) {
                    Object* obj = get_lattice(std::make_tuple(x, y, z));
                    if (!obj->isempty()) {
                        // Start a new cluster
                        flood_fill(idx);
                    } else {
                        visited[idx] = true;
                    }
                }
            }
        }
    }
    clusters_valid = true; // Mark clusters as valid
}

std::vector<int> BOX::cluster_size() {
    if (!clusters_valid) {
        build_clusters();
    }
    // Append length of cluster_indices to cluster_starts to mark the end
    //cluster_starts.push_back(cluster_indices.size());
    // Compute sizes as differences between consecutive cluster_starts
    std::vector<int> sizes;
    sizes.reserve(cluster_starts.size()-1);
    for (size_t i = 0; i < cluster_starts.size() - 1; ++i) {
        int size = cluster_starts[i + 1] - cluster_starts[i];
        sizes.push_back(size);
    }
    sizes.push_back(cluster_indices.size()-cluster_starts.back());
    return sizes;
}

double BOX::average_cluster_size() {
    std::vector<int> sizes = cluster_size();
    if (sizes.empty()) {
        return 0.0;
    }
    double total_size = std::accumulate(sizes.begin(), sizes.end(), 0.0);
    return total_size / sizes.size();
}

double BOX::compute_av_Nneigh() const {
    double total_neighbors = 0.0;
    int occupied_sites = 0;

    for (int x = 0; x < size; ++x) {
        for (int y = 0; y < size; ++y) {
            for (int z = 0; z < size; ++z) {
                Object* obj = get_lattice(std::make_tuple(x, y, z));
                if (!obj->isempty()) {
                    int neighbor_count = 0;
                    auto neighbors = get_neighbors(std::make_tuple(x, y, z));
                    for (const auto& neighbor : neighbors) {
                        Object* neighbor_obj = get_lattice(neighbor);
                        if (!neighbor_obj->isempty()) {
                            ++neighbor_count;
                        }
                    }
                    total_neighbors += neighbor_count;
                    ++occupied_sites;
                }
            }
        }
    }

    return occupied_sites > 0 ? total_neighbors / occupied_sites : 0.0;
}

const std::vector<int>& BOX::get_cluster_indices() {
    if (!clusters_valid) {
        build_clusters();
    }
    return cluster_indices;
}

const std::vector<int>& BOX::get_cluster_starts() {
    if (!clusters_valid) {
        build_clusters();
    }
    return cluster_starts;
}

size_t BOX::get_cluster_indices_size() {
    if (!clusters_valid) {
        build_clusters();
    }
    return cluster_indices.size();
}

size_t BOX::get_cluster_starts_size() {
    if (!clusters_valid) {
        build_clusters();
    }
    return cluster_starts.size();
}
std::tuple<int,int,int> BOX::random_free_site(){
    static std::mt19937 rng(std::random_device{}());
    while(true){
        std::uniform_int_distribution<int> dist(0, lattice.size() - 1);
        int idx = dist(rng);
        if(lattice[idx]->isempty()){return to_xyz(idx,size);}
    }
}