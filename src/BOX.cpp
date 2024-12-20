// BOX.cpp
#include "BOX.h"
#include "ComputeLocalEnergy.h"
#include "Objects.h"
#include "Utilities.h"
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
#include <set>

BOX::BOX(int size_, int nobjects, const std::vector<std::vector<float>>& Interactions,double Evalence_,std::mt19937& rng_)
    : size(size_), E(Interactions), Evalence(Evalence_), rng(rng_){
    if (!is_power_of_two(size)){
        throw std::invalid_argument("Size must be a power of 2 for bitmasking optimization.");
    }
    int total_sites = size * size * size;
    lattice.resize(total_sites);

    // Initialize lattice with Empty objects
    for (int idx = 0; idx < total_sites; ++idx) {
        lattice[idx] = std::make_shared<Empty>(idx);
    }

    // Initialize objects array
    objects.reserve(nobjects);
}

void BOX::create_new_DHH1(int index){
    auto new_dhh1 = std::make_shared<DHH1>(index);
    lattice[index] = new_dhh1;
    objects.push_back(new_dhh1);
}

// BOX.cpp

void BOX::add_RNA(int length) {
    int start_index = random_free_site();

    // Build monomer positions
    std::vector<int> monomer_positions;
    monomer_positions.reserve(length);
    monomer_positions.push_back(start_index);

    std::unordered_set<int> monomer_set;
    monomer_set.insert(start_index);

    for (int i = 1; i < length; ++i) {
        std::vector<int> neighbors = get_neighbors(monomer_positions[i - 1]);
        std::shuffle(neighbors.begin(), neighbors.end(), rng);

        bool found = false;
        for (const auto& neighbor_idx : neighbors) {
            if (get_lattice(neighbor_idx)->isempty() && monomer_set.find(neighbor_idx) == monomer_set.end()) {
                monomer_positions.push_back(neighbor_idx);
                monomer_set.insert(neighbor_idx);
                found = true;
                break;
            }
        }

        if (!found) {
            // Clean up and throw exception
            for (const auto& monomer_idx : monomer_positions) {
                lattice[monomer_idx] = std::make_shared<Empty>(monomer_idx);
            }
            throw std::runtime_error("Unable to place RNA polymer due to dead end");
        }
    }

    // Create RNA monomer objects
    auto monomers = std::make_shared<std::vector<int>>(monomer_positions);

    for (int idx = 0; idx < monomers->size(); ++idx) {
        auto rna = std::make_shared<RNA>(monomers, idx);
        objects.push_back(rna);
        lattice[monomers->at(idx)] = rna;
    }
}


std::shared_ptr<Object> BOX::get_lattice(int index) const {
    return lattice[index];
}

void BOX::set_lattice(int index, std::shared_ptr<Object> obj) {
    lattice[index] = obj;
}

/*inline float BOX::compute_local_energy(int index) const {
    auto obj = get_lattice(index);
    if (obj->isempty()) {
        return 0.0f;
    }

    return obj->compute_local_energy(*this);
}*/

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

    //static std::mt19937 rng(std::random_device{}());
    std::shuffle(indices.begin(), indices.end(), rng);

    indices.resize(N);
    return indices;
}

int BOX::random_free_site(){
    //static std::mt19937 rng(std::random_device{}());
    while(true){
        std::uniform_int_distribution<int> dist(0, lattice.size() - 1);
        int idx = dist(rng);
        if(lattice[idx]->isempty()){return idx;}
    }
}

void BOX::check_consistencty(){
    for(auto& object: objects){
        if(lattice.at(object->getPosition())!=object){
            throw std::length_error("object position and lattice inconsistent");
        }
    }
    std::set<int> poses;
    for(auto& object: objects){
        if(std::find(poses.begin(),poses.end(),object->getPosition())==poses.end()){
            poses.insert(object->getPosition());
        }
        else{
            throw std::length_error("two objects at the same position");
        }
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