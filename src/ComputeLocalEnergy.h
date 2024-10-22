// ComputeLocalEnergy.h
#ifndef COMPUTE_LOCAL_ENERGY_H
#define COMPUTE_LOCAL_ENERGY_H

#include "BOX.h"
#include "Objects.h"

inline std::vector<int> BOX::get_neighbors(int index) const {
    int x, y, z;
    std::tie(x, y, z) = to_xyz(index,size);
    int mask = size - 1;

    std::vector<int> neighbors;
    neighbors.reserve(26);

    for (int dx = -1; dx <= 1; ++dx) {
        int nx = (x + dx) & mask;
        for (int dy = -1; dy <= 1; ++dy) {
            int ny = (y + dy) & mask;
            for (int dz = -1; dz <= 1; ++dz) {
                int nz = (z + dz) & mask;
                if (dx == 0 && dy == 0 && dz == 0) {
                    continue;
                }
                int neighbor_idx = to_single_index(nx, ny, nz,size);
                neighbors.push_back(neighbor_idx);
            }
        }
    }
    return neighbors;
}

// Inline definition of BOX::compute_local_energy
inline float BOX::compute_local_energy(int index) const {
    auto obj = get_lattice(index);
    if (obj->isempty()) {
        return 0.0f;
    }
    return obj->compute_local_energy(*this);
}

// Inline definitions of Object's compute_local_energy methods

// For RNA
inline float RNA::compute_local_energy(const BOX& box) const {
    float local_energy = 0.0f;
    //for (const auto& idx : monomers) {
    auto neighbors = box.get_neighbors(monomers->at(index));
    for (const auto& nidx : neighbors) {
        auto neighbor_obj = box.get_lattice(nidx);
        local_energy -= box.E[Index()][neighbor_obj->Index()];
    }
    //}
    return local_energy;
}

// For DHH1
inline float DHH1::compute_local_energy(const BOX& box) const {
    auto neighbors = box.get_neighbors(getPosition());
    float local_energy = 0.0f;
    bool has_neigh = false;
    for (const auto& nidx : neighbors) {
        auto neighbor_obj = box.get_lattice(nidx);
        local_energy -= box.E[Index()][neighbor_obj->Index()];
        if (neighbor_obj->Index() == 1) {
            has_neigh = true;
        }
    }
    if (has_neigh) {
        local_energy -= box.Evalence;
    }
    return local_energy;
}

// For Empty
inline float Empty::compute_local_energy(const BOX& box) const {
    return 0.0f; // Empty objects have zero energy
}

#endif // COMPUTE_LOCAL_ENERGY_H
