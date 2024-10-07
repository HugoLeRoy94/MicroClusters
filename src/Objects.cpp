// Objects.cpp
#include "Objects.h"
#include "BOX.h"
#include <random>
#include <stdexcept>


// Object Class Implementation
Object::Object(const std::tuple<int, int, int>& position_) : position(position_) {}
Object::~Object() {}

const std::tuple<int, int, int>& Object::getPosition() const { return position; }
std::vector<std::tuple<int, int, int>> Object::get_positions() const {
    return {getPosition()};
}
void Object::setPosition(const std::tuple<int, int, int>& value) { position = value; }
bool Object::isempty() const { return false; }
int Object::Index() const { return 0; }
std::tuple<int, int, int> Object::get_site_to_exchange(const BOX& box) const { return position; }

// Empty Class Implementation
Empty::Empty(const std::tuple<int, int, int>& position_) : Object(position_) {}
Empty::~Empty() {}
bool Empty::isempty() const { return true; }
int Empty::Index() const { return 0; }
float Empty::compute_local_energy(const BOX& box) const {
    // Empty objects have zero energy
    return 0.0f;
}


// RNA Class Implementation
RNA::RNA(std::vector<std::tuple<int,int,int>>& monomers_)
    : Object(monomers_[0]) {
        monomers = monomers_;
}

RNA::~RNA() {}
int RNA::Index() const { return 2; }
bool RNA::isconnected(int idx, const BOX& box) const {
    // Ensure idx is within bounds
    if (idx < 0 || idx >= monomers.size()) {
        throw std::out_of_range("Index out of bounds in RNA::isconnected");
    }

    const auto& current_monomer = monomers[idx];
    // Get neighbors of the current monomer
    std::vector<std::tuple<int, int, int>> neighbors = box.get_neighbors(current_monomer);

    // Check connections
    if (idx > 0) {
        const auto& prev_monomer = monomers[idx - 1];
        // If previous monomer is not a neighbor, return false
        if (std::find(neighbors.begin(), neighbors.end(), prev_monomer) == neighbors.end()) {
            return false;
        }
    }

    if (idx < monomers.size() - 1) {
        const auto& next_monomer = monomers[idx + 1];
        // If next monomer is not a neighbor, return false
        if (std::find(neighbors.begin(), neighbors.end(), next_monomer) == neighbors.end()) {
            return false;
        }
    }

    // All checks passed; monomer is connected properly
    return true;
}

bool RNA::would_be_connected_after_move(int idx, const std::tuple<int, int, int>& new_pos) const {
    int size = monomers.size();

    // Check previous monomer
    if (idx > 0) {
        int dist = chebyshev_distance(new_pos, monomers[idx - 1]);
        if (dist > 1) {
            return false;
        }
    }

    // Check next monomer
    if (idx < size - 1) {
        int dist = chebyshev_distance(new_pos, monomers[idx + 1]);
        if (dist > 1) {
            return false;
        }
    }

    return true;
}


/*
std::tuple<int, int, int> RNA::get_site_to_exchange(const BOX& box) const {
    // Initialize random number generator
    static std::mt19937 rng(std::random_device{}());

    // Step 1: Randomly select a monomer index idx
    std::uniform_int_distribution<int> dist(0, monomers.size() - 1);
    int idx = dist(rng);

    // Get current position of monomer idx
    auto current_pos = monomers[idx];

    // Get neighbors of current_pos
    std::vector<std::tuple<int, int, int>> neighbors = box.get_neighbors(current_pos);

    // Shuffle neighbors to randomize the order
    std::shuffle(neighbors.begin(), neighbors.end(), rng);

    // For each neighbor, check if move is valid
    for (const auto& neighbor_pos : neighbors) {
        Object* neighbor_obj = box.get_lattice(neighbor_pos);

        if (neighbor_obj->isempty()) {
            // Empty site; check if moving keeps RNA connected
            if (isconnected(idx, neighbor_pos)) {
                // Return the site and monomer index as needed
                return neighbor_pos;
            }
        } else if (neighbor_obj->Index() == 2) {
            // Neighbor is another RNA monomer; check swapping
            RNA* neighbor_rna = dynamic_cast<RNA*>(neighbor_obj);

            if (neighbor_rna == nullptr) {
                continue; // Skip if not an RNA object
            }

            // Get index of neighbor monomer in its RNA
            int neighbor_idx = neighbor_rna->get_monomer_index(neighbor_pos);

            if (neighbor_idx == -1) {
                continue; // Monomer not found in its RNA
            }

            // Simulate swapping positions
            // Check if both RNAs remain connected after swap
            bool this_rna_connected = true;
            bool neighbor_rna_connected = true;

            // For this RNA
            if (idx > 0) {
                this_rna_connected &= (chebyshev_distance(neighbor_pos, monomers[idx - 1]) <= 1);
            }
            if (idx < monomers.size() - 1) {
                this_rna_connected &= (chebyshev_distance(neighbor_pos, monomers[idx + 1]) <= 1);
            }

            // For neighbor RNA
            if (neighbor_idx > 0) {
                neighbor_rna_connected &= (chebyshev_distance(current_pos, neighbor_rna->monomers[neighbor_idx - 1]) <= 1);
            }
            if (neighbor_idx < neighbor_rna->monomers.size() - 1) {
                neighbor_rna_connected &= (chebyshev_distance(current_pos, neighbor_rna->monomers[neighbor_idx + 1]) <= 1);
            }

            if (this_rna_connected && neighbor_rna_connected) {
                // Valid swap
                return neighbor_pos;
            }
        } else {
            // Neighboring site contains another object; skip
            continue;
        }
    }

    // No valid site found
    throw std::runtime_error("No valid site found to exchange");
}

*/

int RNA::get_monomer_index(const std::tuple<int, int, int>& position) const {
    auto it = std::find(monomers.begin(), monomers.end(), position);
    if (it != monomers.end()) {
        return static_cast<int>(std::distance(monomers.begin(), it));
    } else {
        return -1; // Not found
    }
}

float RNA::compute_local_energy(const BOX& box) const {
    // Implement the local energy calculation for RNA objects
    auto neighbors = box.get_neighbors(position);
    float local_energy = 0.0f;    
    for (const auto& nxyz : neighbors) {
        Object* neighbor_obj = box.get_lattice(nxyz);
        local_energy -= box.E[Index()][neighbor_obj->Index()];
    }

    return local_energy;
}

std::vector<std::tuple<int, int, int>> RNA::get_positions() const {
    return monomers;
}

std::vector<std::tuple<int, int, int>> RNA::get_positions() const {
    return monomers;
}


// DHH1 Class Implementation
DHH1::DHH1(const std::tuple<int, int, int>& position_) : Object(position_) {}
DHH1::~DHH1() {}
int DHH1::Index() const { return 1; }

std::tuple<int, int, int> DHH1::get_site_to_exchange(const BOX& box) const {
    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, box.size - 1);

    int x, y, z;
    do {
        x = dist(rng);
        y = dist(rng);
        z = dist(rng);
    } while (std::make_tuple(x, y, z) == position);

    return std::make_tuple(x, y, z);
}
float DHH1::compute_local_energy(const BOX& box) const {
    // Implement the local energy calculation for DHH1 objects
    auto neighbors = box.get_neighbors(position);
    float local_energy = 0.0f;
    int neigh_count(0);
    for (const auto& nxyz : neighbors) {
        Object* neighbor_obj = box.get_lattice(nxyz);
        local_energy -= box.E[Index()][neighbor_obj->Index()];
        if(!neighbor_obj->isempty()){
            neigh_count+=1;
        }        
    }
    if(neigh_count>0){local_energy-=box.Evalence;}

    return local_energy;
}


double distance(std::tuple<int,int,int> site1, std::tuple<int,int,int> site2){
    return sqrt(pow(std::get<0>(site1)-std::get<0>(site2),2)+
    pow(std::get<1>(site1)-std::get<1>(site2),2)+
    pow(std::get<2>(site1)-std::get<2>(site2),2));
}

int chebyshev_distance(const std::tuple<int, int, int>& pos1, const std::tuple<int, int, int>& pos2) {
    return std::max({std::abs(std::get<0>(pos1) - std::get<0>(pos2)),
                     std::abs(std::get<1>(pos1) - std::get<1>(pos2)),
                     std::abs(std::get<2>(pos1) - std::get<2>(pos2))});
}