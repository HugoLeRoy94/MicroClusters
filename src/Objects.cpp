// Objects.cpp
#include "Objects.h"
#include "BOX.h"
#include <random>
#include <stdexcept>


namespace{
    struct concrete_DHH1: public DHH1{
        concrete_DHH1(int position_): DHH1(position_){}
    };
    struct concrete_RNA: public RNA{
        concrete_RNA(std::vector<int>positions_):RNA(positions_){}
    };
    struct concrete_Empty: public Empty{
        concrete_Empty(int position_): Empty(position_){}
    };
}


// Object Class Implementation
Object::Object(int position_) : position(position_) {}
Object::~Object() {}

int Object::getPosition() const { return position; }
std::vector<int> Object::get_positions() const {
    return {getPosition()};
}
void Object::setPosition(int value) { position = value; }
bool Object::isempty() const { return false; }
int Object::Index() const { return 0; }
int Object::get_site_to_exchange(const BOX& box) const { return position; }

// Empty Class Implementation
Empty::Empty(int position_) : Object(position_) {}

std::shared_ptr<Empty> Empty::make_shared_ptr(int position_){
    return std::make_shared<concrete_Empty>(position_);
}
Empty::~Empty() {}
bool Empty::isempty() const { return true; }
int Empty::Index() const { return 0; }
float Empty::compute_local_energy(const BOX& box) const {
    // Empty objects have zero energy
    return 0.0f;
}

// RNA Class Implementation
RNA::RNA(std::vector<int>& monomers_)
    : Object(monomers_[0]) {
        monomers = monomers_;
}

std::shared_ptr<RNA> RNA::make_shared_ptr(std::vector<int>& monomers_){
    return std::make_shared<concrete_RNA>(monomers_);
}

RNA::~RNA() {}

int RNA::Index() const { return 2; }

bool RNA::isconnected(int idx, const BOX& box) const {
    // Ensure idx is within bounds
    if (idx < 0 || idx >= monomers.size()) {
        throw std::out_of_range("Index out of bounds in RNA::isconnected");
    }

    int current_monomer = monomers[idx];
    // Get neighbors of the current monomer
    std::vector<int> neighbors = box.get_neighbors(current_monomer);

    // Check connections
    if (idx > 0) {
        int prev_monomer = monomers[idx - 1];
        // If previous monomer is not a neighbor, return false
        if (std::find(neighbors.begin(), neighbors.end(), prev_monomer) == neighbors.end()) {
            return false;
        }
    }

    if (idx < monomers.size() - 1) {
        int next_monomer = monomers[idx + 1];
        // If next monomer is not a neighbor, return false
        if (std::find(neighbors.begin(), neighbors.end(), next_monomer) == neighbors.end()) {
            return false;
        }
    }

    // All checks passed; monomer is connected properly
    return true;
}

bool RNA::would_be_connected_after_move(int idx, int new_pos, int L) const {
    int size = static_cast<int>(monomers.size());

    // Check previous monomer
    if (idx > 0) {
        int dist = chebyshev_distance(new_pos, monomers[idx - 1], L);
        if (dist > 1) {
            return false;
        }
    }

    // Check next monomer
    if (idx < size - 1) {
        int dist = chebyshev_distance(new_pos, monomers[idx + 1], L);
        if (dist > 1) {
            return false;
        }
    }

    return true;
}

int RNA::get_monomer_index(int position) const {
    auto it = std::find(monomers.begin(), monomers.end(), position);
    if (it != monomers.end()) {
        return static_cast<int>(std::distance(monomers.begin(), it));
    } else {
        return -1; // Not found
    }
}

float RNA::compute_local_energy(const BOX& box) const {
    float local_energy = 0.0f;
    for (const auto& idx : monomers) {
        auto neighbors = box.get_neighbors(idx);
        for (const auto& nidx : neighbors) {
            auto neighbor_obj = box.get_lattice(nidx);
            local_energy -= box.E[Index()][neighbor_obj->Index()];
        }
    }
    return local_energy;
}

std::vector<int> RNA::get_positions() const {
    return monomers;
}

// DHH1 Class Implementation
DHH1::DHH1(int position_) : Object(position_) {}

std::shared_ptr<DHH1> DHH1::make_shared_ptr(int position_){
    return std::make_shared<concrete_DHH1>(position_);
}
DHH1::~DHH1() {}
int DHH1::Index() const { return 1; }

int DHH1::get_site_to_exchange(const BOX& box) const {
    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, box.size * box.size * box.size - 1);

    int idx;
    do {
        idx = dist(rng);
    } while (idx == position);

    return idx;
}
float DHH1::compute_local_energy(const BOX& box) const {
    auto neighbors = box.get_neighbors(position);
    float local_energy = 0.0f;
    int neigh_count(0);
    for (const auto& nidx : neighbors) {
        auto neighbor_obj = box.get_lattice(nidx);
        local_energy -= box.E[Index()][neighbor_obj->Index()];
        if(!neighbor_obj->isempty()){
            neigh_count+=1;
        }        
    }
    if(neigh_count>0){local_energy-=box.Evalence;}

    return local_energy;
}

int chebyshev_distance(int pos1, int pos2, int L) {
    int x1, y1, z1;
    std::tie(x1, y1, z1) = to_xyz(pos1, L);
    int x2, y2, z2;
    std::tie(x2, y2, z2) = to_xyz(pos2, L);

    int dx = std::abs(x1 - x2);
    int dy = std::abs(y1 - y2);
    int dz = std::abs(z1 - z2);

    // Account for periodic boundary conditions
    dx = std::min(dx, L - dx);
    dy = std::min(dy, L - dy);
    dz = std::min(dz, L - dz);

    return std::max({dx, dy, dz});
}
