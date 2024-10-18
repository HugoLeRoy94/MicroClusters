// Objects.cpp
#include "Utilities.h"
#include "Objects.h"
#include "BOX.h"
#include "ComputeLocalEnergy.h"
#include <random>
#include <stdexcept>
#include <vector>


// Object Class Implementation
Object::Object(int position_) : position(position_) {}
Object::~Object() {}

int Object::getPosition() const { return position; }
std::vector<int> Object::get_positions() const {
    return {getPosition()};
}
void Object::setPosition(int old_pos, int new_pos) { position = new_pos; }
bool Object::isempty() const { return false; }
int Object::Index() const { return 0; }

// Empty Class Implementation
Empty::Empty(int position_) : Object(position_) {}

Empty::~Empty() {}
bool Empty::isempty() const { return true; }
int Empty::Index() const { return 0; }
void Empty::get_swap_site_candidate(std::vector<std::shared_ptr<Object>>& object1,
                                std::vector<std::shared_ptr<Object>>& object2,
                                std::vector<int>& sites1,
                                std::vector<int>& sites2,
                                 const BOX& box, std::mt19937& rng){return;}
bool Empty::would_be_connected_after_move(int idx, int new_pos, int L) const {return true;}
int Empty   ::get_monomer_index(int site) const{return 0;};
// RNA Class Implementation
RNA::RNA(std::vector<int>& monomers_)
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


std::vector<int> RNA::get_positions() const {
    return monomers;
}

void RNA::setPosition(int old_pos, int new_pos){
    int idx1 = get_monomer_index(old_pos);
            if (idx1 != -1) {
                monomers[idx1] = new_pos;
            }
}

void RNA::get_swap_site_candidate(std::vector<std::shared_ptr<Object>>& object1,
                                std::vector<std::shared_ptr<Object>>& object2,
                                std::vector<int>& sites1,
                                std::vector<int>& sites2,
                                 const BOX& box, std::mt19937& rng) {    
    // single monomer move
    // Get neighbors of the site
    object1.push_back(shared_from_this());
    // select a random site
    std::uniform_int_distribution<int> pos_dist(0, monomers.size() - 1);
    int pos_idx = pos_dist(rng);
    sites1.push_back(monomers[pos_idx]);

    std::vector<int> neighbors = box.get_neighbors(monomers[pos_idx]);    
    // Randomly shuffle the neighbors
    std::shuffle(neighbors.begin(), neighbors.end(), rng);
    // Iterate over neighbors to find a valid candidate
    for (const auto& neighbor_idx : neighbors) {
        auto neighbor_obj = box.get_lattice(neighbor_idx);
        if (neighbor_obj.get() != this) {
            sites2.push_back(neighbor_idx);
            object2.push_back(neighbor_obj);
        }
    }    
}

// DHH1 Class Implementation
DHH1::DHH1(int position_) : Object(position_){}

DHH1::~DHH1() {}
int DHH1::Index() const { return 1; }

void DHH1::get_swap_site_candidate(std::vector<std::shared_ptr<Object>>& object1,
                                    std::vector<std::shared_ptr<Object>>& object2,
                                    std::vector<int>& sites1,
                                    std::vector<int>& sites2, const BOX& box, std::mt19937& rng){
    object1.push_back(shared_from_this());
    sites1.push_back(position);
    // Initialize RNG
    //static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, box.lattice.size() - 1);
    int idx;
    do {
        idx = dist(rng);
    } while (idx == position);
    sites2.push_back(idx);
    object2.push_back(box.get_lattice(idx));
}


bool DHH1::would_be_connected_after_move(int idx, int new_pos, int L) const {return true;}

int DHH1::get_monomer_index(int site) const{return 0;};