// Objects.cpp
#include "Utilities.h"
#include "Objects.h"
#include "BOX.h"
#include "ComputeLocalEnergy.h"
#include <random>
#include <stdexcept>
#include <vector>
#include <set>
#include <iostream>


// Object Class Implementation
Object::Object(int position_) : position(position_) {moved=false;}
Object::~Object() {}

int Object::getPosition() const { return position; }
std::vector<int> Object::get_positions() const {
    return {getPosition()};
}
void Object::setPosition(int new_pos) { position = new_pos; }
bool Object::isempty() const { return false; }
int Object::Index() const { return 0; }
int Object::get_monomer_index() const{return 0;};
void Object::print_position(int size) const{print_xyz(position,size);std::cout<<"\n";}

bool Object::isconnected(const BOX& box) const {return true;}

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

// RNA Class Implementation
RNA::RNA(std::shared_ptr<std::vector<int>> monomers_, int index_)
    : Object((*monomers_)[index_]), monomers(monomers_), index(index_) {}

RNA::~RNA() {}

int RNA::Index() const { return 2; }

int RNA::get_monomer_index() const {
        return index;
}

/*bool RNA::isconnected(int idx, const BOX& box) const {
    // Ensure idx is within bounds
    if (idx < 0 || idx >= monomers.size()) {
        std::cout<<idx<<std::endl;
        throw std::out_of_range("Index out of bounds in RNA::isconnected");
    }
    int current_monomer = monomers[idx];
    // Get neighbors of the current monomer
    std::vector<int> neighbors = box.get_neighbors(current_monomer);

    // Check connections
    if (idx > 0) {        
        if(chebyshev_distance(monomers[idx],monomers[idx-1],box.size)>1){return false;}
    }

    if (idx < monomers.size() - 1) {
        if(chebyshev_distance(monomers[idx],monomers[idx+1],box.size)>1){return false;}
    }
    // All checks passed; monomer is connected properly
    return true;
}*/

bool RNA::isconnected(const BOX& box) const {
    // Check connectivity with adjacent monomers
    // Check previous monomer
    if (index > 0) {
        if (chebyshev_distance(monomers->at(index), monomers->at(index-1), box.size) > 1) {
            return false;
        }
    }
    // Check next monomer
    if (index < monomers->size() - 1) {
        if (chebyshev_distance(monomers->at(index), monomers->at(index+1), box.size) > 1) {
            return false;
        }
    }

    return true;
}

std::vector<int> RNA::get_positions() const {
    return (*monomers);
}

void RNA::print_position(int size)const{
    std::cout<<" - ";
    for(auto& monomer : (*monomers)){        
        print_xyz(monomer,size);
        std::cout<<" - ";
    }std::cout<<"\n";
}

void RNA::setPosition(int new_pos){
    //int idx1 = get_monomer_index(old_pos);
    position=new_pos;
    if (index != -1) {
        monomers->at(index) = new_pos;
    }
    else{
        std::cout<<"cannot move the monomer"<<"\n";
    }
}

void RNA::get_swap_site_candidate(std::vector<std::shared_ptr<Object>>& object1,
                                std::vector<std::shared_ptr<Object>>& object2,
                                std::vector<int>& sites1,
                                std::vector<int>& sites2,
                                 const BOX& box, std::mt19937& rng) {
    // Decide whether to move a single monomer or the whole RNA
    std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
    double move_type_prob = prob_dist(rng);

    if (move_type_prob < 0.) {
        // Single monomer move (50% probability)
        // single monomer move
        // Get neighbors of the site
        object1.push_back(shared_from_this());
        
        sites1.push_back(monomers->at(index));

        std::vector<int> neighbors = box.get_neighbors(monomers->at(index));
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
    else{
        // select a direction by selecting a random neighbor index
        std::uniform_int_distribution<int> prob_dist(0, 25);
        int idx = prob_dist(rng);
        for(int i =0; i<monomers->size();i++){
            object1.push_back(box.get_lattice(monomers->at(i)));
            sites1.push_back(monomers->at(i));
            //for(auto& it: box.get_neighbors(monomers->at(i))){
            //    int x,y,z;
            //    std::tie(x,y,z)= to_xyz(it,box.size);
            //    std::cout<<x<<" "<<y<<" "<<z<<"\n";
            //}            
            int site2(box.get_neighbors(monomers->at(i))[idx]);
            sites2.push_back(site2);
            object2.push_back(box.get_lattice(site2));}
            //std::cout<<idx<<" \n";
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


