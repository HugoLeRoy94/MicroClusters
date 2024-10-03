// Objects.cpp
#include "Objects.h"
#include "BOX.h"
#include <random>

// Object Class Implementation
Object::Object(const std::tuple<int, int, int>& position_) : position(position_) {}
Object::~Object() {}

const std::tuple<int, int, int>& Object::getPosition() const { return position; }
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
RNA::RNA(int length, const std::vector<int>& x, const std::vector<int>& y, const std::vector<int>& z)
    : Object(std::make_tuple(0, 0, 0)) {
    // Implement RNA initialization
}
RNA::~RNA() {}
int RNA::Index() const { return 2; }
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
