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

// RNA Class Implementation
RNA::RNA(int length, const std::vector<int>& x, const std::vector<int>& y, const std::vector<int>& z)
    : Object(std::make_tuple(0, 0, 0)) {
    // Implement RNA initialization
}
RNA::~RNA() {}
int RNA::Index() const { return 2; }


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
