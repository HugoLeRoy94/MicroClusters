// Objects.cpp
#include "Objects.h"
#include "BOX.h"
#include <random>

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
