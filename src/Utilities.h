// Utilities.h
#ifndef UTILITIES_H
#define UTILITIES_H

#include <tuple>
#include <stdexcept>
#include <algorithm>


// Utility functions

inline bool is_power_of_two(int n) {
    return n > 0 && (n & (n - 1)) == 0;
}

inline int compute_n_bits(int size) {
    int n_bits = 0;
    int temp_size = size;
    while (temp_size >>= 1) ++n_bits;
    return n_bits;
}

inline int to_single_index(int x, int y, int z, int size) {
    if (!is_power_of_two(size)) {
        throw std::invalid_argument("Size must be a power of 2 for bitmasking optimization.");
    }
    int mask = size - 1;
    int n_bits = compute_n_bits(size);
    return ((x & mask) << (2 * n_bits)) | ((y & mask) << n_bits) | (z & mask);
}

inline std::tuple<int, int, int> to_xyz(int index, int size) {
    if (!is_power_of_two(size)) {
        throw std::invalid_argument("Size must be a power of 2 for bitmasking optimization.");
    }
    int mask = size - 1;
    int n_bits = compute_n_bits(size);
    int x = (index >> (2 * n_bits)) & mask;
    int y = (index >> n_bits) & mask;
    int z = index & mask;
    return std::make_tuple(x, y, z);
}

// Utilities.h (add at the end)

inline int chebyshev_distance(int pos1, int pos2, int size) {
    int x1, y1, z1;
    std::tie(x1, y1, z1) = to_xyz(pos1, size);
    int x2, y2, z2;
    std::tie(x2, y2, z2) = to_xyz(pos2, size);

    int dx = std::abs(x1 - x2);
    int dy = std::abs(y1 - y2);
    int dz = std::abs(z1 - z2);

    // Account for periodic boundary conditions
    dx = std::min(dx, size - dx);
    dy = std::min(dy, size - dy);
    dz = std::min(dz, size - dz);

    return std::max({dx, dy, dz});
}


#endif // UTILITIES_H
