// BOX.h
#ifndef BOX_H
#define BOX_H

#include <vector>
#include <tuple>
#include "Objects.h"

class BOX {
public:
    int size;
    std::vector<Object*> lattice;  // 1D vector to represent 3D lattice
    std::vector<Object*> objects;
    std::vector<std::vector<float>> E;  // Interaction matrix

    BOX(int size_, int nobjects, const std::vector<std::vector<float>>& Interactions);
    ~BOX();

    void set_lattice(const std::tuple<int, int, int>& site, Object* object);
    Object* get_lattice(const std::tuple<int, int, int>& site) const;
    float compute_local_energy(const std::tuple<int, int, int>& xyz) const;
    float total_energy() const;
    std::vector<std::tuple<int, int, int>> get_neighbors(const std::tuple<int, int, int>& xyz) const;
    bool has_free_neighbor(const std::tuple<int, int, int>& xyz) const;

    // Additional methods
    std::tuple<std::vector<int>, std::vector<int>> build_clusters();
    std::vector<int> cluster_size();
    int compute_av_Nneigh() const;
};

int to_single_index(int x, int y, int z, int L);
std::tuple<int, int, int> to_xyz(int index, int L);

std::vector<std::array<int, 3>> generate_unique_triplets(int N, int L);

#endif // BOX_H
