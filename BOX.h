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

    BOX(int size_, int nobjects,const std::vector<std::vector<float>>& Interactions);
    ~BOX();

    void create_new_DHH1(const std::tuple<int, int, int>& site, int object_idx);
    Object* get_lattice(const std::tuple<int, int, int>& site) const;
    void swap(const std::tuple<int, int, int>& site1, const std::tuple<int, int, int>& site2);
    float compute_local_energy(const std::tuple<int, int, int>& xyz) const;
    float total_energy() const;
    std::vector<std::tuple<int, int, int>> get_neighbors(const std::tuple<int, int, int>& xyz) const;
    bool has_free_neighbor(const std::tuple<int, int, int>& xyz) const;

    // Additional methods
    // Function to build clusters and update member variables
    void build_clusters();

    // Getters for cluster data
    const std::vector<int>& get_cluster_indices();
    const std::vector<int>& get_cluster_starts();

    size_t get_cluster_indices_size();
    size_t get_cluster_starts_size();

    std::vector<int> cluster_size();
    double average_cluster_size();
    double compute_av_Nneigh() const;
private:
    int npolymers;
    // Disallow copying
    BOX(const BOX&) = delete;
    BOX& operator=(const BOX&) = delete;

    // Cluster data
    std::vector<int> cluster_indices;
    std::vector<int> cluster_starts;
    bool clusters_valid = false;
};

int to_single_index(int x, int y, int z, int L);
std::tuple<int, int, int> to_xyz(int index, int L);
std::vector<std::array<int, 3>> generate_unique_triplets(int N, int L);

#endif // BOX_H
