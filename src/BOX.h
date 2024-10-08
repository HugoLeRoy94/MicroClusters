// BOX.h
#ifndef BOX_H
#define BOX_H

#include <vector>
#include <tuple>
#include <memory>
#include "Objects.h"

class BOX {
public:
    int size;
    std::vector<std::shared_ptr<Object>> lattice;  // 1D vector to represent 3D lattice
    std::vector<std::shared_ptr<Object>> objects;
    std::vector<std::vector<float>> E;  // Interaction matrix
    double Evalence;
    bool clusters_valid = false;

    BOX(int size_, int nobjects,const std::vector<std::vector<float>>& Interactions,double Evalence_);

    void create_new_DHH1(int index);
    std::shared_ptr<RNA> add_RNA(int length);
    std::shared_ptr<Object> get_lattice(int index) const;
    void set_lattice(int index, std::shared_ptr<Object> obj);
    void swap(int idx1, int idx2);
    float compute_local_energy(int index) const;
    float total_energy() const;
    std::vector<int> get_neighbors(int index) const;
    bool has_free_neighbor(int index) const;

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

    std::vector<int> generate_unique_indices(int N);

private:
    int npolymers;
    // Disallow copying
    BOX(const BOX&) = delete;
    BOX& operator=(const BOX&) = delete;

    // Cluster data
    std::vector<int> cluster_indices;
    std::vector<int> cluster_starts;    

    int random_free_site();
};

int to_single_index(int x, int y, int z, int L);
std::tuple<int, int, int> to_xyz(int index, int L);

#endif // BOX_H
