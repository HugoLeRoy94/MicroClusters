// BOX.h
#ifndef BOX_H
#define BOX_H

#include <vector>
#include <tuple>
#include <memory>
#include <random>

class Object;
class RNA;
//class DHH1;

class BOX {
public:
    int size;
    std::vector<std::shared_ptr<Object>> lattice;  // 1D vector to represent 3D lattice
    std::vector<std::shared_ptr<Object>> objects;
    std::vector<std::vector<float>> E;  // Interaction matrix
    double Evalence;
    bool clusters_valid = false;

    BOX(int size_, int nobjects,const std::vector<std::vector<float>>& Interactions,double Evalence_,std::mt19937& rng_);

    void create_new_DHH1(int index);
    std::shared_ptr<RNA> add_RNA(int length);
    std::shared_ptr<Object> get_lattice(int index) const;
    void set_lattice(int index, std::shared_ptr<Object> obj);
    float compute_local_energy(int index) const;
    /*inline float compute_local_energy(int index) const {
        auto obj = get_lattice(index);
        if (obj->isempty()) {
            return 0.0f;
        }
        return obj->compute_local_energy(*this);
    }*/
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
    std::mt19937& rng;  // Reference to the RNG from MC
};
/*
inline int compute_n_bits(int size) {
    int n_bits = 0;
    int temp_size = size;
    while (temp_size >>= 1) ++n_bits;
    return n_bits;
}
inline int to_single_index(int x, int y, int z, int mask, int n_bits) {
    return ((x & mask) << (2 * n_bits)) | ((y & mask) << n_bits) | (z & mask);
}
inline bool is_power_of_two(int n) {
    return n > 0 && (n & (n - 1)) == 0;
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
}*/
#endif // BOX_H
