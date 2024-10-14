#include "MC.h"
#include <iostream>
#include <algorithm>

int main() {
    int size = 128;
    int nparticles = 210;
    int npolymers = 0;
    int lpolymer = 0;
    float temperature = 0.1f;
    double Evalence=0;

    // Define interaction matrix E (3x3 matrix for example)
    std::vector<std::vector<float>> interactions = {
        {0.0f, 0.0f, 0.0f},
        {0.0f, 2.5f, 0.0f},
        {0.0f, 0.0f, 0.0f}
    };

    MC simulation(size, nparticles, npolymers, lpolymer, interactions,Evalence, temperature,651);
    for (int i=0;i<10;i++){
    std::cout<<i<<std::endl;
    simulation.monte_carlo_steps(pow(10,6));

    }
    // Compute cluster sizes
    //std::vector<int> cluster_sizes = simulation.box.cluster_size();
    //std::cout << "Cluster sizes:\n";
    //for (const auto& size : cluster_sizes) {
    //    std::cout << size << " ";
    //}
    //std::cout << std::endl;
//
    //// Compute average number of occupied neighbors
    //double avg_neighbors = simulation.box.compute_av_Nneigh();
    //std::cout << "Average number of occupied neighbors per occupied site: " << avg_neighbors << std::endl;

    return 0;
}
