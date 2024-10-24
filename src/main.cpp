#include "MC.h"
#include <iostream>
#include <algorithm>

int main() {
    int size = 32;
    int nparticles = 0;
    int npolymers = 1;
    int lpolymer = 2;
    float temperature = 0.1f;
    double Evalence=0;

    // Define interaction matrix E (3x3 matrix for example)
    std::vector<std::vector<float>> interactions = {
        {0.0f, 0.0f, 0.0f},
        {0.0f, 2.5f, 0.0f},
        {0.0f, 0.0f, 0.0f}
    };
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Define the range for the random integer (e.g., 1 to 100)
    std::uniform_int_distribution<> distr(1, 100);

    // Generate a random integer
    int randomInt = distr(gen);
    std::cout<<randomInt<<"\n";
    MC simulation(size, nparticles, npolymers, lpolymer, interactions,Evalence, temperature,13);
    for (int i=0;i<1;i++){
    std::cout<<i<<std::endl;
    simulation.monte_carlo_steps(pow(10,0));
    }
    std::cout<<simulation.get_RNA_positions().size()<<"\n";
    for(auto& it1 : simulation.get_RNA_positions()){
        for(auto& it: it1){
            std::cout<<it<<"-";
        }
        std::cout<<std::endl;
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
