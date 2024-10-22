// MC.h
#ifndef MC_H
#define MC_H

#include "BOX.h"
#include "Objects.h"
#include "Move.h"
#include <vector>

class MC {
public:
    int npolymers;
    int lpolymer;
    int nparticles;
    float T;  // Temperature
    std::vector<std::vector<float>> E;  // Interaction matrix
    std::vector<int> get_DHH1_positions() const;
    std::vector<std::vector<int>> get_RNA_positions() const;

    BOX box;

    MC(int size, int nparticles_, int npolymers_, int lpolymer_, const std::vector<std::vector<float>>& interactions,double Evalence_, float temperature, int seed);

    void generate_polymers(int npolymers, int lpolymer);
    void generate_particles(int nparticles);
    bool monte_carlo_step();
    std::vector<bool> monte_carlo_steps(int steps);
    double get_energy()const;

    private:
        std::vector<int> weights;
        std::mt19937 rng;
        std::uniform_int_distribution<int> dist;
};

#endif // MC_H
