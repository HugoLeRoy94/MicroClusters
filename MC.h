// MC.h
#ifndef MC_H
#define MC_H

#include "BOX.h"
#include "Objects.h"
#include <vector>

class MC {
public:
    int npolymers;
    int lpolymer;
    int nparticles;
    float T;  // Temperature
    std::vector<std::vector<float>> E;  // Interaction matrix

    BOX box;

    MC(int size, int nparticles_, int npolymers_, int lpolymer_, const std::vector<std::vector<float>>& interactions,double Evalence_, float temperature);

    void generate_polymers(int npolymers, int lpolymer);
    bool add_random_poly(int lpolymer);
    void generate_particles(int nparticles);
    bool monte_carlo_step();
    std::vector<bool> monte_carlo_steps(int steps);
};

#endif // MC_H
