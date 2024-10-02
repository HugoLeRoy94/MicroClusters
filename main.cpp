// main.cpp
#include "MC.h"
#include <iostream>
#include <algorithm>

int main() {
    int size = 10;
    int nparticles = 100;
    int npolymers = 0;
    int lpolymer = 0;
    float temperature = 1.0f;

    // Define interaction matrix E (3x3 matrix for example)
    std::vector<std::vector<float>> interactions = {
        {0.0f, 1.0f, 2.0f},
        {1.0f, 0.0f, 1.5f},
        {2.0f, 1.5f, 0.0f}
    };

    MC simulation(size, nparticles, npolymers, lpolymer, interactions, temperature);

    int steps = 1000;
    auto success = simulation.monte_carlo_steps(steps);

    // Output the results
    int accepted_moves = std::count(success.begin(), success.end(), true);
    std::cout << "Accepted moves: " << accepted_moves << "/" << steps << std::endl;

    return 0;
}