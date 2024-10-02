// MC.cpp
#include "MC.h"
#include <random>
#include <algorithm>

MC::MC(int size, int nparticles_, int npolymers_, int lpolymer_, const std::vector<std::vector<float>>& interactions, float temperature)
    : nparticles(nparticles_), npolymers(npolymers_), lpolymer(lpolymer_), E(interactions), T(temperature), box(size, nparticles_ + npolymers_, interactions) {
    generate_polymers(npolymers, lpolymer);
    generate_particles(nparticles);
}

void MC::generate_polymers(int npolymers, int lpolymer) {
    int npoly = 0;
    while (npoly < npolymers) {
        if (add_random_poly(lpolymer)) {
            ++npoly;
        }
    }
}

bool MC::add_random_poly(int lpolymer) {
    // Implement the function to add a random polymer of length lpolymer
    // This is non-trivial and depends on how you represent polymers
    // For now, return false as a placeholder
    return false;
}

void MC::generate_particles(int nparticles) {
    auto particles = generate_unique_triplets(nparticles, box.size - 1);
    for (int idx = 0; idx < particles.size(); ++idx) {
        int x = particles[idx][0];
        int y = particles[idx][1];
        int z = particles[idx][2];

        DHH1* new_dhh1 = new DHH1(std::make_tuple(x, y, z));
        box.set_lattice(std::make_tuple(x, y, z), new_dhh1);
        box.objects[npolymers + idx] = new_dhh1;  // Polymers must be generated first
    }
}

bool MC::monte_carlo_step() {
    // Perform a single Monte Carlo step using the Metropolis algorithm with site exchange

    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, box.objects.size() - 1);

    int idx = dist(rng);
    Object* object1 = box.objects[idx];

    DHH1* dhh1_object = dynamic_cast<DHH1*>(object1);
    if (!dhh1_object) {
        return false;  // Skip if not DHH1 object
    }

    std::tuple<int, int, int> site1 = object1->getPosition();
    std::tuple<int, int, int> site2 = dhh1_object->get_site_to_exchange(box);

    Object* object2 = box.get_lattice(site2);

    // Compute energy before the move
    float initial_energy = box.compute_local_energy(site1) + box.compute_local_energy(site2);

    // Swap the particles
    box.set_lattice(site1, object2);
    box.set_lattice(site2, object1);

    // Compute energy after the move
    float final_energy = box.compute_local_energy(site1) + box.compute_local_energy(site2);

    // Calculate energy difference
    float delta_e = final_energy - initial_energy;

    // Decide whether to accept the move
    std::uniform_real_distribution<float> rand_dist(0.0f, 1.0f);
    if (delta_e > 0 && rand_dist(rng) >= std::exp(-delta_e / T)) {
        // Reject the move (revert)
        box.set_lattice(site1, object1);
        box.set_lattice(site2, object2);
        return false;  // Move was rejected
    }
    return true;  // Move was accepted
}

std::vector<bool> MC::monte_carlo_steps(int steps) {
    std::vector<bool> success(steps);
    for (int step = 0; step < steps; ++step) {
        success[step] = monte_carlo_step();
    }
    return success;
}