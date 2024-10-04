// MC.cpp
#include "MC.h"
#include <random>
#include <algorithm>

MC::MC(int size, int nparticles_, int npolymers_, int lpolymer_, const std::vector<std::vector<float>>& interactions,double Evalence_, float temperature)
    : nparticles(nparticles_), npolymers(npolymers_), lpolymer(lpolymer_), E(interactions), T(temperature), box(size, nparticles_ + npolymers_, interactions,Evalence_) {
    generate_polymers(npolymers, lpolymer);
    generate_particles(nparticles);
}

void MC::generate_polymers(int npolymers, int lpolymer) {
    //int npoly = 0;
    //while (npoly < npolymers) {
    //    if (add_random_poly(lpolymer)) {
    //        ++npoly;
    //    }
    //}
    return;
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

        box.create_new_DHH1(std::make_tuple(x, y, z),npolymers+idx);
        //DHH1* new_dhh1 = new DHH1(std::make_tuple(x, y, z));
        //box.set_lattice(std::make_tuple(x, y, z), new_dhh1);
        //box.objects[npolymers + idx] = new_dhh1;  // Polymers must be generated first        
    }
}

bool MC::monte_carlo_step() {
    // Perform a single Monte Carlo step using the Metropolis algorithm with site exchange

    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, box.objects.size() - 1);

    int idx = dist(rng);
    Object* object1 = box.objects[idx];    

    std::tuple<int, int, int> site1 = object1->getPosition();
    std::tuple<int, int, int> site2 = object1->get_site_to_exchange(box);
    Object* object2 = box.get_lattice(site2);
    //Object* object2 = object1->get_object_to_exchange(box);

    // Compute energy before the move
    float initial_energy = box.compute_local_energy(site1) + box.compute_local_energy(site2);

    // Swap the particles
    box.swap(site1, site2);

    // Compute energy after the move
    float final_energy = box.compute_local_energy(site1) + box.compute_local_energy(site2);

    // Calculate energy difference
    float delta_e = final_energy - initial_energy;

    // Decide whether to accept the move
    std::uniform_real_distribution<float> rand_dist(0.0f, 1.0f);
    if (delta_e > 0 && rand_dist(rng) >= std::exp(-delta_e / T)) {
        // Reject the move (revert)
        box.swap(site1, site2);
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

double MC::get_energy()const{
    return box.total_energy();
}