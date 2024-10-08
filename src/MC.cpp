// MC.cpp
#include "MC.h"
#include <random>
#include <algorithm>
#include <stdexcept>
#include <iostream>

MC::MC(int size, int nparticles_, int npolymers_, int lpolymer_, const std::vector<std::vector<float>>& interactions,double Evalence_, float temperature)
    : nparticles(nparticles_), npolymers(npolymers_), lpolymer(lpolymer_), E(interactions), T(temperature), box(size, nparticles_ + npolymers_, interactions,Evalence_) {
    generate_polymers(npolymers, lpolymer);
    std::cout<<"polymer successfully generated"<<std::endl;
    generate_particles(nparticles);
    std::cout<<"particles successfully added"<<std::endl;
}

void MC::generate_polymers(int npolymers, int lpolymer) {
    int npoly = 0;
    int counter(0);    
    while (npoly < npolymers) {
        counter++;
        try{
            box.add_RNA(lpolymer);
            ++npoly;
        }
        catch(std::exception& e){}
        if(counter>pow(box.size,3)){
            throw std::runtime_error("cannot find a place to place a polymer");
        }
    }
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

        box.create_new_DHH1(std::make_tuple(x, y, z));
        //DHH1* new_dhh1 = new DHH1(std::make_tuple(x, y, z));
        //box.set_lattice(std::make_tuple(x, y, z), new_dhh1);
        //box.objects[npolymers + idx] = new_dhh1;  // Polymers must be generated first        
    }
}

bool MC::monte_carlo_step() {
    static std::mt19937 rng(std::random_device{}());

    // Adjust the weight to make the probability of picking a polymer proportional to its length
    std::vector<int> weights(box.objects.size());
    for (int i = 0; i < box.objects.size(); ++i) {
        if (i < npolymers) {
            weights[i] = lpolymer;
        } else {
            weights[i] = 1;
        }
    }

    std::discrete_distribution<int> dist(weights.begin(), weights.end());

    int counter = 0;
    while (counter < pow(box.size, 3)) {
        counter++;
        int idx = dist(rng);
        Object* object1 = box.objects[idx];

        // Get positions associated with object1
        std::vector<std::tuple<int, int, int>> positions = object1->get_positions();

        // Randomly select one of the positions
        std::uniform_int_distribution<int> pos_dist(0, positions.size() - 1);
        int pos_idx = pos_dist(rng);
        auto site1 = positions[pos_idx];

        // Get neighbors of site1
        std::vector<std::tuple<int, int, int>> neighbors = box.get_neighbors(site1);
        std::shuffle(neighbors.begin(), neighbors.end(), rng);

        for (const auto& site2 : neighbors) {
            Object* object2 = box.get_lattice(site2);

            // Create a move proposal
            Move move(object1, site1, object2, site2, MoveType::Swap);

            // Validate the move
            if (!move.validate(box)) {
                continue;
            }

            // Compute energy before the move
            float initial_energy = box.compute_local_energy(site1) + box.compute_local_energy(site2);

            // Apply the move
            move.apply(box);

            // Compute energy after the move
            float final_energy = box.compute_local_energy(site1) + box.compute_local_energy(site2);

            // Calculate energy difference
            float delta_e = final_energy - initial_energy;

            // Decide whether to accept the move
            std::uniform_real_distribution<float> rand_dist(0.0f, 1.0f);
            if (delta_e > 0 && rand_dist(rng) >= std::exp(-delta_e / T)) {
                // Reject the move (revert)
                move.revert(box);
                return false;  // Move was rejected
            }
            return true;  // Move was accepted
        }
    }
//    return false; // No valid move found
    throw std::out_of_range("No valid move found");
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