// MC.cpp
#include "MC.h"
#include <random>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include "ComputeLocalEnergy.h"

MC::MC(int size, int nparticles_, int npolymers_, int lpolymer_, const std::vector<std::vector<float>>& interactions,double Evalence_, float temperature, int seed)
    : nparticles(nparticles_), npolymers(npolymers_), lpolymer(lpolymer_), E(interactions), T(temperature), box(size, nparticles_ + npolymers_, interactions,Evalence_,rng) {        
    rng.seed(seed);
    generate_polymers(npolymers, lpolymer);
    std::cout<<"polymer successfully generated"<<std::endl;
    generate_particles(nparticles);
    std::cout<<"particles successfully added"<<std::endl;

    // Adjust the weight to make the probability of picking a polymer proportional to its length    
    weights.resize(box.objects.size());
    for (size_t i = 0; i < box.objects.size(); ++i) {
        if (i < static_cast<size_t>(npolymers)) {
            weights[i] = lpolymer;
        } else {
            weights[i] = 1;
        }
    }
    dist = std::discrete_distribution<int>(weights.begin(), weights.end());
}
// Function to get positions of all DHH1 particles
std::vector<int> MC::get_DHH1_positions() const {
    std::vector<int> positions;
    for (const auto& obj : box.objects) {
        if (obj->Index() == 1) { // DHH1 objects have Index() == 1
            positions.push_back(obj->getPosition());
        }
    }
    return positions;
}

// Function to get positions of all RNA polymers
std::vector<std::vector<int>> MC::get_RNA_positions() const {
    std::vector<std::vector<int>> all_positions;
    for (const auto& obj : box.objects) {
        if (obj->Index() == 2) { // RNA objects have Index() == 2
            all_positions.push_back(obj->get_positions());
        }
    }
    return all_positions;
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

void MC::generate_particles(int nparticles) {
    std::vector<int> empty_indices;
    int total_sites = box.size * box.size * box.size;
    
    // Collect all empty lattice indices
    for (int idx = 0; idx < total_sites; ++idx) {
        if (box.get_lattice(idx)->isempty()) {
            empty_indices.push_back(idx);
        }
    }

    // Check if there are enough empty sites
    if (nparticles > empty_indices.size()) {
        throw std::runtime_error("Not enough empty positions to place all DHH1 particles.");
    }

    // Shuffle the empty indices
    //static std::mt19937 rng(std::random_device{}());
    std::shuffle(empty_indices.begin(), empty_indices.end(), rng);

    // Place DHH1 particles in the first nparticles empty indices
    for (int i = 0; i < nparticles; ++i) {
        int idx = empty_indices[i];
        box.create_new_DHH1(idx);
    }
}

bool MC::monte_carlo_step() {
    //static std::mt19937 rng(std::random_device{}());

    

    int counter = 0;
    while (counter < pow(box.size, 3)){
        counter++;
        int idx = dist(rng);
        auto object1 = box.objects[idx];
        // Get positions associated with object1
        std::vector<int> positions = object1->get_positions();

        // Randomly select one of the positions
        std::uniform_int_distribution<int> pos_dist(0, positions.size() - 1);
        int pos_idx = pos_dist(rng);
        int site1 = positions[pos_idx];

        //while(true){
        // Get neighbors of site1
        // Get a random swap site candidate from object1
        //# need to find a way to return -1 when the RNA has tried all of tis neighbors.
        int site2 = object1->get_swap_site_candidate(site1, box,rng);
        if (site2 == -1) {
            continue; // No valid candidate found, try again
        }
            // start from a random point, and go back the the begining of candidates if the end is crossed
            auto object2 = box.get_lattice(site2);

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
        //}
    }
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
