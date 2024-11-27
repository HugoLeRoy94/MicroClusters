// MC.cpp
#include "MC.h"
#include "Objects.h"
#include <random>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <set>
#include "ComputeLocalEnergy.h"

MC::MC(int size, int nparticles_, int npolymers_, int lpolymer_, const std::vector<std::vector<float>>& interactions,double Evalence_, float temperature, int seed, double diff_moves_ratio)
    : nparticles(nparticles_), npolymers(npolymers_), lpolymer(lpolymer_), E(interactions), T(temperature), box(size, nparticles_ + npolymers_, interactions,Evalence_,rng) {        
    RNA::diff_moves_ratio = diff_moves_ratio;
    rng.seed(seed);
    generate_polymers(npolymers, lpolymer);
    std::cout<<"polymer successfully generated"<<std::endl;
    generate_particles(nparticles);
    std::cout<<"particles successfully added"<<std::endl;

    energy = compute_energy();
    // Adjust the weight to make the probability of picking a polymer proportional to its length    
    //weights.resize(box.objects.size());
    //for (size_t i = 0; i < box.objects.size(); ++i) {
    //    if (i < static_cast<size_t>(npolymers)) {
    //        weights[i] = lpolymer;
    //    } else {
    //        weights[i] = 1;
    //    }
    //}
    //dist = std::discrete_distribution<int>(weights.begin(), weights.end());
    dist = std::uniform_int_distribution<int>(0,box.objects.size()-1);
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
    std::set<std::vector<int>> all_positions;
    for (const auto& obj : box.objects) {
        if (obj->Index() == 2) { // RNA objects have Index() == 2
            all_positions.insert(obj->get_positions());
        }
    }
    std::vector<std::vector<int>> all_positions_vec(all_positions.begin(),all_positions.end());
    return all_positions_vec;
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
        //std::cout<<counter<<" "<<pow(box.size, 3)<<std::endl;
        counter++;
        int idx = dist(rng);
        auto object1 = box.objects[idx];        
        std::vector<int> positions = object1->get_positions();

        // Randomly select one of the positions
        //std::uniform_int_distribution<int> pos_dist(0, positions.size() - 1);
        //int pos_idx = pos_dist(rng);
        //int site1 = positions[pos_idx];

        std::vector<std::shared_ptr<Object>> objects1;
        std::vector<std::shared_ptr<Object>> objects2;
        std::vector<int> sites1;        
        std::vector<int> sites2;        

        object1->get_swap_site_candidate(objects1,objects2, sites1,sites2, box, rng);
        if (sites2.back() == -1) {
            //for(auto it : sites2){
            //    std::cout<<it<<std::endl;
            //}
            //for(auto it : sites1){
            //    std::cout<<it<<std::endl;
            //}
            continue; // No valid candidate found, try again
        }

            // Create a move proposal
            Move move(objects1, objects2,sites1, sites2);
            // Compute energy before the move
            float initial_energy(0);
            for( int i = 0; i<sites1.size();i++){
                initial_energy += box.compute_local_energy(sites1[i]) + box.compute_local_energy(sites2[i]);
            }
            //std::cout<<"before \n object 1:\n";
            //for(auto& object: objects1){object->print_position(box.size);}
            //std::cout<<"object 2:\n";
            //for(auto& object: objects2){object->print_position(box.size);}
            //std::cout<<"\n";
            // Apply the move
            move.apply(box);
            //std::cout<<"after \n object 1:\n";
            //for(auto& object: objects1){object->print_position(box.size);}
            //std::cout<<"object 2:\n";
            //for(auto& object//: objects2){object->print_position(box.size);}
            //std::cout<<"\n";

                        // Validate the move
            if (!move.validate(box)) {
                //std::cout<<"move invalidate------------------------------------------"<<std::endl;
                move.revert(box);
                continue;
            }
            //else{
            //    std::cout<<"move validated"<<std::endl;
            //}   

            // Compute energy after the move
            float final_energy(0);
            for( int i = 0; i<sites1.size();i++){
                final_energy += box.compute_local_energy(sites1[i]) + box.compute_local_energy(sites2[i]);
            }

            // Calculate energy difference
            float delta_e = final_energy - initial_energy;

            //box.check_consistencty();

            // Decide whether to accept the move
            std::uniform_real_distribution<float> rand_dist(0.0f, 1.0f);
            if (delta_e > 0 && rand_dist(rng) >= std::exp(-delta_e / T)) {
                // Reject the move (revert)
                move.revert(box);
                return false;  // Move was rejected
            }
            energy += delta_e;
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

double MC::compute_energy()const{
    return box.total_energy();
}

double MC::get_energy()const{
    return energy;
}