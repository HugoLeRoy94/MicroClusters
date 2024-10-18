#include "Move.h"
#include <stdexcept>

Move::Move(const std::vector<std::shared_ptr<Object>>& objs1,
            const std::vector<std::shared_ptr<Object>>& objs2,
            const std::vector<int>& old_pos,
            const std::vector<int>& new_pos)
    : objects1(objs1), objects2(objs2), sites1(old_pos), sites2(new_pos) {
    // Initialization if needed
    if(objs1.size() != old_pos.size() or objs1.size() !=old_pos.size()){
        throw std::length_error("invalid set of new position for the size of the object");
    }

}

bool Move::validate(const BOX& box) const {
    // Perform connectivity checks and any other validations
    // For RNA objects, check if they remain connected after the move
            for(int i = 0; i < objects1.size();i++){
                if(!objects1[i]->would_be_connected_after_move(
                    objects1[i]->get_monomer_index(sites1[i]),
                    sites2[i],
                    box.size) &&
                    !objects2[i]->would_be_connected_after_move(
                    objects2[i]->get_monomer_index(sites2[i]),
                    sites1[i],
                    box.size)
                    ){
                        return false;
                    }
            }
    // Additional validations can be added here

    return true; // Move is valid
    
}

void Move::apply(BOX& box) {
            for(int i =0; i < objects1.size();i++){
                swap_sites(sites1[i],sites2[i],objects1[i],objects2[i],box);
            }        
    // Invalidate cluster data if necessary
    box.clusters_valid = false;
}


void Move::revert(BOX& box){        
            for(int i = 0; i<objects1.size();i++){
                swap_sites(sites2[i],sites1[i],objects1[i],objects2[i],box);
            }
    // Invalidate cluster data if necessary
    box.clusters_valid = false;
}


void Move::swap_sites(int site1, int site2,std::shared_ptr<Object> obj1, std::shared_ptr<Object> obj2, BOX& box) {
    // Swap objects in the lattice
    box.set_lattice(site1, obj2);
    box.set_lattice(site2, obj1);
    // Update positions in the objects
    obj1->setPosition(site1,site2); // prev site, next site
    obj2->setPosition(site2,site1);
}
