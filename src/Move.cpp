#include "Move.h"

Move::Move(Object* obj1, const std::tuple<int, int, int>& pos1,
           Object* obj2, const std::tuple<int, int, int>& pos2,
           MoveType type)
    : object1(obj1), object2(obj2), site1(pos1), site2(pos2), move_type(type) {
    // Initialization if needed
}

bool Move::validate(const BOX& box) const {
    // Perform connectivity checks and any other validations
    // For RNA objects, check if they remain connected after the move

    // Check for object1
    if (RNA* rna1 = dynamic_cast<RNA*>(object1)) {
        int idx1 = rna1->get_monomer_index(site1);
        if (idx1 != -1) {
            if (!rna1->would_be_connected_after_move(idx1, site2)) {
                return false;
            }
        }
    }

    // Check for object2 (if not Empty)
    if (object2 && !object2->isempty()) {
        if (RNA* rna2 = dynamic_cast<RNA*>(object2)) {
            int idx2 = rna2->get_monomer_index(site2);
            if (idx2 != -1) {
                if (!rna2->would_be_connected_after_move(idx2, site1)) {
                    return false;
                }
            }
        }
    }

    // Additional validations can be added here

    return true; // Move is valid
}

void Move::apply(BOX& box) {
    // Swap the objects in the lattice
    box.set_lattice(site1, object2);
    box.set_lattice(site2, object1);

    // Update the positions of the objects
    if (RNA* rna1 = dynamic_cast<RNA*>(object1)) {
        // Update the monomer position in rna1's monomers vector
        int idx1 = rna1->get_monomer_index(site1);
        if (idx1 != -1) {
            rna1->monomers[idx1] = site2;
        }
    } else {
        // For non-RNA objects
        object1->setPosition(site2);
    }

    if (object2 && !object2->isempty()) {
        if (RNA* rna2 = dynamic_cast<RNA*>(object2)) {
            // Update the monomer position in rna2's monomers vector
            int idx2 = rna2->get_monomer_index(site2);
            if (idx2 != -1) {
                rna2->monomers[idx2] = site1;
            }
        } else {
            object2->setPosition(site1);
        }
    } else if (object2 && object2->isempty()) {
        // If object2 is Empty, update its position
        object2->setPosition(site1);
    }

    // Invalidate cluster data if necessary
    box.clusters_valid = false;
}


void Move::revert(BOX& box) {
    // Swap back the objects in the lattice
    box.set_lattice(site1, object1);
    box.set_lattice(site2, object2);

    // Revert the positions of the objects
    if (RNA* rna1 = dynamic_cast<RNA*>(object1)) {
        // Revert the monomer position in rna1's monomers vector
        int idx1 = rna1->get_monomer_index(site2);
        if (idx1 != -1) {
            rna1->monomers[idx1] = site1;
        }
    } else {
        // For non-RNA objects
        object1->setPosition(site1);
    }

    if (object2 && !object2->isempty()) {
        if (RNA* rna2 = dynamic_cast<RNA*>(object2)) {
            // Revert the monomer position in rna2's monomers vector
            int idx2 = rna2->get_monomer_index(site1);
            if (idx2 != -1) {
                rna2->monomers[idx2] = site2;
            }
        } else {
            object2->setPosition(site2);
        }
    } else if (object2 && object2->isempty()) {
        // If object2 is Empty, revert its position
        object2->setPosition(site2);
    }

    // Invalidate cluster data if necessary
    box.clusters_valid = false;
}
