#include "Move.h"
#include <iostream>
#include <stdexcept>
#include <tuple>

Move::Move(const std::vector<std::shared_ptr<Object>>& objs1,
            const std::vector<std::shared_ptr<Object>>& objs2,
            const std::vector<int>& old_pos,
            const std::vector<int>& new_pos)
    : objects1(objs1), objects2(objs2), sites1(old_pos), sites2(new_pos) {
    // Initialization if needed
    if(objs1.size() != old_pos.size() or objs1.size() !=old_pos.size()){
        throw std::length_error("invalid set of new position for the size of the object");
    }
    for(int i = 0;i <objs2.size(); i++){
        rflow[new_pos[i]] = old_pos[i]-new_pos[i];
    }
}

bool Move::validate(const BOX& box) const {
    // Perform connectivity checks and any other validations
    // For RNA objects, check if they remain connected after the move
            for(int i = 0; i < objects1.size();i++){
                //int x,y,z;
                //std::tie(x, y, z) = to_xyz(sites1[i],box.size);
                //std::cout<<x<<" "<<y<<" "<<z<<"\n";
                //std::tie(x, y, z) = to_xyz(sites2[i],box.size);
                //std::cout<<x<<" "<<y<<" "<<z<<"\n"<<"\n";
                if(!objects1[i]->isconnected(box) ||!objects2[i]->isconnected(box)){
                        return false;
                    }
            }
    // Additional validations can be added here

    return true; // Move is valid
    
}

void Move::apply(BOX& box) {
    // first moves the objects 1
    for(int i =0; i < objects1.size();i++){
        move(objects1[i],sites2[i],box);
        objects1[i]->moved = true;
    }
    // moves the second object according to the flow matrix
    /*for(auto& object : objects2){
        if(object->moved){continue;}
        int new_pos(object->getPosition());
        try{
            while(true){
                new_pos += rflow.at(new_pos);
            }
        } catch(const std::out_of_range& e){}
        move(object,new_pos,box);
    }*/
    for (auto& object : objects2) {
        if (object->moved) {
            continue;
        }
        int new_pos = object->getPosition();
        while (true) {
            auto it = rflow.find(new_pos);
            if (it != rflow.end()) {
                new_pos += it->second;
            } else {
                break;
            }
        }
        move(object, new_pos, box);
    }
    // Invalidate cluster data if necessary
    box.clusters_valid = false;
    for(auto& object: objects1){object->moved=false;}
}


void Move::revert(BOX& box){
    for(int i = 0; i<objects1.size();i++){
        //swap_sites(sites2[i],sites1[i],objects1[i],objects2[i],box);
        move(objects1[i],sites1[i],box);
        move(objects2[i],sites2[i],box);
    }
    // Invalidate cluster data if necessary
    box.clusters_valid = false;
}


void Move::swap_sites(int site1, int site2,std::shared_ptr<Object> obj1, std::shared_ptr<Object> obj2, BOX& box) {
    
    // Swap objects in the lattice
    box.set_lattice(site1, obj2);
    box.set_lattice(site2, obj1);
    // Update positions in the objects
    obj1->setPosition(site2); // prev site, next site
    obj2->setPosition(site1);
}
void Move::move(std::shared_ptr<Object> object, int site, BOX& box){
    box.set_lattice(site,object);
    object->setPosition(site);
}