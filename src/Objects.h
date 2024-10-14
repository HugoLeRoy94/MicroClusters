// Objects.h
#ifndef OBJECTS_H
#define OBJECTS_H

#include <tuple>
#include <vector>
#include <memory>
#include <algorithm>  // for std::shuffle
#include <random>
#include "Utilities.h"



// Forward declaration of BOX class
class BOX;

class Object {
    friend class BOX;  // Declare BOX as a friend class
protected:
    int position;

    // Protected constructor and destructor
    Object(int position_);    

    // Disallow copying and assignment
    Object(const Object&) = delete;
    Object& operator=(const Object&) = delete;

public:    
    virtual ~Object();
    int getPosition() const;
    virtual std::vector<int> get_positions() const;
    virtual void setPosition(int value);

    virtual int get_swap_site_candidate(int idx,const BOX& box, std::mt19937& rng) const;

    virtual bool isempty() const;
    virtual int Index() const;

    virtual float compute_local_energy(const BOX& box) const = 0;
};

class Empty : public Object {
    friend class BOX;  // Declare BOX as a friend class
protected:
    Empty(int position_);
    static std::shared_ptr<Empty> make_shared_ptr(int position_);  
public:
    virtual ~Empty();
    virtual bool isempty() const override;
    virtual int Index() const override;
    virtual float compute_local_energy(const BOX& box) const override;
};

class RNA : public Object {
    friend class BOX;  // Declare BOX as a friend class
    friend class Move;
protected:
    RNA(std::vector<int>& monomers_);
    static std::shared_ptr<RNA> make_shared_ptr(std::vector<int>& monomers_);
    bool isconnected(int idx, const BOX& box) const;
    bool would_be_connected_after_move(int idx, int new_pos, int L) const;    
    virtual int get_swap_site_candidate(int idx, const BOX& box, std::mt19937& rng) const override;
    std::vector<int> monomers;
    int get_monomer_index(int position) const;

public:
    virtual ~RNA();
    virtual int Index() const override;
    virtual float compute_local_energy(const BOX& box) const override;
    /*inline float compute_local_energy(const BOX& box) const override{
    float local_energy = 0.0f;
    for (const auto& idx : monomers) {
        auto neighbors = box.get_neighbors(idx);
        for (const auto& nidx : neighbors) {
            auto neighbor_obj = box.get_lattice(nidx);
            local_energy -= box.E[Index()][neighbor_obj->Index()];
        }
    }
    return local_energy;
    }*/
    //virtual void setPosition(int value) override;
    virtual std::vector<int> get_positions() const override;
};

class DHH1 : public Object {
    friend class BOX;  // Declare BOX as a friend class
protected:
    DHH1(int position_);    
    static std::shared_ptr<DHH1> make_shared_ptr(int position_);
public:
    virtual ~DHH1();
    virtual int Index() const override;
    virtual int get_swap_site_candidate(int idx,const BOX& box, std::mt19937& rng) const override;
    virtual float compute_local_energy(const BOX& box) const override;
    /*inline float compute_local_energy(const BOX& box) const override {
        auto neighbors = box.get_neighbors(getPosition());
        float local_energy = 0.0f;
        bool has_neigh = false;
        for (const auto& nidx : neighbors) {
            auto neighbor_obj = box.get_lattice(nidx);
            local_energy -= box.E[Index()][neighbor_obj->Index()];
            if (neighbor_obj->Index() == 1) {
                has_neigh = true;
            }
        }
        if (has_neigh) {
            local_energy -= box.Evalence;
        }
        return local_energy;
    }*/
};

#endif // OBJECTS_H
