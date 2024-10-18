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

class Object :public std::enable_shared_from_this<Object>{
    friend class BOX;  // Declare BOX as a friend class
    // Protected constructor and destructor    
protected:
    int position;    

    // Disallow copying and assignment
    Object(const Object&) = delete;
    Object& operator=(const Object&) = delete;    

public:
Object(int position_);
    virtual ~Object();
    int getPosition() const;
    virtual std::vector<int> get_positions() const;
    virtual void setPosition(int old_pos, int new_pos );

    virtual void get_swap_site_candidate(std::vector<std::shared_ptr<Object>>& object1,
                                        std::vector<std::shared_ptr<Object>>& object2,
                                        std::vector<int>& sites1,
                                        std::vector<int>& sites2,
                                        const BOX& box, std::mt19937& rng) =0;

    virtual bool isempty() const;
    virtual int Index() const;

    virtual float compute_local_energy(const BOX& box) const = 0;

    virtual bool would_be_connected_after_move(int idx, int new_pos, int L) const = 0;
    virtual int get_monomer_index(int site) const = 0;
    
};

class Empty : public Object{
    friend class BOX;  // Declare BOX as a friend class
protected:
public:
    Empty(int position_);
    virtual ~Empty();
    virtual bool isempty() const override;
    virtual int Index() const override;
    virtual float compute_local_energy(const BOX& box) const override;
    virtual bool would_be_connected_after_move(int idx, int new_pos, int L) const override;
    virtual void get_swap_site_candidate(std::vector<std::shared_ptr<Object>>& object1,
                                        std::vector<std::shared_ptr<Object>>& object2,
                                        std::vector<int>& sites1,
                                        std::vector<int>& sites2,
                                         const BOX& box, std::mt19937& rng) override;
    virtual int get_monomer_index(int site) const;
};

class RNA : public Object{
    friend class BOX;  // Declare BOX as a friend class
    friend class Move;
protected:    
    bool isconnected(int idx, const BOX& box) const;
    virtual void get_swap_site_candidate(std::vector<std::shared_ptr<Object>>& object1,
                                        std::vector<std::shared_ptr<Object>>& object2,
                                        std::vector<int>& sites1,
                                        std::vector<int>& sites2,
                                         const BOX& box, std::mt19937& rng) override;
    std::vector<int> monomers;
    
public:
    RNA(std::vector<int>& monomers_);
    virtual ~RNA();
    virtual int Index() const override;
    virtual float compute_local_energy(const BOX& box) const override;
    virtual void setPosition(int old_pos, int new_pos) override;
    virtual std::vector<int> get_positions() const override;
    int get_monomer_index(int position) const;
    bool would_be_connected_after_move(int idx, int new_pos, int L) const;
};

class DHH1 : public Object {
    friend class BOX;  // Declare BOX as a friend class
protected:    
public:
    DHH1(int position_);
    virtual ~DHH1();
    virtual int Index() const override;
    virtual void get_swap_site_candidate(std::vector<std::shared_ptr<Object>>& object1,
                                        std::vector<std::shared_ptr<Object>>& object2,
                                        std::vector<int>& sites1,
                                        std::vector<int>& sites2,
                                        const BOX& box, std::mt19937& rng) override;
    virtual float compute_local_energy(const BOX& box) const override;
    virtual bool would_be_connected_after_move(int idx, int new_pos, int L) const override;
    virtual int get_monomer_index(int site) const;
    };

#endif // OBJECTS_H
