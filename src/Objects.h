// Objects.h
#ifndef OBJECTS_H
#define OBJECTS_H

#include <tuple>
#include <vector>
#include <memory>
#include <algorithm>  // for std::shuffle
#include <random>



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

    virtual int get_site_to_exchange(const BOX& box) const;

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
    //virtual int get_site_to_exchange(const BOX& box) const override;
    std::vector<int> monomers;
    int get_monomer_index(int position) const;

public:
    virtual ~RNA();
    virtual int Index() const override;
    virtual float compute_local_energy(const BOX& box) const override;
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
    virtual int get_site_to_exchange(const BOX& box) const override;
    virtual float compute_local_energy(const BOX& box) const override;
};

int chebyshev_distance(int pos1, int pos2, int L);

#endif // OBJECTS_H
