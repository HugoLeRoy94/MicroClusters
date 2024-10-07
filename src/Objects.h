// Objects.h
#ifndef OBJECTS_H
#define OBJECTS_H

#include <tuple>
#include <vector>
#include <algorithm>  // for std::shuffle
#include <random>

// Forward declaration of BOX class
class BOX;

class Object {
    friend class BOX;  // Declare BOX as a friend class
protected:
    std::tuple<int, int, int> position;

    // Protected constructor and destructor
    Object(const std::tuple<int, int, int>& position_);
    virtual ~Object();

    // Disallow copying and assignment
    Object(const Object&) = delete;
    Object& operator=(const Object&) = delete;


public:
    const std::tuple<int, int, int>& getPosition() const;
    virtual std::vector<std::tuple<int, int, int>> get_positions() const;
    virtual void setPosition(const std::tuple<int, int, int>& value);

    virtual std::tuple<int, int, int> get_site_to_exchange(const BOX& box) const;

    virtual bool isempty() const;
    virtual int Index() const;

    virtual float compute_local_energy(const BOX& box) const = 0;
};

class Empty : public Object {
    friend class BOX;  // Declare BOX as a friend class
protected:
    Empty(const std::tuple<int, int, int>& position_);
    virtual ~Empty();

public:
    virtual bool isempty() const override;
    virtual int Index() const override;
    virtual float compute_local_energy(const BOX& box) const override;
};

class RNA : public Object {
    friend class BOX;  // Declare BOX as a friend class
    friend class Move;
protected:
    RNA(std::vector<std::tuple<int,int,int>>& monomers_);
    virtual ~RNA();
    bool isconnected(int idx, const BOX& box) const;
    bool would_be_connected_after_move(int idx, const std::tuple<int, int, int>& new_pos) const;
    virtual std::tuple<int, int, int> get_site_to_exchange(const BOX& box) const override;
    std::vector<std::tuple<int,int,int>> monomers;
    int get_monomer_index(const std::tuple<int, int, int>& position) const;

public:
    virtual int Index() const override;
    virtual float compute_local_energy(const BOX& box) const override;
    virtual void setPosition(const std::tuple<int, int, int>& value) override;
    virtual std::vector<std::tuple<int, int, int>> RNA::get_positions() const override;

};

class DHH1 : public Object {
    friend class BOX;  // Declare BOX as a friend class
protected:
    DHH1(const std::tuple<int, int, int>& position_);
    virtual ~DHH1();

public:
    virtual int Index() const override;
    virtual std::tuple<int, int, int> get_site_to_exchange(const BOX& box) const override;
    virtual float compute_local_energy(const BOX& box) const override;
};

double distance(std::tuple<int,int,int> site1, std::tuple<int,int,int> site2);
int chebyshev_distance(const std::tuple<int, int, int>& pos1, const std::tuple<int, int, int>& pos2);

#endif // OBJECTS_H
