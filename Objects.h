// Objects.h
#ifndef OBJECTS_H
#define OBJECTS_H

#include <tuple>
#include <vector>
#include <random>

// Forward declaration of BOX class
class BOX;

class Object {
protected:
    std::tuple<int, int, int> position;

public:
    Object(const std::tuple<int, int, int>& position_) : position(position_) {}
    virtual ~Object() {}

    const std::tuple<int, int, int>& getPosition() const { return position; }
    void setPosition(const std::tuple<int, int, int>& value) { position = value; }

    virtual bool isempty() const { return false; }
    virtual int Index() const { return 0; }
};

class Empty : public Object {
public:
    Empty(const std::tuple<int, int, int>& position_) : Object(position_) {}
    virtual bool isempty() const override { return true; }
    virtual int Index() const override { return 0; }
};

class RNA : public Object {
public:
    RNA(int length, const std::vector<int>& x, const std::vector<int>& y, const std::vector<int>& z)
        : Object(std::make_tuple(0, 0, 0)) {
        // Implementation of RNA object
    }

    virtual int Index() const override { return 2; }
};

class DHH1 : public Object {
public:
    DHH1(const std::tuple<int, int, int>& position_) : Object(position_) {}

    virtual int Index() const override { return 1; }

    std::tuple<int, int, int> get_site_to_exchange(const BOX& box) const;
};

#endif // OBJECTS_H
