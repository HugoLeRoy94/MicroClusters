// Objects.h
#ifndef OBJECTS_H
#define OBJECTS_H

#include <tuple>
#include <vector>
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
    void setPosition(const std::tuple<int, int, int>& value);

    virtual std::tuple<int, int, int> get_site_to_exchange(const BOX& box) const;

    virtual bool isempty() const;
    virtual int Index() const;
};

class Empty : public Object {
    friend class BOX;  // Declare BOX as a friend class
protected:
    Empty(const std::tuple<int, int, int>& position_);
    virtual ~Empty();

public:
    virtual bool isempty() const override;
    virtual int Index() const override;
};

class RNA : public Object {
    friend class BOX;  // Declare BOX as a friend class
protected:
    RNA(int length, const std::vector<int>& x, const std::vector<int>& y, const std::vector<int>& z);
    virtual ~RNA();

public:
    virtual int Index() const override;
};

class DHH1 : public Object {
    friend class BOX;  // Declare BOX as a friend class
protected:
    DHH1(const std::tuple<int, int, int>& position_);
    virtual ~DHH1();

public:
    virtual int Index() const override;
    virtual std::tuple<int, int, int> get_site_to_exchange(const BOX& box) const override;
};

#endif // OBJECTS_H
