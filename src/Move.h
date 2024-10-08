#ifndef MOVE_H
#define MOVE_H

#include <tuple>
#include <memory>
#include "Objects.h"
#include "BOX.h"

enum class MoveType {
    Swap,
    // Other move types can be added here
};

class Move {
public:
    Move(std::shared_ptr<Object> obj1, int pos1,
         std::shared_ptr<Object> obj2, int pos2,
         MoveType type);

    bool validate(const BOX& box) const;
    void apply(BOX& box);
    void revert(BOX& box);

    // Additional methods as needed

private:
    std::shared_ptr<Object> object1;
    std::shared_ptr<Object> object2;
    int site1;
    int site2;
    MoveType move_type;

    // Store any additional state needed to revert the move
};

#endif // MOVE_H
