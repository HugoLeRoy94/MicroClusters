#include <tuple>
#include "Objects.h"
#include "BOX.h"


enum class MoveType {
    Swap,
    // Other move types can be added here
};

class Move {
public:
    Move(Object* obj1, const std::tuple<int, int, int>& pos1,
         Object* obj2, const std::tuple<int, int, int>& pos2,
         MoveType type);

    bool validate(const BOX& box) const;
    void apply(BOX& box);
    void revert(BOX& box);

    // Additional methods as needed

private:
    Object* object1;
    Object* object2;
    std::tuple<int, int, int> site1;
    std::tuple<int, int, int> site2;
    MoveType move_type;

    // Store any additional state needed to revert the move
};
