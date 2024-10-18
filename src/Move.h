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
    Move(const std::vector<std::shared_ptr<Object>>& objs1,
        const std::vector<std::shared_ptr<Object>>& objs2,
         const std::vector<int>& old_pos,
         const std::vector<int>& new_pos);

    bool validate(const BOX& box) const;
    void apply(BOX& box);
    void revert(BOX& box);

    // Additional methods as needed

private:
    void swap_sites(int site1, int site2,std::shared_ptr<Object> obj1, std::shared_ptr<Object> obj2,BOX& box);

    std::vector<std::shared_ptr<Object>> objects1;
    std::vector<std::shared_ptr<Object>> objects2;
    std::vector<int> sites1;
    std::vector<int> sites2;
    MoveType move_type;
    // Store any additional state needed to revert the move
};


#endif // MOVE_H