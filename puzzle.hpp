#ifndef PUZZLE_HPP
#define PUZZLE_HPP

#include <string>
#include <iostream>
#include <vector>
#include <cstdint>
#include <bitset>
#include <unordered_set>

constexpr int TILE_BITS = 4;

constexpr size_t BITS_GRID_8 = 36;   // 9 tiles * 4 bits
constexpr size_t BITS_GRID_15 = 64;  // 16 tiles * 4 bits

constexpr int MOVE_LEFT  = -1;
constexpr int MOVE_RIGHT = 1;
constexpr int MOVE_UP_8  = 3;
constexpr int MOVE_DOWN_8 = -3;
constexpr int MOVE_UP_15  = 4;
constexpr int MOVE_DOWN_15 = -4;

template <size_t BITS_GRID>
class Puzzle {
public:
    explicit Puzzle(const std::vector<std::vector<int>>& initial_states)
        : states(initial_states), 
          GOAL_PUZZLE(create_goal_state()) {
        bool is_8_puzzle = (BITS_GRID == BITS_GRID_8);
        max_pos = is_8_puzzle ? 8 : 15;
        grid_size = is_8_puzzle ? 3 : 4;
    }

    bool solve_bfs();
    bool solve_idfs();
    bool solve_astar();
    bool solve_iastar();
    bool solve_gbfs();

private:
    int max_pos;
    int grid_size;
    std::vector<std::vector<int>> states;
    std::bitset<BITS_GRID> GOAL_PUZZLE;

    std::bitset<BITS_GRID> create_goal_state() const;
    std::bitset<BITS_GRID> vector_to_bitset(const std::vector<int>& grid_vec) const;

    std::vector<std::bitset<BITS_GRID>> expand(const std::bitset<BITS_GRID>& grid) const;
};

#endif // PUZZLE_HPP
