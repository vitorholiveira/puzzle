#ifndef PUZZLE_HPP
#define PUZZLE_HPP
#include <string>
#include <iostream>
#include <vector>
#include <cstdint>
#include <bitset>

#define BFS "-bfs"
#define IDFS "-idfs"
#define ASTAR "-astar"
#define IDASTAR "idastar"
#define GBFS "-gbfs"
#define NUM_BITS 4
#define EIGHT 8
#define FIFTEEN 15

struct Puzzle {
    Puzzle(std::vector<std::vector<std::bitset<NUM_BITS>>> &states) : states(states) {
        type = (states[0].size() == 9) ? EIGHT : FIFTEEN;
    }
    int type;
    std::vector<std::vector<std::bitset<NUM_BITS>>> states;
};

class PuzzleSolver {
public:

    // static bfs_eight

    // static idfs_eight

    // static astar_eight

    // static iastar_eight

    // static gbfs_eight

    // static astar_fifteen

};

#endif // PUZZLE_HPP