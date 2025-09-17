#ifndef PUZZLE_HPP
#define PUZZLE_HPP
#include <string>
#include <iostream>
#include <vector>

#define BFS "-bfs"
#define IDFS "-idfs"
#define ASTAR "-astar"
#define IDASTAR "idastar"
#define GBFS "-gbfs"
#define EIGHT 8
#define FIFTEEN 15

struct Puzzle {
    Puzzle(std::vector<std::vector<int>> states) : states(states) {
        type = (states[0].size() == 9) ? EIGHT : FIFTEEN;
    }
    int type;
    std::vector<std::vector<int>> states;
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