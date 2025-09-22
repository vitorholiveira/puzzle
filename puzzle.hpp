#ifndef PUZZLE_HPP
#define PUZZLE_HPP

#include <string>
#include <iostream>
#include <vector>
#include <cstdint>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <functional>
#include <chrono>
#include <stack>

#define BFS "-bfs"
#define IDFS "-idfs"
#define ASTAR "-astar"
#define IDASTAR "-idastar"
#define GBFS "-gbfs"

constexpr int TILE_BITS = 4;

constexpr int MOVE_LEFT  = -1;
constexpr int MOVE_RIGHT = 1;
constexpr int MOVE_UP_8  = 3;
constexpr int MOVE_DOWN_8 = -3;
constexpr int MOVE_UP_15  = 4;
constexpr int MOVE_DOWN_15 = -4;

struct SearchStatistics {
    u_int32_t expanded;
    u_int32_t solution_depth;
    double elapsed_seconds;
    double avg_heuristic;
    u_int32_t start_heuristic;
    SearchStatistics(u_int32_t expanded = 0, u_int32_t solution_depth = 0,  double elapsed_seconds = 0.0,
                     double avg_heuristic = 0.0, u_int32_t start_heuristic = 0) :
                     expanded(expanded), solution_depth(solution_depth), elapsed_seconds(elapsed_seconds),
                     avg_heuristic(avg_heuristic), start_heuristic(start_heuristic) { print(); }

    void print() const {
        std::cout << expanded << ","
                  << solution_depth << ","
                  << std::to_string(elapsed_seconds) << ","
                  << std::to_string(avg_heuristic) << ","
                  << start_heuristic << std::endl;
    }
};

class Puzzle {
public:
    Puzzle(const std::vector<std::vector<int>>& initial_states)
        : states(initial_states) {
        int num_tiles = initial_states[0].size();
        bool is_8_puzzle = (num_tiles == 9);
        max_pos = is_8_puzzle ? 8 : 15;
        grid_size = is_8_puzzle ? 3 : 4;
        precompute_manhattan_table();
    }

    void solve(const std::string& algorithm);

    std::vector<std::vector<int>> manhattan_table;
    int max_pos;
    int grid_size;
    u_int64_t goal;

    std::vector<std::vector<int>> states;

    bool solve_bfs(const u_int64_t& start);
    bool solve_idfs(const u_int64_t& start);
    bool solve_astar(const u_int64_t& start);
    bool solve_idastar(const u_int64_t& start);
    bool solve_gbfs(const u_int64_t& start);
    int recursive_dls(const u_int64_t& current_state, 
                      const u_int64_t& parent_state,
                      u_int32_t limit, 
                      u_int32_t current_depth,
                      u_int32_t& n_expanded,
                      u_int32_t& solution_depth);
    void precompute_manhattan_table();
    u_int16_t manhattan_distance(const u_int64_t& state);

    u_int64_t create_goal_state() const;
    u_int64_t vector_to_state(const std::vector<int>& grid_vec) const;
    void print_state(const u_int64_t& state) const;

    std::vector<u_int64_t> expand(const u_int64_t& state) const;


};

#endif // PUZZLE_HPP
