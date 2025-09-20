#include "puzzle.hpp"

// OK
template<size_t BITS_GRID>
void Puzzle<BITS_GRID>::solve(const std::string& algorithm) {
    goal = create_goal_state();
    for(const auto& state : states) {
        if(state.size() != max_pos + 1) continue;
        u_int64_t start = vector_to_state(state);
        if(algorithm == BFS) {
            solve_bfs(start);
        } else if (algorithm == IDFS) {
            solve_idfs(start);
        } else if (algorithm == ASTAR) {
            solve_astar(start);
        } else if (algorithm == IDASTAR) {
            solve_iastar(start);
        } else if (algorithm == GBFS) {
            solve_gbfs(start);
        }
    }
}

/*
    UTILS
*/

template<size_t BITS_GRID>
u_int64_t Puzzle<BITS_GRID>::vector_to_state(const std::vector<int>& v) const {
    u_int64_t state = 0;
    for (size_t pos = 0; pos < v.size(); ++pos) {
        u_int64_t tile = static_cast<u_int64_t>(v[pos]) & 0xF; // 4 bits
        state |= (tile << (pos * 4));
    }
    return state;
}

// Goal state: tiles [0,1,2,3,4,5,6,7,8]
template<size_t BITS_GRID>
u_int64_t Puzzle<BITS_GRID>::create_goal_state() const {
    u_int64_t goal = 0;
    for (u_int8_t pos = 0; pos < max_pos + 1; ++pos) {
        goal |= (static_cast<u_int64_t>(pos) << (pos * 4));
    }
    return goal;
}

template<size_t BITS_GRID>
void Puzzle<BITS_GRID>::print_state(const u_int64_t& state) const {
    u_int64_t tile;
    for (u_int8_t pos = 0; pos < max_pos + 1; ++pos) {
        if((pos % grid_size) == 0)
            std::cout << std::endl;
        tile = (state >> (pos * 4)) & 0xF;
        std::cout << tile << " ";
    }
    std::cout << std::endl << "--------" << std::endl;
}

template<size_t BITS_GRID>
std::vector<u_int64_t> Puzzle<BITS_GRID>::expand(const u_int64_t& state) const {
    std::vector<u_int64_t> children;

    int blank_pos = -1;
    for (u_int8_t pos = 0; pos < max_pos + 1; ++pos) {
        if (((state >> (pos * 4)) & 0xF) == 0) {
            blank_pos = pos;
            break;
        }
    }

    if (blank_pos == -1) return children;

    u_int8_t row = blank_pos / grid_size;
    u_int8_t col = blank_pos % grid_size;

    std::vector<int> moves = {-grid_size, -1, 1, grid_size}; // Up, Left, Right, Down
    std::vector<bool> valid_moves = {
        (row > 0),             // Up
        (col > 0),             // Left
        (col < grid_size - 1), // Right
        (row < grid_size - 1)  // Down
    };

    for (int i = 0; i < 4; ++i) {
        if (!valid_moves[i]) continue;
        int new_pos = blank_pos + moves[i];
        
        u_int64_t new_state = state;
        
        // Extract tile to be moved
        u_int64_t tile_value = (state >> (new_pos * 4)) & 0xF;
        
        // Clear both the new tile position and the blank's old position
        u_int64_t mask = ~(
            (static_cast<u_int64_t>(0xF) << (new_pos * 4)) |
            (static_cast<u_int64_t>(0xF) << (blank_pos * 4))
        );
        new_state &= mask;
        
        // Place the tile in the blank's old position
        new_state |= tile_value << (blank_pos * 4);
        
        children.push_back(new_state);
    }
    return children;
}

template<size_t BITS_GRID>
u_int16_t Puzzle<BITS_GRID>::manhattan_distance(const u_int64_t& state) const {
    u_int16_t total = 0;
    int pos;
    for (pos = 0; pos < max_pos + 1; ++pos) {
        int tile = (state >> (pos * TILE_BITS)) & 0xF;  // extract 4 bits
        if (tile == 0) continue; // skip blank

        int cur_row  = pos / grid_size;
        int cur_col  = pos % grid_size;
        int goal_row = tile / grid_size;
        int goal_col = tile % grid_size;


        total += static_cast<u_int16_t>(std::abs(cur_row - goal_row) + std::abs(cur_col - goal_col));
        //std::cout << std::endl << "tile: " << static_cast<int>(tile) << " pos: " << static_cast<int>(pos) << " h = " << std::abs(cur_row - goal_row) + std::abs(cur_col - goal_col) << std::endl;
    }
    //std::cout << "total: " << total << " pos: " << pos << std::endl;
    return total;
}

/*
    PUZZLE SOLVERS
*/
template<size_t BITS_GRID>
bool Puzzle<BITS_GRID>::solve_bfs(const u_int64_t& start) const {
    std::unordered_map<u_int64_t, u_int32_t> expansion_counts;
    std::queue<u_int64_t> frontier;
    u_int32_t n_expanded = 0;

    frontier.push(start);
    expansion_counts[start] = 0;

    auto start_time = std::chrono::high_resolution_clock::now();

    if (start == goal) {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        double seconds = duration.count() / 1'000'000.0;

        std::cout << "Solution found with BFS:" << std::endl;
        std::cout << n_expanded << "," << expansion_counts[start] << "," << seconds << "," << 0 << "," << static_cast<int>(manhattan_distance(start)) << std::endl;
        return true;
    }

    while (!frontier.empty()) {
        u_int64_t current = frontier.front();
        frontier.pop();

        n_expanded++;

        //std::cout << "Exploring node:\n";
        //print_state(current);

        // Expand neighbors
        for (u_int64_t child : expand(current)) {
            //std::cout << "Expanding to:\n";
            //print_state(child);
            // If we already saw this child, skip
            if (expansion_counts.find(child) != expansion_counts.end())
                continue;

            expansion_counts[child] = expansion_counts[current] + 1;

            if (child == goal) {
                auto end_time = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
                double seconds = duration.count() / 1'000'000.0;

                std::cout << "Solution found with BFS:" << std::endl;
                std::cout << n_expanded << "," << expansion_counts[child] << "," << seconds << "," << 0 << "," << static_cast<int>(manhattan_distance(start)) << std::endl;
                return true;
            }

            frontier.push(child);
        }
    }

    std::cout << "Nenhuma solução encontrada.\n";
    return false;
}

template<size_t BITS_GRID>
bool Puzzle<BITS_GRID>::solve_idfs(const u_int64_t& start) const {
    int n_expanded = 0;

    auto start_time = std::chrono::high_resolution_clock::now();

    int depth_limit = 0;
    // The main IDFS loop: iterate through increasing depth limits.
    while (true) {
        depth_limit++;
        // The stack stores pairs of {state, current_depth}
        std::stack<std::pair<u_int64_t, int>> frontier;
        // Map to store parent-child relationships for path reconstruction.
        std::unordered_map<u_int64_t, u_int64_t> parent_map;
        // Map to track the shallowest depth a node was visited at in this DLS run.
        // This prevents cycles and redundant exploration within one iteration.
        std::unordered_map<u_int64_t, int> visited_depth;

        frontier.push({start, 0});
        visited_depth[start] = 0;
        
        bool solution_found_this_iteration = false;

        while (!frontier.empty()) {
            auto [current, depth] = frontier.top();
            frontier.pop();
            n_expanded++;


            // Goal check
            if (current == goal) {
                solution_found_this_iteration = true;
                break; // Exit the while loop
            }

            // If at the depth limit, do not expand this node.
            if (depth >= depth_limit) {
                continue;
            }

            // Expand neighbors (children)
            // This is optional but ensures consistency.
            for (u_int64_t& child : expand(current)) {
                int new_depth = depth + 1;
                // If we haven't seen this child before, or if we found a shorter path to it,
                // add it to the frontier.
                if (visited_depth.find(child) == visited_depth.end() || visited_depth[child] > new_depth) {
                    visited_depth[child] = new_depth;
                    parent_map[child] = current;
                    frontier.push({child, new_depth});
                }
            }
        }
        // Iterative DLS Ends

        if (solution_found_this_iteration) {
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
            double seconds = duration.count() / 1'000'000.0;

            // Reconstruct the path from the parent map
            std::vector<u_int64_t> path;
            u_int64_t path_node = goal;
            // The check 'path_node != start' handles the case where start == goal
            while (parent_map.count(path_node) && path_node != start) {
                path.push_back(path_node);
                path_node = parent_map[path_node];
            }
            path.push_back(start);
            std::reverse(path.begin(), path.end());

            std::cout << "Solution found with IDFS:" << std::endl;
            std::cout << n_expanded << "," << path.size() << "," << seconds << "," << 0 << "," << static_cast<int>(manhattan_distance(start)) << std::endl;

            return true;
        }
    }

    // This part is generally unreachable if a solution exists.
    std::cout << "No solution found.\n";
    return false;
}

template<size_t BITS_GRID>
bool Puzzle<BITS_GRID>::solve_astar(const u_int64_t& start) const{
    // TODO: Setup nodes
    std::cout << start << std::endl;
    // TODO: implement actual ASTAR search using expand()
    return false;
}

template<size_t BITS_GRID>
bool Puzzle<BITS_GRID>::solve_iastar(const u_int64_t& start) const{
    // TODO: Setup nodes
    std::cout << start << std::endl;
    // TODO: implement actual IASTAR search using expand()
    return false;
}

template<size_t BITS_GRID>
bool Puzzle<BITS_GRID>::solve_gbfs(const u_int64_t& start) const{
    // A map to store states we've seen and a pair:
    // - First element is the path length;
    // - Second element is the sum of all heuristic value.
    std::unordered_map<u_int64_t, std::pair<u_int32_t, u_int64_t>> expansion_counts;

    // Use a min-priority queue to store states to explore.
    // The pair stores <heuristic_value, state_representation>
    // std::greater<std::pair<u_int32_t, u_int64_t>> makes it a min-heap based on the first element (the heuristic).
    std::priority_queue<std::pair<u_int16_t, u_int64_t>,
                        std::vector<std::pair<u_int16_t, u_int64_t>>,
                        std::greater<std::pair<u_int16_t, u_int64_t>>> frontier;

    u_int32_t n_expanded = 0;

    u_int16_t start_heuristic = manhattan_distance(start);

    frontier.push({start_heuristic, start});
    expansion_counts[start] = {0, static_cast<u_int64_t>(start_heuristic)};

    auto start_time = std::chrono::high_resolution_clock::now();

    while (!frontier.empty()) {
        // Get the state with the lowest heuristic value (the "best" one).
        // The 'top()' method returns a reference to the top element.
        u_int64_t current = frontier.top().second;
        u_int16_t current_heuristic = frontier.top().first;

        // Then, remove it from the frontier.
        frontier.pop();

        std::cout << "==> Exploring node:\n";
        print_state(current);
        std::cout << "heuristic value: " << static_cast<int>(current_heuristic) << std::endl << std::endl;

        // If we've already expanded this state, skip it to avoid redundant work.
        if (expansion_counts.find(current) != expansion_counts.end() && expansion_counts[current].first < n_expanded) {
            continue;
        }
        
        n_expanded++;

        // Expand neighbors
        for (u_int64_t child : expand(current)) {
            // Check if we've seen this child before.
            if (expansion_counts.find(child) != expansion_counts.end()) continue; 

            std::cout << "Expanding to:\n";
            print_state(child);
            u_int16_t h = manhattan_distance(child);
            std::cout << "heuristic value: " << static_cast<int>(h) << std::endl << std::endl;

            // If not, calculate its heuristic and add to the frontier.
            u_int16_t child_heuristic = manhattan_distance(child);
            frontier.push({child_heuristic, child});
            expansion_counts[child] = {expansion_counts[current].first + 1, expansion_counts[current].second + static_cast<u_int64_t>(child_heuristic)};

            if (child == goal) {
                auto end_time = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
                double seconds = duration.count() / 1'000'000.0;

                std::cout << "Solution found with GBFS:" << std::endl;
                std::cout << n_expanded << "," << expansion_counts[current].first << "," << seconds << "," <<float(expansion_counts[current].second) / expansion_counts[current].first << "," << start_heuristic << std::endl;
                return true;
            }
        }
    }

    std::cout << "Nenhuma solução encontrada.\n";
    return false;
}

// --- Explicit instantiations (required if you use .cpp for templates) ---
template class Puzzle<BITS_GRID_8>;
template class Puzzle<BITS_GRID_15>;