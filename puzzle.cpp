#include "puzzle.hpp"
#include "bucket.hpp"

// OK
void Puzzle::solve(const std::string& algorithm) {
    goal = create_goal_state();
    for(const auto& state : states) {
        if(static_cast<int>(state.size()) != max_pos + 1) continue;
        u_int64_t start = vector_to_state(state);
        if(algorithm == BFS) {
            solve_bfs(start);
        } else if (algorithm == IDFS) {
            solve_idfs(start);
        } else if (algorithm == ASTAR) {
            solve_astar(start);
        } else if (algorithm == IDASTAR) {
            //solve_idastar(start);
        } else if (algorithm == GBFS) {
            solve_gbfs(start);
        }
    }
}

/*
    UTILS
*/

u_int64_t Puzzle::vector_to_state(const std::vector<int>& v) const {
    u_int64_t state = 0;
    for (size_t pos = 0; pos < v.size(); ++pos) {
        u_int64_t tile = static_cast<u_int64_t>(v[pos]) & 0xF; // 4 bits
        state |= (tile << (pos * 4));
    }
    return state;
}

// Goal state: tiles [0,1,2,3,4,5,6,7,8]
u_int64_t Puzzle::create_goal_state() const {
    u_int64_t goal = 0;
    for (u_int8_t pos = 0; pos < max_pos + 1; ++pos) {
        goal |= (static_cast<u_int64_t>(pos) << (pos * 4));
    }
    return goal;
}

void Puzzle::print_state(const u_int64_t& state) const {
    u_int64_t tile;
    for (u_int8_t pos = 0; pos < max_pos + 1; ++pos) {
        if((pos % grid_size) == 0)
            std::cout << std::endl;
        tile = (state >> (pos * 4)) & 0xF;
        std::cout << tile << " ";
    }
    std::cout << std::endl << "--------" << std::endl;
}

std::vector<u_int64_t> Puzzle::expand(const u_int64_t& state) const {
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

u_int16_t Puzzle::manhattan_distance(const u_int64_t& state) {
    int h = 0;
    for (int pos = 0; pos < grid_size * grid_size; pos++) {
        int tile = (state >> (4 * pos)) & 0xF;
        h += manhattan_table[tile][pos];
    }
    return h;
}

void Puzzle::precompute_manhattan_table() {
    manhattan_table.resize(grid_size * grid_size);
    for (int tile = 0; tile < grid_size * grid_size; tile++) {
        manhattan_table[tile].resize(grid_size * grid_size);
        for (int pos = 0; pos < grid_size * grid_size; pos++) {
            if (tile == 0) {
                manhattan_table[tile][pos] = 0; // empty tile
            } else {
                int target_pos = tile; // tile i should be at position i
                int curr_row = pos / grid_size, curr_col = pos % grid_size;
                int target_row = target_pos / grid_size, target_col = target_pos % grid_size;
                manhattan_table[tile][pos] = abs(curr_row - target_row) + abs(curr_col - target_col);
            }
        }
    }
}
/*
    PUZZLE SOLVERS
*/
bool Puzzle::solve_bfs(const u_int64_t& start) {
    // Stores the g value of the already visited nodes
    std::unordered_map<u_int64_t, u_int32_t> visited;
    std::queue<u_int64_t> frontier;
    u_int32_t n_expanded = 0;

    frontier.push(start);
    visited[start] = 0;

    auto start_time = std::chrono::high_resolution_clock::now();

    if (start == goal) {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        double seconds = duration.count() / 1'000'000.0;
        int h_start = manhattan_distance(start);

        std::cout << n_expanded << "," << visited[start] << "," << seconds << "," << 0 << "," << h_start << std::endl;
        return true;
    }

    while (!frontier.empty()) {
        u_int64_t current = frontier.front();
        frontier.pop();

        n_expanded++;

        // Expand neighbors
        for (const auto& child : expand(current)) {
            // If we already saw this child, skip
            if (visited.find(child) != visited.end())
                continue;

            visited[child] = visited[current] + 1;

            if (child == goal) {
                auto end_time = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
                double seconds = duration.count() / 1'000'000.0;
                int h_start = manhattan_distance(start);

                std::cout << n_expanded << "," << visited[child] << "," << seconds << "," << 0 << "," << h_start << std::endl;
                return true;
            }

            frontier.push(child);
        }
    }

    return false;
}

bool Puzzle::solve_idfs(const u_int64_t& start) {
    u_int32_t n_expanded = 0;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Verifica se já começamos no goal
    if (start == goal) {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        double seconds = duration.count() / 1'000'000.0;
        
        SearchStatistics(n_expanded, 0, seconds, 0, manhattan_distance(start));
        return true;
    }

    for (u_int32_t depth_limit = 1; depth_limit <= 50; depth_limit++) {
        u_int32_t solution_depth = 0;
        
        int result = recursive_dls(start, 0, depth_limit, 0, n_expanded, solution_depth);
        
        if (result == 1) {
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
            double seconds = duration.count() / 1'000'000.0;
            
            SearchStatistics(n_expanded, solution_depth, seconds, 0, manhattan_distance(start));
            return true;
        }
        
        if (result == 0) {
            break;
        }
    }
    
    return false;
}

int Puzzle::recursive_dls(
    const u_int64_t& current_state, 
    const u_int64_t& parent_state,
    u_int32_t limit, 
    u_int32_t current_depth,
    u_int32_t& n_expanded,
    u_int32_t& solution_depth) {
    
    // Goal test
    if (current_state == goal) {
        solution_depth = current_depth;
        return 1;
    }
    
    if (limit == 0) {
        return -1;
    }
    
    n_expanded++;
    
    bool cutoff_occurred = false;
    
    for (u_int64_t child : expand(current_state)) {
        if (child == parent_state) {
            continue;
        }
        
        int result = recursive_dls(child, current_state, limit - 1, current_depth + 1, n_expanded, solution_depth);
        
        if (result == -1) {
            cutoff_occurred = true;
        } else if (result == 1) {
            return 1;
        }
    }
    
    if (cutoff_occurred) {
        return -1;
    } else {
        return 0;
    }
}


bool Puzzle::solve_astar(const u_int64_t& start) {
    if (start == goal) {
        SearchStatistics stats(0, 0, 0.0, 0.0, 0);
        return true; // Already solved
    }
    
    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();
    
    BucketQueue<AstarNode> open_set;
    std::unordered_set<u_int64_t> closed_set;
    std::unordered_map<u_int64_t, AstarNode> all_nodes;
    
    int start_h = manhattan_distance(start);
    AstarNode start_node(start, 0, start_h);
    
    open_set.push(start_node);
    all_nodes[start] = start_node;
    
    int nodes_expanded = 0;
    int solution_depth = 0;
    double total_heuristic = 0.0;
    int heuristic_count = 0;
    
    while (!open_set.empty()) {
        AstarNode current = open_set.pop();
        
        if (closed_set.count(current.state)) {
            continue; // Already processed
        }
        
        closed_set.insert(current.state);
        nodes_expanded++;
        
        // Track heuristic values for average calculation
        total_heuristic += current.h;
        heuristic_count++;
        
        if (current.state == goal) {
            // Calculate elapsed time
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
            double elapsed_seconds = duration.count() / 1000000.0;
            
            // Get solution depth from the goal node
            solution_depth = current.g;
            
            // Calculate average heuristic
            double avg_heuristic = (heuristic_count > 0) ? total_heuristic / heuristic_count : 0.0;
            
            // Create and store statistics
            SearchStatistics stats(nodes_expanded, solution_depth, elapsed_seconds, avg_heuristic, start_h);
            
            return true;
        }
        
        for (const auto& child : expand(current.state)) {
            if (closed_set.count(child)) {
                continue;
            }
            
            int new_g = current.g + 1;
            int new_h = manhattan_distance(child);
            
            auto it = all_nodes.find(child);
            if (it == all_nodes.end() || new_g < it->second.g) {
                AstarNode successor_node(child, new_g, new_h, current.state);
                open_set.push(successor_node);
                all_nodes[child] = successor_node;
            }
        }
    }
    
    return false; // No solution found
}


// bool Puzzle::solve_idastar(const u_int64_t& start) {
//     struct Frame {
//         u_int64_t state;
//         int g;
//         int h;
//         size_t child_index;
//         std::vector<u_int64_t> children;
//     };

//     auto start_time = std::chrono::high_resolution_clock::now();

//     int bound = manhattan_distance(start);

//     while (true) {
//         u_int32_t n_expanded = 0;
//         double h_sum = 0.0;
//         int min_cutoff = std::numeric_limits<int>::max();

//         std::stack<Frame> stack;
//         stack.push({start, 0, manhattan_distance(start), 0, {}});

//         while (!stack.empty()) {
//             Frame &node = stack.top();
//             int f = node.g + node.h;

//             // Cutoff check
//             if (f > bound) {
//                 if (f < min_cutoff) min_cutoff = f;
//                 stack.pop();
//                 continue;
//             }

//             // Expansion counting (only first time we see the node)
//             if (node.child_index == 0) {
//                 n_expanded++;
//                 h_sum += node.h;
//             }

//             // Goal check
//             if (node.state == goal) {
//                 auto end_time = std::chrono::high_resolution_clock::now();
//                 auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
//                 double seconds = duration.count() / 1'000'000.0;
//                 double avg_h = (n_expanded > 0) ? (h_sum / n_expanded) : 0.0;

//                 std::cout << n_expanded << "," << node.g << "," << seconds << "," << avg_h << "," << manhattan_distance(start) << std::endl;
//                 return true;
//             }

//             // Expand children lazily
//             if (node.children.empty()) {
//                 node.children = expand(node.state);
//             }

//             if (node.child_index < node.children.size()) {
//                 u_int64_t child = node.children[node.child_index++];
//                 int h_child = manhattan_distance(child);
//                 stack.push({child, node.g + 1, h_child, 0, {}});
//             } else {
//                 stack.pop();
//             }
//         }

//         if (min_cutoff == std::numeric_limits<int>::max()) {
//             return false;
//         }

//         bound = min_cutoff; // increase cutoff
//     }
// }

bool Puzzle::solve_gbfs(const u_int64_t& start) {
    using Node = std::tuple<u_int64_t, u_int32_t, u_int64_t, u_int64_t>; // state, g, insertion order, parent
    std::unordered_set<u_int64_t> visited;
    u_int32_t n_expanded = 0;
    u_int64_t insertion_order = 0; // controla o desempate LIFO

    // Comparador para priority_queue
    auto cmp = [&](const Node& a,
        const Node& b) {
        u_int32_t h_a = static_cast<u_int32_t>(manhattan_distance(std::get<0>(a)));
        u_int32_t h_b = static_cast<u_int32_t>(manhattan_distance(std::get<0>(b)));
        
        if (h_a != h_b) return h_a > h_b;
        if (std::get<1>(a) != std::get<1>(b)) return std::get<1>(a) < std::get<1>(b);
        return std::get<2>(a) < std::get<2>(b);
    };

    std::priority_queue<Node, std::vector<Node>, decltype(cmp)> frontier(cmp);

    auto start_time = std::chrono::high_resolution_clock::now();

    // nó inicial
    int heuristic_calls = 1;
    int heuristic_sum = manhattan_distance(start);
    frontier.push({start, 0, insertion_order++, 0});

    while (frontier.size() > 0) {
        auto [current, g, order, parent] = frontier.top();
        frontier.pop();

        if (visited.find(current) != visited.end()) continue;
        visited.insert(current);

        if(current == goal) {
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
            double seconds = duration.count() / 1'000'000.0;
            u_int64_t h_start = static_cast<u_int64_t>(manhattan_distance(start));
            auto avg_h = float(heuristic_sum) / heuristic_calls;
            SearchStatistics(n_expanded, g, seconds, avg_h, h_start);
            return true;
        }

        n_expanded++;
        for (u_int64_t child : expand(current)) {
            if(child == parent) continue;
            frontier.push({child, g + 1, insertion_order++, current});
            heuristic_calls++;
            heuristic_sum += manhattan_distance(child);
        }
    }
    std::cout << "No solution found." << std::endl;
    return false;
}


