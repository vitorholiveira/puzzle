#include "puzzle.hpp"

// OK
void Puzzle::solve(const std::string& algorithm) {
    goal = create_goal_state();
    for(const auto& state : states) {
        if(static_cast<int>(state.size()) != max_pos + 1) continue;
        u_int64_t start = vector_to_state(state);
        if(algorithm == BFS) {
            solve_bfs(start);
        } else if (algorithm == IDFS) {
            //solve_idfs(start);
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
    }
    return total;
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

// bool Puzzle::solve_idfs(const u_int64_t& start) {
//     int n_expanded = 0;
//     auto start_time = std::chrono::high_resolution_clock::now();

//     for (int depth_limit = 0; depth_limit < 99999; depth_limit++) {
//         std::stack<std::tuple<u_int64_t, u_int16_t, size_t>> frontier;
//         std::unordered_map<u_int64_t, u_int64_t> parent_map;
//         std::unordered_set<u_int64_t> visited;
        
//         frontier.push({start, 0, 0});
//         bool solution_found = false;

//         while (!frontier.empty()) {
//             auto [current, depth, child_index] = frontier.top();
            
//             if (child_index == 0) {
//                 if (current == goal) {
//                     solution_found = true;
//                     break;
//                 }
                
//                 if (visited.count(current) || depth >= depth_limit) {
//                     visited.erase(current);
//                     frontier.pop();
//                     continue;
//                 }
                
//                 visited.insert(current);
//                 n_expanded++;
//             }
            
//             std::vector<u_int64_t> children = expand(current);
            
//             if (child_index < children.size()) {
//                 u_int64_t child = children[child_index];
//                 std::get<2>(frontier.top()) = child_index + 1;
                
//                 if (visited.find(child) == visited.end()) {
//                     parent_map[child] = current;
//                     frontier.push({child, static_cast<u_int16_t>(depth + 1), 0});
//                 }
//             } else {
//                 visited.erase(current);
//                 frontier.pop();
//             }
//         }

//         if (solution_found) {
//             // Reconstruir caminho e imprimir resultado
//             auto end_time = std::chrono::high_resolution_clock::now();
//             auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
//             double seconds = duration.count() / 1'000'000.0;

//             std::vector<u_int64_t> path;
//             u_int64_t path_node = goal;
//             while (parent_map.count(path_node) && path_node != start) {
//                 path.push_back(path_node);
//                 path_node = parent_map[path_node];
//             }
//             path.push_back(start);
//             std::reverse(path.begin(), path.end());

//             std::cout << n_expanded << "," << path.size() - 1 << "," << seconds << "," << 0 << "," << static_cast<int>(manhattan_distance(start)) << std::endl;

//             return true;
//         }
//     }
//     return false;
// }

bool Puzzle::solve_astar(const u_int64_t& start) {
    using Node = std::tuple<u_int64_t, u_int32_t, u_int64_t>; // state, g, insertion order
    std::unordered_map<u_int64_t, u_int64_t> parent;          // for path reconstruction
    std::unordered_map<u_int64_t, u_int32_t> g_score;         // best g found so far (discovered nodes)
    std::unordered_set<u_int64_t> closed;                     // expanded states
    u_int32_t n_expanded = 0;
    u_int64_t insertion_order = 0;

    // Priority by f = g + h, then lower h, then LIFO
    auto cmp = [&](const Node& a, const Node& b) {
        auto f = [&](const Node& n) {
            return std::get<1>(n) + (u_int32_t)manhattan_distance(std::get<0>(n));
        };
        u_int32_t f_a = f(a), f_b = f(b);
        if (f_a != f_b) return f_a > f_b;                     // lower f first
        u_int32_t h_a = (u_int32_t)manhattan_distance(std::get<0>(a));
        u_int32_t h_b = (u_int32_t)manhattan_distance(std::get<0>(b));
        if (h_a != h_b) return h_a > h_b;                     // then lower h
        return std::get<2>(a) < std::get<2>(b);               // then newer first (LIFO)
    };

    std::priority_queue<Node, std::vector<Node>, decltype(cmp)> frontier(cmp);

    auto start_time = std::chrono::high_resolution_clock::now();

    // initialize
    g_score[start] = 0;
    parent[start]  = start;
    frontier.push({start, 0, insertion_order++});

    while (!frontier.empty()) {
        auto [current, g, order] = frontier.top();
        frontier.pop();

        if (closed.count(current)) continue;  // already expanded with best cost
        closed.insert(current);

        // Goal test when popping ensures optimality
        if (current == goal) {
            auto end_time = std::chrono::high_resolution_clock::now();
            double seconds = std::chrono::duration_cast<std::chrono::microseconds>(
                                  end_time - start_time).count() / 1'000'000.0;

            // average heuristic over all discovered nodes (keys in g_score)
            double avg_h = 0.0;
            if (!g_score.empty()) {
                unsigned long long heuristic_sum = 0ULL;
                for (const auto &p : g_score) {
                    heuristic_sum += static_cast<u_int64_t>(manhattan_distance(p.first));
                }
                avg_h = double(heuristic_sum) / double(g_score.size());
            }

            u_int32_t h_start = static_cast<u_int32_t>(manhattan_distance(start));
            std::cout << n_expanded << "," << g << "," << seconds << "," << avg_h << "," << h_start << std::endl;
            return true;
        }

        n_expanded++;

        for (const auto& child : expand(current)) {
            u_int32_t g_new = g + 1;
            if (!g_score.count(child) || g_new < g_score[child]) {
                g_score[child] = g_new;
                parent[child]  = current;
                frontier.push({child, g_new, insertion_order++});
            }
        }
    }

    return false; // no solution
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
    using Node = std::tuple<u_int64_t, u_int32_t, u_int64_t>; // state, g, insertion order
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

    Node final_state;

    auto start_time = std::chrono::high_resolution_clock::now();

    // nó inicial
    frontier.push({start, 0, insertion_order++});
    visited.insert(start);

    if (start == goal) {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        double seconds = duration.count() / 1'000'000.0;
        int h_start = static_cast<int>(manhattan_distance(start));

        std::cout << n_expanded << "," << 0 << "," << seconds << "," << h_start << "," << h_start << std::endl;
        return true;
    }
    bool found_solution = false;

    while (!found_solution) {
        auto [current, g, order] = frontier.top();
        frontier.pop();
        n_expanded++;
        for (u_int64_t child : expand(current)) {
            if (visited.find(child) != visited.end()) continue;
            visited.insert(child);
            if (child == goal) {
                found_solution = true;
                final_state = {child, g + 1, insertion_order++};
            }
            frontier.push({child, g + 1, insertion_order++});
        }
    }

    if(found_solution){
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        double seconds = duration.count() / 1'000'000.0;
        u_int64_t h_start = static_cast<u_int64_t>(manhattan_distance(start));
        u_int64_t heuristic_sum = 0;
        for (const auto &state : visited) {
            heuristic_sum += static_cast<u_int64_t>(manhattan_distance(state));
        }

        std::cout << n_expanded << "," << std::get<1>(final_state) << "," << seconds << "," << float(heuristic_sum) / visited.size() << "," << h_start << std::endl;
        return true;
    }

    return false;
}


