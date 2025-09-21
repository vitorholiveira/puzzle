#include "puzzle.hpp"

// OK
template<size_t BITS_GRID>
void Puzzle<BITS_GRID>::solve(const std::string& algorithm) {
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
            solve_idastar(start);
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
        // The stack stores pairs of {state, current_depth}
        std::stack<std::pair<u_int64_t, u_int16_t>> frontier;
        // Map to store parent-child relationships for path reconstruction.
        std::unordered_map<u_int64_t, u_int64_t> parent_map;
        // Map to track the shallowest depth a node was visited at in this DLS run.
        // This prevents cycles and redundant exploration within one iteration.
        std::unordered_map<u_int64_t, u_int16_t> visited_depth;

        frontier.push({start, 0});
        visited_depth[start] = 0;
        
        bool solution_found_this_iteration = false;
        

        while (!frontier.empty()) {
            auto [current, depth] = frontier.top();
            frontier.pop();        

            // std::cout << "===> " << depth << std::endl;
            // std::cout << "Exploring node:\n";
            // print_state(current);


            // Goal check
            if (current == goal) {
                solution_found_this_iteration = true;
                break; // Exit the while loop
            }

            // If at the depth limit, do not expand this node.
            if (depth >= depth_limit) {
                continue;
            }

            n_expanded++;
            
            std::vector<u_int64_t> mock_children = expand(current);
            std::reverse(mock_children.begin(), mock_children.end());
            // Expand neighbors (children)
            // This is optional but ensures consistency.
            for (const u_int64_t& child : mock_children) {
                int new_depth = depth + 1;
                // std::cout << "Exploring to:\n";
                // print_state(child);
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
            std::cout << n_expanded << "," << path.size() - 1 << "," << seconds << "," << 0 << "," << static_cast<int>(manhattan_distance(start)) << std::endl;

            return true;
        }
        depth_limit++;
    }

    // This part is generally unreachable if a solution exists.
    std::cout << "No solution found.\n";
    return false;
}

template<size_t BITS_GRID>
bool Puzzle<BITS_GRID>::solve_astar(const u_int64_t& start) const {
    using State = u_int64_t;
    using Cost = u_int32_t;

    struct Node {
        State state;
        Cost g; // cost so far
        Cost f; // g + h

        bool operator>(const Node& other) const {
            return f > other.f; // min-heap
        }
    };

    std::unordered_map<State, Cost> g_costs;
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> frontier;

    g_costs[start] = 0;
    frontier.push({start, 0, static_cast<Cost>(manhattan_distance(start))});

    u_int32_t n_expanded = 0;
    double h_sum = 0.0;  // track heuristic sum

    auto start_time = std::chrono::high_resolution_clock::now();

    if (start == goal) {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        double seconds = duration.count() / 1'000'000.0;

        double avg_h = static_cast<int>(manhattan_distance(start)); // only one state
        std::cout << "Solution found with A*:" << std::endl;
        std::cout << n_expanded << "," << 0 << "," << seconds << "," << avg_h << "," << static_cast<int>(manhattan_distance(start)) << std::endl;
        return true;
    }

    while (!frontier.empty()) {
        Node current = frontier.top();
        frontier.pop();

        n_expanded++;
        int h_val = static_cast<int>(manhattan_distance(current.state));
        h_sum += h_val;

        if (current.state == goal) {
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
            double seconds = duration.count() / 1'000'000.0;

            double avg_h = (n_expanded > 0) ? (h_sum / n_expanded) : 0.0;

            std::cout << "Solution found with A*:" << std::endl;
            std::cout << n_expanded << "," << current.g << "," << seconds << "," << avg_h << "," << static_cast<int>(manhattan_distance(start)) << std::endl;
            return true;
        }

        for (State child : expand(current.state)) {
            Cost tentative_g = current.g + 1;

            if (g_costs.find(child) == g_costs.end() || tentative_g < g_costs[child]) {
                g_costs[child] = tentative_g;
                Cost h = static_cast<Cost>(manhattan_distance(child));
                frontier.push({child, tentative_g, tentative_g + h});
            }
        }
    }

    std::cout << "Nenhuma solução encontrada.\n";
    return false;
}

template<size_t BITS_GRID>
bool Puzzle<BITS_GRID>::solve_idastar(const u_int64_t& start) const {
    struct Frame {
        u_int64_t state;
        int g;
        int h;
        size_t child_index;
        std::vector<u_int64_t> children;
    };

    auto start_time = std::chrono::high_resolution_clock::now();

    int bound = manhattan_distance(start);

    while (true) {
        u_int32_t n_expanded = 0;
        double h_sum = 0.0;
        int min_cutoff = std::numeric_limits<int>::max();

        std::stack<Frame> stack;
        stack.push({start, 0, manhattan_distance(start), 0, {}});

        while (!stack.empty()) {
            Frame &node = stack.top();
            int f = node.g + node.h;

            // Cutoff check
            if (f > bound) {
                if (f < min_cutoff) min_cutoff = f;
                stack.pop();
                continue;
            }

            // Expansion counting (only first time we see the node)
            if (node.child_index == 0) {
                n_expanded++;
                h_sum += node.h;
            }

            // Goal check
            if (node.state == goal) {
                auto end_time = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
                double seconds = duration.count() / 1'000'000.0;
                double avg_h = (n_expanded > 0) ? (h_sum / n_expanded) : 0.0;

                std::cout << "Solution found with IDA*:" << std::endl;
                std::cout << n_expanded << "," << node.g << "," << seconds << "," << avg_h << "," << manhattan_distance(start) << std::endl;
                return true;
            }

            // Expand children lazily
            if (node.children.empty()) {
                node.children = expand(node.state);
            }

            if (node.child_index < node.children.size()) {
                u_int64_t child = node.children[node.child_index++];
                int h_child = manhattan_distance(child);
                stack.push({child, node.g + 1, h_child, 0, {}});
            } else {
                stack.pop();
            }
        }

        if (min_cutoff == std::numeric_limits<int>::max()) {
            std::cout << "Nenhuma solução encontrada.\n";
            return false;
        }

        bound = min_cutoff; // increase cutoff
    }
}

template<size_t BITS_GRID>
bool Puzzle<BITS_GRID>::solve_gbfs(const u_int64_t& start) const {
    // expansion_counts guarda o g (número de passos até o nó)
    std::unordered_map<u_int64_t, u_int32_t> expansion_counts;
    u_int32_t n_expanded = 0;
    u_int64_t insertion_order = 0; // controla o desempate LIFO

    // Comparador para priority_queue
    auto cmp = [&](const std::tuple<u_int64_t,u_int32_t,u_int64_t>& a,
                   const std::tuple<u_int64_t,u_int32_t,u_int64_t>& b) {
        // a = (estado, g, ordem)
        u_int64_t state_a = std::get<0>(a);
        u_int32_t g_a = std::get<1>(a);
        u_int64_t order_a = std::get<2>(a);

        u_int64_t state_b = std::get<0>(b);
        u_int32_t g_b = std::get<1>(b);
        u_int64_t order_b = std::get<2>(b);

        u_int32_t h_a = manhattan_distance(state_a);
        u_int32_t h_b = manhattan_distance(state_b);

        if (h_a != h_b) return h_a > h_b;   // menor h tem prioridade
        if (g_a != g_b) return g_a < g_b;   // maior g tem prioridade
        return order_a < order_b;           // LIFO: último inserido tem prioridade
    };

    std::priority_queue<
        std::tuple<u_int64_t,u_int32_t,u_int64_t>, // estado, g, ordem
        std::vector<std::tuple<u_int64_t,u_int32_t,u_int64_t>>,
        decltype(cmp)
    > frontier(cmp);

    auto start_time = std::chrono::high_resolution_clock::now();

    // nó inicial
    expansion_counts[start] = 0;
    frontier.push({start, 0, insertion_order++});

    if (start == goal) {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        double seconds = duration.count() / 1'000'000.0;

        std::cout << "Solution found with GBFS:" << std::endl;
        std::cout << n_expanded << "," << 0 << "," << seconds << "," << 0
                  << "," << static_cast<int>(manhattan_distance(start)) << std::endl;
        return true;
    }

    while (!frontier.empty()) {
        auto [current, g, order] = frontier.top();
        frontier.pop();
        n_expanded++;

        for (u_int64_t child : expand(current)) {
            if (expansion_counts.find(child) != expansion_counts.end())
                continue; // já visitado

            expansion_counts[child] = g + 1;

            if (child == goal) {
                auto end_time = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
                double seconds = duration.count() / 1'000'000.0;

                std::cout << "Solution found with GBFS:" << std::endl;
                std::cout << n_expanded << "," << expansion_counts[child] << ","
                          << seconds << "," << 0 << "," << static_cast<int>(manhattan_distance(start)) << std::endl;
                return true;
            }

            frontier.push({child, g + 1, insertion_order++});
        }
    }

    std::cout << "Nenhuma solução encontrada.\n";
    return false;
}

// --- Explicit instantiations (required if you use .cpp for templates) ---
template class Puzzle<BITS_GRID_8>;
template class Puzzle<BITS_GRID_15>;