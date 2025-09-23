#include "puzzle.hpp"
#include "bucket.hpp"
#define MAX_INT 999999999

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
            solve_idastar(start);
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
        long double seconds = (long double)duration.count() / 1'000'000.0;
        int h_start = manhattan_distance(start);

        SearchStatistics(n_expanded, 0, seconds, 0, h_start);
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
                long double seconds = (long double) duration.count() / 1'000'000.0;
                int h_start = manhattan_distance(start);

                SearchStatistics(n_expanded, visited[child], seconds, 0, h_start);
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
        long double seconds = (long double)duration.count() / 1'000'000.0;
        
        SearchStatistics(n_expanded, 0, seconds, 0, manhattan_distance(start));
        return true;
    }

    for (u_int32_t depth_limit = 1; depth_limit <= 50; depth_limit++) {
        u_int32_t solution_depth = 0;
        
        int result = recursive_dls(start, 0, depth_limit, 0, n_expanded, solution_depth);
        
        if (result == 1) {
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
            long double seconds = (long double)duration.count() / 1'000'000.0;
            
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
    using Node = std::tuple<u_int64_t, u_int32_t, u_int64_t, u_int32_t, u_int64_t>; // state, g, insertion order, h, parent
    std::unordered_set<u_int64_t> visited;                     // expanded states
    u_int32_t n_expanded = 0;
    u_int64_t insertion_order = 0;
    u_int32_t g_new;

    BucketQueue<Node> frontier;
    auto start_time = std::chrono::high_resolution_clock::now();

    // initialize
    u_int32_t start_h = (u_int32_t)manhattan_distance(start);
    frontier.push({start, 0, insertion_order++, start_h, 0}, start_h, start_h);
    int heuristic_calls = 1;
    int heuristic_sum = manhattan_distance(start);

    while (!frontier.empty()) {
        auto [current, g, order, h, parent] = frontier.pop();

        if (visited.count(current)) continue;
        visited.insert(current);

        if (current == goal) {
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
            long double seconds = (long double)duration.count() / 1'000'000.0;

            long double avg_h = (long double)heuristic_sum / (long double)heuristic_calls;
            u_int32_t h_start = static_cast<u_int32_t>(manhattan_distance(start));

            SearchStatistics(n_expanded, g, seconds, avg_h, h_start);
            return true;
        }

        n_expanded++;
        g_new = g + 1;

        for (const auto& child : expand(current)) {
            if (child != parent) {
                u_int32_t h = (u_int32_t)manhattan_distance(child);
                frontier.push({child, g_new, insertion_order++, h, current},g_new+h,h);
                heuristic_calls++;
                heuristic_sum += manhattan_distance(child);
            }
        }
    }
    
    return false; // No solution found
}

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
            long double seconds = (long double)duration.count() / 1'000'000.0;
            u_int64_t h_start = static_cast<u_int64_t>(manhattan_distance(start));
            long double avg_h = (long double)heuristic_sum / (long double)heuristic_calls;
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

bool Puzzle::solve_idastar(const u_int64_t& start) {
    u_int32_t n_expanded = 0;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    int depth_limit = manhattan_distance(start);
    int heuristic_calls = 1;
    int heuristic_sum = manhattan_distance(start);

    while (true) {
        u_int32_t solution_depth = 0;

        int result = recursive_adls(start, 0, depth_limit, 0, n_expanded, solution_depth, heuristic_calls, heuristic_sum);
        
        if (result == -1) {
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
            long double seconds = (long double)duration.count() / 1'000'000.0;
            long double avg_h = (long double)heuristic_sum / (long double)heuristic_calls;

            SearchStatistics(n_expanded, solution_depth, seconds, avg_h, manhattan_distance(start));
            return true;
        }
        
        if (result == MAX_INT) {
            break;
        }

        depth_limit = result;
    }
    
    return false;
}

int Puzzle::recursive_adls(
    const u_int64_t& current_state, 
    const u_int64_t& parent_state,
    u_int32_t limit, 
    u_int32_t current_depth,
    u_int32_t& n_expanded,
    u_int32_t& solution_depth,
    int& heuristic_calls,
    int& heuristic_sum) {

    int heuristic = manhattan_distance(current_state);
    heuristic_calls++;
    heuristic_sum += heuristic;

    u_int32_t f_cost = current_depth + heuristic;
    
    if (f_cost > limit) {
        return f_cost;
    }
    
    if (current_state == goal) {
        solution_depth = current_depth;
        return -1;
    }
    
    n_expanded++;
    int next_limit = MAX_INT;
    
    for (u_int64_t child : expand(current_state)) {
        if (child == parent_state) {
            continue;
        }

        int result = recursive_adls(child, current_state, limit, current_depth + 1, n_expanded, solution_depth, heuristic_calls, heuristic_sum);

        if (result == -1) {
            return -1;
        } else {
            next_limit = std::min(next_limit, result);
        }
    }
    
    return next_limit;
}
