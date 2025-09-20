#include "puzzle.hpp"

// OK
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
std::vector<u_int64_t> Puzzle<BITS_GRID>::expand(const u_int64_t& grid) const {
    std::vector<u_int64_t> children;

    int blank_pos = -1;
    for (u_int8_t pos = 0; pos < max_pos + 1; ++pos) {
        if (((grid >> (pos * 4)) & 0xF) == 0) {
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
        
        u_int64_t new_grid = grid;
        
        // Extract tile to be moved
        u_int64_t tile_value = (grid >> (new_pos * 4)) & 0xF;
        
        // Clear both the new tile position and the blank's old position
        u_int64_t mask = ~(
            (static_cast<u_int64_t>(0xF) << (new_pos * 4)) |
            (static_cast<u_int64_t>(0xF) << (blank_pos * 4))
        );
        new_grid &= mask;
        
        // Place the tile in the blank's old position
        new_grid |= tile_value << (blank_pos * 4);
        
        children.push_back(new_grid);
    }
    return children;
}

// ERROR
template<size_t BITS_GRID>
inline u_int16_t Puzzle<BITS_GRID>::manhattan_distance(const u_int64_t& state) {
    int total = 0;
    for (u_int8_t pos = 0; pos < max_pos; ++pos) {
        u_int8_t tile = (state >> (pos * TILE_BITS)) & 0xF;  // extract 4 bits
        if (tile == 0) continue; // skip blank

        u_int8_t cur_row  = pos / grid_size;
        u_int8_t cur_col  = pos % grid_size;
        u_int8_t goal_row = tile / grid_size;
        u_int8_t goal_col = tile % grid_size;

        total += std::abs(cur_row - goal_row) + std::abs(cur_col - goal_col);
    }
    return total;
}

/*
    PUBLIC METHODS (PUZZLE SOLVERS)
*/
template<size_t BITS_GRID>
bool Puzzle<BITS_GRID>::solve_bfs() {
    std::unordered_map<u_int64_t, u_int32_t> expansion_counts;
    std::queue<u_int64_t> frontier;
    u_int32_t n_expanded = 0;

    u_int64_t start = vector_to_state(states[0]);
    u_int64_t goal  = create_goal_state();

    frontier.push(start);
    expansion_counts[start] = 0;

    auto start_time = std::chrono::high_resolution_clock::now();

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

                std::cout << "Solução encontrada com BFS!\n";
                std::cout << "Número de nós expandidos: " << n_expanded << "\n";
                std::cout << "Comprimento solução ótima: " << expansion_counts[child] << "\n";
                std::cout << "Tempo para a solução: " << seconds << "s\n";
                return true;
            }

            frontier.push(child);
        }
    }

    std::cout << "Nenhuma solução encontrada.\n";
    return false;
}

template<size_t BITS_GRID>
bool Puzzle<BITS_GRID>::solve_idfs() {
    // TODO: Setup nodes

    // TODO: implement actual IDFS search using expand()
    return false;
}

template<size_t BITS_GRID>
bool Puzzle<BITS_GRID>::solve_astar() {
    // TODO: Setup nodes

    // TODO: implement actual ASTAR search using expand()
    return false;
}

template<size_t BITS_GRID>
bool Puzzle<BITS_GRID>::solve_iastar() {
    // TODO: Setup nodes

    // TODO: implement actual IASTAR search using expand()
    return false;
}

template<size_t BITS_GRID>
bool Puzzle<BITS_GRID>::solve_gbfs() {
    // A map to store states we've seen and their path length.
    std::unordered_map<u_int64_t, u_int32_t> expansion_counts;

    // Use a min-priority queue to store states to explore.
    // The pair stores <heuristic_value, state_representation>
    // std::greater<std::pair<u_int32_t, u_int64_t>> makes it a min-heap based on the first element (the heuristic).
    std::priority_queue<std::pair<u_int16_t, u_int64_t>,
                        std::vector<std::pair<u_int16_t, u_int64_t>>,
                        std::greater<std::pair<u_int16_t, u_int64_t>>> frontier;

    u_int32_t n_expanded = 0;

    u_int64_t start = vector_to_state(states[0]);
    u_int64_t goal  = create_goal_state();

    // The heuristic value for the start state
    u_int16_t start_heuristic = manhattan_distance(start);
    u_int64_t heuristic_sum = start_heuristic;

    // Add the starting state to the frontier.
    frontier.push({start_heuristic, start});
    expansion_counts[start] = 0;

    auto start_time = std::chrono::high_resolution_clock::now();

    while (!frontier.empty()) {
        // Get the state with the lowest heuristic value (the "best" one).
        // The 'top()' method returns a reference to the top element.
        u_int64_t current = frontier.top().second;
        // Then, remove it from the frontier.
        frontier.pop();

        // If we've already expanded this state, skip it to avoid redundant work.
        if (expansion_counts.find(current) != expansion_counts.end() && expansion_counts[current] < n_expanded) {
            continue;
        }

        n_expanded++;

        // If we found the goal, we're done.
        if (current == goal) {
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
            double seconds = duration.count() / 1'000'000.0;

            std::cout << "Solução encontrada com GBFS!" << std::endl;
            std::cout << "Número de nós expandidos: " << n_expanded << std::endl;
            std::cout << "Comprimento solução (não garantida como ótima): " << expansion_counts[current] << std::endl;
            std::cout << "Tempo para a solução: " << seconds << "s" << std::endl;
            std::cout << "Valor inicial da heurística: " << start_heuristic << std::endl;
            std::cout << "Valor médio da heurística: " << heuristic_sum / n_expanded << std::endl;
            return true;
        }

        // Expand neighbors
        for (u_int64_t child : expand(current)) {
            // Check if we've seen this child before.
            if (expansion_counts.find(child) != expansion_counts.end()) continue; 

            // If not, calculate its heuristic and add to the frontier.
            u_int16_t child_heuristic = manhattan_distance(child);
            heuristic_sum += child_heuristic;
            frontier.push({child_heuristic, child});
            expansion_counts[child] = expansion_counts[current] + 1;
            // If we found the goal, we're done.
            if (child == goal) {
                auto end_time = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
                double seconds = duration.count() / 1'000'000.0;

                std::cout << "Solução encontrada com GBFS!" << std::endl;
                std::cout << "Número de nós expandidos: " << n_expanded << std::endl;
                std::cout << "Comprimento solução (não garantida como ótima): " << expansion_counts[current] << std::endl;
                std::cout << "Tempo para a solução: " << seconds << "s" << std::endl;
                std::cout << "Valor inicial da heurística: " << start_heuristic << std::endl;
                std::cout << "Valor médio da heurística: " << float(heuristic_sum) / float(n_expanded) << std::endl;
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