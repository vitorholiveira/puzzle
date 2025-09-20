#include "puzzle.hpp"
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
std::vector<u_int8_t> Puzzle<BITS_GRID>::state_to_vector(const u_int64_t& grid) const {
    std::vector<u_int8_t> result(max_pos + 1, 0);

    for (u_int8_t pos = 0; pos < max_pos + 1; ++pos) {
        result[pos] |= (static_cast<u_int64_t>(grid) & (0xF >> (pos * 4)));
        std::cout << result[pos] << std::endl;
    }
    return result;
}

template<size_t BITS_GRID>
void Puzzle<BITS_GRID>::print_state(const u_int64_t& state) const {
    auto vec = state_to_vector(state);
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            std::cout << vec[i * grid_size + j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "--------\n";
}

template<size_t BITS_GRID>
std::vector<u_int64_t> Puzzle<BITS_GRID>::expand(const u_int64_t& grid) const {
    std::vector<u_int64_t> children;

    // Find blank (value = 0)
    int blank_pos = -1;    
    for (u_int8_t pos = 0; pos < max_pos + 1; ++pos) {
        if((0x0 & ((pos * 4) << grid)) == 0x0)
            blank_pos = pos;
    }
    if (blank_pos == -1) return children;

    u_int8_t row = blank_pos / grid_size;
    u_int8_t col = blank_pos % grid_size;

    // Up, Left, Right, Down moves
    std::vector<int> moves = {-grid_size, -1, 1, grid_size};
    std::vector<bool> valid_moves = {
        (row > 0), (col > 0),
        (col < grid_size - 1), (row < grid_size - 1)
    };

    for (int i = 0; i < 4; ++i) {
        if (!valid_moves[i]) continue;
        int new_pos = blank_pos + moves[i];
        if (!(new_pos >= 0 && new_pos <= max_pos)) continue;

        u_int64_t new_grid = grid;

        // Extract tile at new_pos
        int tile_pos = moves[i] + blank_pos;
        u_int64_t tile_value = (grid >> ((tile_pos) * 4)) & 0xF;

        // Clear new tile position
        new_grid &= ~(static_cast<u_int64_t>(0xF) << (tile_pos * 4));
    
        // Place tile in blank position        
        new_grid |= tile_value << (new_pos * 4);
        
        children.push_back(new_grid);
    }
    return children;
}

template<size_t BITS_GRID>
inline int Puzzle<BITS_GRID>::manhattan_distance(const std::bitset<BITS_GRID>& state) {
    u_int64_t raw = state.to_ullong();  // works because 36 <= 64
    int total = 0;

    for (u_int8_t pos = 0; pos < max_pos; ++pos) {
        u_int8_t tile = (raw >> (pos * TILE_BITS)) & 0xF;  // extract 4 bits
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
    std::unordered_map<unsigned long long, u_int32_t> expansion_counts;
    std::queue<unsigned long long> frontier;
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

        std::cout << "Exploring node:\n";
        print_state(current);

        // Expand neighbors
        for (u_int64_t child : expand(current)) {
            std::cout << "Expanding to:\n";
            print_state(child);
            // If we already saw this child, skip
            if (expansion_counts.find(child) != expansion_counts.end())
                continue;

            expansion_counts[child] = expansion_counts[current] + 1;

            if (child == goal) {
                std::cout << "Solução encontrada com BFS!\n";
                std::cout << "Número de nós expandidos: " << n_expanded << "\n";
                std::cout << "Comprimento solução ótima: " << expansion_counts[child] << "\n";

                auto end_time = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
                double seconds = duration.count() / 1'000'000.0;
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
    // TODO: Setup nodes

    // TODO: implement actual GBFS search using expand()
    return false;
}

// --- Explicit instantiations (required if you use .cpp for templates) ---
template class Puzzle<BITS_GRID_8>;
template class Puzzle<BITS_GRID_15>;