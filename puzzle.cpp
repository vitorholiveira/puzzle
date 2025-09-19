#include "puzzle.hpp"

/*
    CREATING GOAL STATES
*/
template<size_t BITS_GRID>
std::bitset<BITS_GRID> Puzzle<BITS_GRID>::create_goal_state() const {
    std::bitset<BITS_GRID> goal;

    for (int pos = 0; pos <= max_pos; pos++) {
        int value = pos;
        for (int bit = 0; bit < TILE_BITS; bit++) {
            if (value & (1 << bit)) {
                goal.set(pos * TILE_BITS + bit);
            }
        }
    }
    // last tile is blank = 0 (already default)

    return goal;
}

/*
    PRIVATE METHODS
*/
template<size_t BITS_GRID>
std::bitset<BITS_GRID> Puzzle<BITS_GRID>::vector_to_bitset(const std::vector<int>& grid_vec) const {
    std::bitset<BITS_GRID> result;

    for (size_t i = 0; i <= max_pos; i++) {
        for (int bit = 0; bit < TILE_BITS; bit++) {
            if (grid_vec[i] & (1 << bit)) {
                result.set(i * TILE_BITS + bit);
            }
        }
    }
    return result;
}

template<size_t BITS_GRID>
std::vector<int> Puzzle<BITS_GRID>::bitset_to_vector(const std::bitset<BITS_GRID>& grid) const {
    std::vector<int> result(max_pos + 1, 0);

    for (int i = 0; i < max_pos + 1; ++i) {
        int value = 0;
        for (int bit = 0; bit < TILE_BITS; ++bit) {
            if (grid[i * TILE_BITS + bit]) {
                value |= (1 << bit);
            }
        }
        result[i] = value;
    }
    return result;
}

template<size_t BITS_GRID>
void Puzzle<BITS_GRID>::print_state(const std::bitset<BITS_GRID>& grid) const {
    auto vec = bitset_to_vector(grid);
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            std::cout << vec[i * grid_size + j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "--------\n";
}

template<size_t BITS_GRID>
std::vector<std::bitset<BITS_GRID>> Puzzle<BITS_GRID>::expand(const std::bitset<BITS_GRID>& grid) const {
    std::vector<std::bitset<BITS_GRID>> children;

    // Find blank (value = 0)
    int blank_pos = -1;
    for (int i = 0; i <= max_pos; ++i) {
        bool is_blank = true;
        for (int bit = 0; bit < TILE_BITS; ++bit) {
            if (grid[i * TILE_BITS + bit]) {
                is_blank = false;
                break;
            }
        }
        if (is_blank) {
            blank_pos = i;
            break;
        }
    }
    if (blank_pos == -1) return children;

    int row = blank_pos / grid_size;
    int col = blank_pos % grid_size;

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

        std::bitset<BITS_GRID> new_grid = grid;

        // Extract tile at new_pos
        std::bitset<TILE_BITS> tile_value;
        for (int bit = 0; bit < TILE_BITS; ++bit) {
            tile_value[bit] = new_grid[new_pos * TILE_BITS + bit];
        }

        // Clear both positions
        for (int bit = 0; bit < TILE_BITS; ++bit) {
            new_grid[blank_pos * TILE_BITS + bit] = 0;
            new_grid[new_pos * TILE_BITS + bit] = 0;
        }

        // Place tile in blank position
        for (int bit = 0; bit < TILE_BITS; ++bit) {
            new_grid[blank_pos * TILE_BITS + bit] = tile_value[bit];
        }

        children.push_back(new_grid);
    }
    return children;
}

/*
    PUBLIC METHODS (PUZZLE SOLVERS)
*/
template<size_t BITS_GRID>
bool Puzzle<BITS_GRID>::solve_bfs() {
    // Setup nodes
    std::unordered_set<std::bitset<BITS_GRID>> visited;
    std::queue<std::bitset<BITS_GRID>> frontier;
    int n_expanded = 0;

    auto start = vector_to_bitset(states[0]);
    auto goal = create_goal_state();

    frontier.push(start);

    while (!frontier.empty()) {
        auto current = frontier.front();
        frontier.pop();

        if (visited.find(current) != visited.end()) {
            continue;
        }
        std::cout << "Exploring node:\n";
        print_state(current);

        visited.insert(current);

        // Expande vizinhos
        n_expanded++;
        for (auto& child : expand(current)) {
            std::cout << "Expanding to:\n";
            print_state(child);

            if (child == goal) {
                std::cout << "Solução encontrada com BFS!" << std::endl;
                std::cout << "Número de nós expandidos: " << n_expanded << std::endl;
                return true;
            }

            if (child != current && visited.find(child) == visited.end()) {
                frontier.push(child);
            }
        }
    }

    std::cout << "Nenhuma solução encontrada." << std::endl;
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