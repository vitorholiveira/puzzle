#include "puzzle.hpp"

/*
    CREATING GOAL STATES
*/
template<size_t BITS_GRID>
std::bitset<BITS_GRID> Puzzle<BITS_GRID>::create_goal_state() const {
    std::bitset<BITS_GRID> goal;

    for (int pos = 0; pos < max_pos; ++pos) {
        int value = pos + 1;
        for (int bit = 0; bit < TILE_BITS; ++bit) {
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

    for (int i = 0; i <= max_pos; ++i) {
        for (int bit = 0; bit < TILE_BITS; ++bit) { // fixed <= to <
            if (grid_vec[i] & (1 << bit)) {
                result.set(i * TILE_BITS + bit);
            }
        }
    }
    return result;
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

    // Up, Down, Left, Right moves
    std::vector<int> moves = {-grid_size, grid_size, -1, 1};
    std::vector<bool> valid_moves = {
        (row > 0), (row < grid_size - 1),
        (col > 0), (col < grid_size - 1)
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
    PUBLIC METHODS (example: BFS)
*/
template<size_t BITS_GRID>
bool Puzzle<BITS_GRID>::solve_bfs() {
    std::unordered_set<std::bitset<BITS_GRID>> visited;

    for (const auto& state : states) {
        auto bits = vector_to_bitset(state);
        visited.insert(bits);
    }

    // TODO: implement actual BFS search using expand()
    return false;
}

// --- Explicit instantiations (required if you use .cpp for templates) ---
template class Puzzle<BITS_GRID_8>;
template class Puzzle<BITS_GRID_15>;
