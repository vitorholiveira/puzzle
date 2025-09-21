#include "puzzle.hpp"

std::vector<std::vector<int>> read_states(int argc, char* argv[]) {
    std::vector<std::vector<int>> states;
    std::vector<int> first_vector;
    first_vector.reserve(16);
    states.push_back(first_vector);
    int state_index = 0;
    std::cout << "Puzzle state " << state_index + 1 << ": ";
    for(int i = 2; i < argc; i++) {
        std::string v = std::string(argv[i]);
        if(v.length() > 1 && (v[1] == ',' || v[2] == ',')){
            std::vector<int> curr_vector;
            curr_vector.reserve(16);
            states.push_back(curr_vector);
            v = (v[1] == ',') ? v.substr(0,1) : v.substr(0,2);
            states[state_index].push_back(stoi(v));
            state_index++;
            std::cout << v << std::endl;
            std::cout << "Puzzle state " << state_index + 1 << ": ";
        } else {
            std::cout << v << " ";
            states[state_index].push_back(stoi(v));
        }
    }
    std::cout << std::endl;
    return states;
}

int main(int argc, char* argv[]) {
    std::string usage_msg = "Usage: ./main {-bfs | -idfs | -astar | -idastar | -gbfs} list_numbers"; 
    if(argc < 2) {
        std::cerr << usage_msg << std::endl;
        return 1;
    }
    std::string algorithm = std::string(argv[1]);
    if(algorithm != BFS   && algorithm != IDFS    &&
       algorithm != ASTAR && algorithm != IDASTAR &&
       algorithm != GBFS) {
        std::cerr << "[ERROR] -> " << usage_msg << std::endl;
        return 2;
    }
    std::cout << "algorithm: " << algorithm << std::endl;

    int len_list_numbers = argc - 2;
    if(((len_list_numbers % 9) != 0 ) && (((len_list_numbers % 16) != 0))){
        std::cerr << "[ERROR] -> " << "You should pass a state of a 8-puzzle or 15-puzzle" << std::endl;
        std::cerr << "[ERROR] -> " << "The list should have a multiple of 9 or a multiple of 16 elements" << std::endl;
        std::cerr << "[ERROR] -> " << usage_msg << std::endl;
    }

    
    std::vector<std::vector<int>> states = read_states(argc, argv);
    int type = states[0].size() - 1;
    if(type == 8) {
        Puzzle<BITS_GRID_8> puzzle(states);
        puzzle.solve(algorithm);
    } else if(type == 15) {
        Puzzle<BITS_GRID_15> puzzle(states);
        if (algorithm == ASTAR) {
            puzzle.solve(algorithm);
        }
    } else {
        std::cerr << "[ERROR] -> " << "There is no implementation to solve the " << type << "-puzzle with the " << algorithm << " algorithm." << std::endl;
    }
    return 0;
}