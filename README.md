# 8-Puzzle and 15-Puzzle Solvers

This project implements and analyzes the performance of several classic search algorithms applied to solving the 8-Puzzle and 15-Puzzle. The analysis focuses on metrics such as expanded nodes, solution length, execution time, and success rate.


## Compile and Run

1.  Compile the program by running the `make` command, which creates the `main` executable.
    ```bash
    make
    ```
2.  Run the solver on a single puzzle using the following format.
    ```bash
    ./main {algorithm} {state}
    ```
      * **`{algorithm}`**: Choose from `-bfs`, `-idfs`, `-astar`, `-idastar`, or `-gbfs`.
      * **`{state}`**: A space-separated list of numbers for the puzzle. The list must have 9 numbers for an 8-puzzle or 16 for a 15-puzzle.

### **Examples**

  * 8-Puzzle:
    ```bash
    ./main -astar 1 2 3 4 5 6 7 8 0
    ```
  * 15-Puzzle:
    ```bash
    ./main -astar 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 0
    ```

## Batch Testing

To run tests on multiple puzzle instances from a file, use the `run.sh` script.

1.  **Configure the script** by editing the variables at the top of `run.sh` to set the algorithm, input/output files, and resource limits (memory and time).
2.  **Execute the script**:
    ```bash
    ./run.sh
    ```

The script will process each puzzle instance, enforce the limits, and save the results to your specified output file.


## Results

### Test Environment

The hardware used for testing was a Dell Inspiron 14 5440 laptop, with the following specifications:

* **Processor:** Intel Core 7 150U × 12
* **RAM:** 32 GB
* **Operating System:** Ubuntu 24.04.3 LTS

### Comparative Performance Table

The table below summarizes the average results obtained for each algorithm in a battery of tests.

| Algorithm | Average Expanded Nodes | Average Solution Length | Average Time (s) | Average Heuristic | Average Initial Heuristic | Success Rate |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `bfs8` | 81,459.54 | 22.16 | 0.022 | 0.00 | 13.88 | 100% |
| `idfs8` | 2,578,290.56 | 22.16 | 0.186 | 0.00 | 13.88 | 100% |
| `gbfs8` | 392.42 | 140.52 | 0.0004 | 6.89 | 13.88 | 100% |
| `astar8` | 895.00 | 22.16 | 0.0002 | 10.03 | 13.88 | 100% |
| `idastar8` | 2,373.03 | 22.16 | 0.0006 | 10.43 | 13.88 | 100% |
| `astar15` | 5,643,956.59 | 52.14 | 3.572 | 25.63 | 36.70 | 90% |
| `idastar15`| 46,561,394.96 | 51.89 | 4.691 | 15.11 | 36.58 | 85% |

### Failure Analysis in the 15-Puzzle

All failures for both algorithms occurred due to the execution time limit.

* **A* Algorithm (`astar15`):** Failed in 10 instances of the 15-Puzzle. The unsolved instances have the following indices: 17, 49, 53, 56, 59, 60, 66, 82, 88, and 92.
* **IDA* Algorithm (`idastar15`):** Failed to solve 15 instances of the 15-Puzzle, which include the following indices: 14, 17, 22, 32, 49, 56, 59, 60, 63, 66, 72, 82, 88, 91, and 92.

## Implementation Details

### State Representation and Base Structure

The puzzle states were represented by a 64-bit unsigned integer (`u_int_64`), where each piece occupies 4 bits. This choice was made because it is more optimized for memory access than more complex structures like `struct`s and vectors, allowing for faster expansion and comparison operations. Bit manipulation requires shift and mask operations, but it eliminates the overhead of pointers and indirect access.

### Heuristic Function and Auxiliary Structures

The Manhattan distance heuristic was optimized using a pre-computed table (`std::vector<std::vector<int>>`), reducing runtime calculations to simple indexed lookups. This choice represents a trade-off: increased memory consumption in exchange for performance gains during node expansion.

### Algorithms

#### Breadth-First Search (BFS)

BFS was implemented with a queue (`std::queue`) to manage the state frontier and a hash table (`std::unordered_map`) to record visited states along with their depths.

#### Iterative Deepening Depth-First Search (IDFS)

For IDFS, a recursive approach was used through the `recursive_dls` function. Each call represents the expansion of a state up to the current depth limit.

#### A*

A* combines accumulated cost and heuristics. The chosen structure to represent the frontier was a custom priority queue (`BucketQueue`), designed to reduce the insertion and removal cost from O(log n) (with `std::priority_queue`) to O(1). The set of expanded states was stored in a hash table (`std::unordered_set`). The decision to use `BucketQueue` is based on A*'s main weakness, which is the large number of memory accesses.

#### Greedy Best-First Search (GBFS)

In GBFS, the frontier was maintained in a priority queue (`std::priority_queue`), which greedily evaluates the next state to expand based only on the lowest heuristic value. A hash set (`std::unordered_set`) was used for visited nodes.

#### Iterative Deepening A* (IDA*)

IDA* was implemented recursively, similar to IDFS, but with the difference that the depth limit is defined by the sum of the cost *g* and the heuristic *h* (f = g + h). This choice eliminates the need for priority queues, keeping memory consumption extremely low. The main challenge lies in the recursion's return logic, which needs to signal search cutoffs and dynamically adjust the limit for the next iteration.