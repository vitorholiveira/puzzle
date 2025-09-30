/*
 * bucket.hpp
 *
 * Defines the BucketQueue data structure, a custom priority queue
 * based on "buckets" to organize elements by cost/heuristic.
 * This implementation offers more efficient insertions and removals (close to O(1))
 * compared to std::priority_queue, being especially useful in search algorithms
 * like A*.
 */
#ifndef BUCKETQUEUE_HPP
#define BUCKETQUEUE_HPP

#include <vector>
#include <deque>
#include <stdexcept>
#include <climits>

inline constexpr std::size_t MAX_F_VALUE = 200;

template <typename T>
struct BucketQueue {
private:
    // Use deque instead of stack for better memory layout and cache performance
    using Container = std::deque<T>;
    
    // Pre-allocate all buckets to avoid dynamic resizing
    std::vector<std::vector<Container>> buckets;
    std::size_t min_f;
    std::size_t min_h_for_f;
    
    // Bit vector to quickly check if f-bucket has any elements
    std::vector<bool> f_has_elements;
    // Bit vector to quickly check if specific (f,h) bucket has elements
    std::vector<std::vector<bool>> fh_has_elements;

public:
    BucketQueue() : buckets(MAX_F_VALUE), min_f(MAX_F_VALUE), min_h_for_f(MAX_F_VALUE),
                    f_has_elements(MAX_F_VALUE, false), fh_has_elements(MAX_F_VALUE) {
        // Pre-allocate all inner vectors to avoid dynamic allocation during push
        for (std::size_t f = 0; f < MAX_F_VALUE; ++f) {
            buckets[f].resize(MAX_F_VALUE);
            fh_has_elements[f].resize(MAX_F_VALUE, false);
        }
    }

    inline void push(const T& node, std::size_t f, std::size_t h) noexcept {
        // Bounds checking removed for performance - caller must ensure valid indices
        buckets[f][h].push_back(node);
        
        // Update bit vectors
        f_has_elements[f] = true;
        fh_has_elements[f][h] = true;
        
        // Update minimums with branchless comparison
        if (f < min_f || (f == min_f && h < min_h_for_f)) {
            min_f = f;
            min_h_for_f = h;
        }
    }

    inline T pop() {
        // Find next valid bucket using bit vectors for fast lookup
        while (min_f < MAX_F_VALUE) {
            if (f_has_elements[min_f]) {
                // Find minimum h in current f bucket
                while (min_h_for_f < MAX_F_VALUE && !fh_has_elements[min_f][min_h_for_f]) {
                    ++min_h_for_f;
                }
                
                if (min_h_for_f < MAX_F_VALUE) {
                    T node = std::move(buckets[min_f][min_h_for_f].back());
                    buckets[min_f][min_h_for_f].pop_back();
                    
                    // Update bit vectors if bucket became empty
                    if (buckets[min_f][min_h_for_f].empty()) {
                        fh_has_elements[min_f][min_h_for_f] = false;
                        
                        // Check if entire f-bucket is now empty
                        bool f_empty = true;
                        for (std::size_t h = 0; h < MAX_F_VALUE; ++h) {
                            if (fh_has_elements[min_f][h]) {
                                f_empty = false;
                                break;
                            }
                        }
                        if (f_empty) {
                            f_has_elements[min_f] = false;
                        }
                        
                        // Find next non-empty h in current f
                        ++min_h_for_f;
                    }
                    
                    return node;
                }
            }
            
            // Move to next f bucket
            ++min_f;
            min_h_for_f = 0;
        }
        
        throw std::runtime_error("Queue is empty");
    }

    inline bool empty() const noexcept {
        return min_f >= MAX_F_VALUE;
    }
    
    // Additional optimized methods
    inline void reserve_capacity(std::size_t capacity) {
        for (std::size_t f = 0; f < MAX_F_VALUE; ++f) {
            for (std::size_t h = 0; h < MAX_F_VALUE; ++h) {
                buckets[f][h].reserve(capacity / (MAX_F_VALUE * MAX_F_VALUE) + 1);
            }
        }
    }
    
    inline std::size_t size() const noexcept {
        std::size_t total = 0;
        for (std::size_t f = 0; f < MAX_F_VALUE; ++f) {
            if (f_has_elements[f]) {
                for (std::size_t h = 0; h < MAX_F_VALUE; ++h) {
                    if (fh_has_elements[f][h]) {
                        total += buckets[f][h].size();
                    }
                }
            }
        }
        return total;
    }
};

#endif // BUCKETQUEUE_HPP
