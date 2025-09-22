#ifndef BUCKETQUEUE_HPP
#define BUCKETQUEUE_HPP

#include <vector>
#include <stack>
#include <stdexcept>
#include <algorithm>

inline constexpr std::size_t MAX_F_VALUE = 200;

template <typename T>
struct BucketQueue {
    std::vector<std::vector<std::stack<T>>> buckets;
    std::size_t min_f;

    BucketQueue() : buckets(MAX_F_VALUE), min_f(MAX_F_VALUE) {}

    void push(const T& node);
    T pop();
    bool empty() const;
};

template <typename T>
void BucketQueue<T>::push(const T& node) {
    if (node.f < buckets.size()) {
        if (buckets[node.f].empty()) {
            buckets[node.f].resize(MAX_F_VALUE);
        }
        buckets[node.f][node.h].push(node);
        min_f = std::min(min_f, static_cast<std::size_t>(node.f));
    }
}

template <typename T>
T BucketQueue<T>::pop() {
    while (min_f < buckets.size()) {
        if (!buckets[min_f].empty()) {
            for (std::size_t h = 0; h < buckets[min_f].size(); ++h) {
                if (!buckets[min_f][h].empty()) {
                    T node = buckets[min_f][h].top();
                    buckets[min_f][h].pop();

                    bool found_min = false;
                    for (std::size_t f = min_f; f < buckets.size() && !found_min; ++f) {
                        for (std::size_t h_check = 0; h_check < buckets[f].size(); ++h_check) {
                            if (!buckets[f][h_check].empty()) {
                                min_f = f;
                                found_min = true;
                                break;
                            }
                        }
                    }
                    if (!found_min) min_f = MAX_F_VALUE;
                    return node;
                }
            }
        }
        ++min_f;
    }
    throw std::runtime_error("Queue is empty");
}

template <typename T>
bool BucketQueue<T>::empty() const {
    return min_f >= MAX_F_VALUE;
}

#endif // BUCKETQUEUE_HPP
