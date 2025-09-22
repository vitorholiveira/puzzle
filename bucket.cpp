#include "bucket.hpp"

void BucketQueue::push(const AstarNode& node) {
    if (node.f < buckets.size()) {
        // Within each f-bucket, we maintain h-buckets
        if (buckets[node.f].empty()) {
            buckets[node.f].resize(MAX_F_VALUE);
        }
        buckets[node.f][node.h].push(node);
        min_f = std::min(min_f, node.f);
    }
}

AstarNode BucketQueue::pop() {
    // Find minimum f-value
    while (min_f < buckets.size()) {
        if (!buckets[min_f].empty()) {
            // Find minimum h-value within this f-bucket
            for (auto h = 0; h < buckets[min_f].size(); h++) {
                if (!buckets[min_f][h].empty()) {
                    AstarNode node = buckets[min_f][h].top();
                    buckets[min_f][h].pop();
                    
                    // Update min_f if this was the last node with current min_f
                    bool found_min = false;
                    for (auto f = min_f; f < buckets.size() && !found_min; f++) {
                        for (auto h_check = 0; h_check < buckets[f].size(); h_check++) {
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
        min_f++;
    }
    throw std::runtime_error("Queue is empty");
}

bool BucketQueue::empty() const {
    return min_f >= MAX_F_VALUE;
}