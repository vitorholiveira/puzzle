
#ifndef BUCKETQUEUE_HPP
#define BUCKETQUEUE_HPP

#include <string>
#include <vector>
#include <queue>
#include <chrono>
#include <stack>
#include <iostream>

static constexpr int MAX_F_VALUE = 200; // Adjust based on expected max f-value

struct AstarNode {
    u_int64_t state;
    int g;  // cost from start
    int h;  // heuristic cost to goal
    int f;  // g + h
    u_int64_t parent;
    int move; // move that led to this state

    AstarNode() : state(0), g(0), h(0), f(0), parent(0), move(0) {} // novo
    AstarNode(u_int64_t s, int g_val, int h_val, u_int64_t p = 0, int m = -1) 
        : state(s), g(g_val), h(h_val), f(g_val + h_val), parent(p), move(m) {}
};

struct AstarBucketQueue {
    std::vector<std::vector<std::stack<AstarNode>>> buckets;
    int min_f;
    
    AstarBucketQueue() : buckets(MAX_F_VALUE), min_f(MAX_F_VALUE) {}
    void push(const AstarNode& node);
    AstarNode pop();
    bool empty() const;
};

#endif // BUCKETQUEUE_HPP