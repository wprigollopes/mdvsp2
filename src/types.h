#ifndef TYPES_H
#define TYPES_H

#include <vector>

// Set to 1 when JV (Jonker-Volgenant) LAP solver files are available
#ifndef HAS_JV
#define HAS_JV 0
#endif

constexpr unsigned int ZMIN = 0;
constexpr unsigned int ITER_MIN = 30;
constexpr float OMEGA_MIN = 0.9f;
constexpr long long BIGM = 1000000LL;
constexpr long long INFEASIBLE_VALUE = 100000000LL;
constexpr int CG_TYPE = 1; // 1: with relaxation; 2: without
constexpr double NUM_STABILIZED = 1e6;

struct Arc {
    int to;
    long long base_cost;
};

#endif
