#ifndef SHORTEST_PATH_H
#define SHORTEST_PATH_H

#include "types.h"
#include <ilcplex/ilocplex.h>
#include <vector>

struct SubproblemGraph {
    int num_trips;
    int source; // = num_trips
    int sink;   // = num_trips + 1
    std::vector<std::vector<Arc>> adj;

    // Build graph topology and base costs from cost matrix for a given depot.
    // Call once after reductions; reuse across CG iterations.
    void build(const IloArray<IloNumArray> &costMatrix, int depots, int depot);
};

// DAG shortest path: O(V+E). Requires trips sorted by start time (DAG property).
// Returns distance to sink; negative means an improving column was found.
// Predecessor array is filled for path reconstruction.
long long dag_shortest_path(const SubproblemGraph &graph,
                            std::vector<int> &pred,
                            const IloNumArray &pi,
                            const IloNumArray &sigma,
                            int depot);

// SPFA (Shortest Path Faster Algorithm) with SLF heuristic.
// Validation/fallback: works on any graph, no DAG assumption.
long long spfa_shortest_path(const SubproblemGraph &graph,
                             std::vector<int> &pred,
                             const IloNumArray &pi,
                             const IloNumArray &sigma,
                             int depot);

#endif
