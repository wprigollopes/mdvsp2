#ifndef SHORTEST_PATH_H
#define SHORTEST_PATH_H

#include "types.h"
#include <ilcplex/ilocplex.h>
#include <vector>

// Shared trip-to-trip arcs, built once and referenced by all depot subproblems.
// Avoids rebuilding O(N^2) trip-trip edges D times.
struct TripArcs {
    int num_trips;
    std::vector<std::vector<Arc>> adj; // adj[i] = arcs from trip i to other trips
    void build(const IloArray<IloNumArray> &costMatrix, int depots);
};

// Per-depot subproblem graph. References shared trip-trip arcs and stores
// only depot-specific source/sink arcs.
struct SubproblemGraph {
    int num_trips;
    int source; // = num_trips
    int sink;   // = num_trips + 1
    const TripArcs *trip_arcs; // shared, not owned
    std::vector<Arc> source_arcs;      // depot -> trip arcs
    std::vector<long long> sink_costs; // sink_costs[i] = trip i -> depot cost (-1 if infeasible)

    void build(const IloArray<IloNumArray> &costMatrix, int depots, int depot,
               const TripArcs &shared);
};

// DAG shortest path — O(V+E). Uses pre-allocated dist/pred workspace.
long long dag_shortest_path(const SubproblemGraph &graph, std::vector<int> &pred,
                            std::vector<long long> &dist,
                            const IloNumArray &pi, const IloNumArray &sigma, int depot);

// SPFA shortest path — fallback. Uses pre-allocated dist/pred workspace.
long long spfa_shortest_path(const SubproblemGraph &graph, std::vector<int> &pred,
                             std::vector<long long> &dist,
                             const IloNumArray &pi, const IloNumArray &sigma, int depot);

#endif
