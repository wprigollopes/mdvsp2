#include "shortest_path.h"
#include <algorithm>
#include <deque>

void TripArcs::build(const IloArray<IloNumArray> &costMatrix, int depots)
{
    num_trips = costMatrix.getSize() - depots;
    adj.assign(num_trips, std::vector<Arc>());

    for (int i = 0; i < num_trips; i++) {
        adj[i].reserve(32);
        for (int j = 0; j < num_trips; j++) {
            if (costMatrix[i + depots][j + depots] != -1)
                adj[i].push_back({j, (long long)costMatrix[i + depots][j + depots]});
        }
    }
}

void SubproblemGraph::build(const IloArray<IloNumArray> &costMatrix,
                            int depots, int depot, const TripArcs &shared)
{
    num_trips = shared.num_trips;
    source = num_trips;
    sink = num_trips + 1;
    trip_arcs = &shared;

    source_arcs.clear();
    source_arcs.reserve(num_trips);
    sink_costs.assign(num_trips, -1);

    for (int i = 0; i < num_trips; i++) {
        if (costMatrix[depot][i + depots] != -1)
            source_arcs.push_back({i, (long long)costMatrix[depot][i + depots]});
        if (costMatrix[i + depots][depot] != -1)
            sink_costs[i] = (long long)costMatrix[i + depots][depot];
    }
}

long long dag_shortest_path(const SubproblemGraph &graph, std::vector<int> &pred,
                            std::vector<long long> &dist,
                            const IloNumArray &pi, const IloNumArray &sigma, int depot)
{
    int N = graph.num_trips;
    int src = graph.source;
    int snk = graph.sink;

    std::fill(dist.begin(), dist.end(), INFEASIBLE_VALUE);
    dist[src] = 0;
    std::fill(pred.begin(), pred.end(), -1);

    // Relax arcs from source (depot -> trips)
    long long sigma_dual = (long long)sigma[depot];
    for (const auto &arc : graph.source_arcs) {
        long long rc = arc.base_cost - sigma_dual;
        if (dist[arc.to] > rc) {
            dist[arc.to] = rc;
            pred[arc.to] = src;
        }
    }

    // Relax arcs from each trip node in topological order (0..N-1)
    for (int i = 0; i < N; i++) {
        if (dist[i] >= INFEASIBLE_VALUE)
            continue;
        long long pi_dual = (long long)pi[i];

        // Trip-to-trip arcs (shared across depots)
        for (const auto &arc : graph.trip_arcs->adj[i]) {
            long long rc = arc.base_cost - pi_dual;
            if (dist[arc.to] > dist[i] + rc) {
                dist[arc.to] = dist[i] + rc;
                pred[arc.to] = i;
            }
        }

        // Trip-to-sink arc (depot-specific)
        if (graph.sink_costs[i] >= 0) {
            long long rc = graph.sink_costs[i] - pi_dual;
            if (dist[snk] > dist[i] + rc) {
                dist[snk] = dist[i] + rc;
                pred[snk] = i;
            }
        }
    }

    return dist[snk];
}

long long spfa_shortest_path(const SubproblemGraph &graph, std::vector<int> &pred,
                             std::vector<long long> &dist,
                             const IloNumArray &pi, const IloNumArray &sigma, int depot)
{
    int N = graph.num_trips;
    int src = graph.source;
    int snk = graph.sink;

    std::fill(dist.begin(), dist.end(), INFEASIBLE_VALUE);
    dist[src] = 0;
    std::fill(pred.begin(), pred.end(), -1);

    std::vector<bool> in_queue(N + 2, false);
    std::deque<int> queue;
    queue.push_back(src);
    in_queue[src] = true;

    while (!queue.empty()) {
        int u = queue.front();
        queue.pop_front();
        in_queue[u] = false;

        // Determine which dual to subtract
        long long dual = 0;
        if (u == src)
            dual = (long long)sigma[depot];
        else if (u < N)
            dual = (long long)pi[u];

        // Edge relaxation with SLF heuristic
        auto relax = [&](int to, long long base_cost) {
            long long rc = base_cost - dual;
            if (dist[to] > dist[u] + rc) {
                dist[to] = dist[u] + rc;
                pred[to] = u;
                if (!in_queue[to]) {
                    if (!queue.empty() && dist[to] < dist[queue.front()])
                        queue.push_front(to);
                    else
                        queue.push_back(to);
                    in_queue[to] = true;
                }
            }
        };

        if (u == src) {
            for (const auto &arc : graph.source_arcs)
                relax(arc.to, arc.base_cost);
        } else if (u < N) {
            for (const auto &arc : graph.trip_arcs->adj[u])
                relax(arc.to, arc.base_cost);
            if (graph.sink_costs[u] >= 0)
                relax(snk, graph.sink_costs[u]);
        }
    }

    return dist[snk];
}
