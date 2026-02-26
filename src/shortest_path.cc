#include "shortest_path.h"
#include <deque>

void SubproblemGraph::build(const IloArray<IloNumArray> &costMatrix,
                            int depots, int depot)
{
    num_trips = costMatrix.getSize() - depots;
    source = num_trips;
    sink = num_trips + 1;
    adj.assign(num_trips + 2, std::vector<Arc>());

    // Reserve typical capacity to avoid repeated allocations
    adj[source].reserve(num_trips);
    for (int i = 0; i < num_trips; i++)
        adj[i].reserve(32);

    for (int i = 0; i < num_trips; i++) {
        // Trip-to-trip arcs
        for (int j = 0; j < num_trips; j++) {
            if (costMatrix[i + depots][j + depots] == -1)
                continue;
            adj[i].push_back({j, (long long)costMatrix[i + depots][j + depots]});
        }

        // Source (depot) -> trip i
        if (costMatrix[depot][i + depots] != -1) {
            adj[source].push_back({i, (long long)costMatrix[depot][i + depots]});
        }

        // Trip i -> sink (back to depot)
        if (costMatrix[i + depots][depot] != -1) {
            adj[i].push_back({sink, (long long)costMatrix[i + depots][depot]});
        }
    }
}

long long dag_shortest_path(const SubproblemGraph &graph,
                            std::vector<int> &pred,
                            const IloNumArray &pi,
                            const IloNumArray &sigma,
                            int depot)
{
    int N = graph.num_trips;
    int src = graph.source;
    int snk = graph.sink;

    std::vector<long long> dist(N + 2, INFEASIBLE_VALUE);
    dist[src] = 0;
    pred.assign(N + 2, -1);

    // Relax arcs from source (depot -> trips)
    long long sigma_dual = (long long)sigma[depot];
    for (const auto &arc : graph.adj[src]) {
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
        for (const auto &arc : graph.adj[i]) {
            long long rc = arc.base_cost - pi_dual;
            if (dist[arc.to] > dist[i] + rc) {
                dist[arc.to] = dist[i] + rc;
                pred[arc.to] = i;
            }
        }
    }

    return dist[snk];
}

long long spfa_shortest_path(const SubproblemGraph &graph,
                             std::vector<int> &pred,
                             const IloNumArray &pi,
                             const IloNumArray &sigma,
                             int depot)
{
    int N = graph.num_trips;
    int src = graph.source;
    int snk = graph.sink;

    std::vector<long long> dist(N + 2, INFEASIBLE_VALUE);
    dist[src] = 0;
    pred.assign(N + 2, -1);

    std::vector<bool> in_queue(N + 2, false);
    std::deque<int> queue;
    queue.push_back(src);
    in_queue[src] = true;

    while (!queue.empty()) {
        int u = queue.front();
        queue.pop_front();
        in_queue[u] = false;

        // Determine which dual to subtract for arcs leaving u
        long long dual = 0;
        if (u == src)
            dual = (long long)sigma[depot];
        else if (u < N)
            dual = (long long)pi[u];

        for (const auto &arc : graph.adj[u]) {
            long long rc = arc.base_cost - dual;
            if (dist[arc.to] > dist[u] + rc) {
                dist[arc.to] = dist[u] + rc;
                pred[arc.to] = u;
                if (!in_queue[arc.to]) {
                    // SLF heuristic: push to front if label < front label
                    if (!queue.empty() && dist[arc.to] < dist[queue.front()])
                        queue.push_front(arc.to);
                    else
                        queue.push_back(arc.to);
                    in_queue[arc.to] = true;
                }
            }
        }
    }

    return dist[snk];
}
