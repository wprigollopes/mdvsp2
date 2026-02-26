#include "reduction.h"
#include <algorithm>
#include <cmath>
#include <utility>

using namespace std;

// R1: Keep sqrt(N) cheapest arcs per row and per column.
// Rewritten with nth_element: O(N^2) instead of original O(N^2 * sqrt(N)).
// Writes directly to flat assignMatrix (no intermediate reduceItems matrix).
void reduction1(const IloArray<IloNumArray> &costMatrix, int depots, int nodes,
                int matrixSize, vector<char> &assignMatrix)
{
    unsigned int keep = ceil(sqrt((float)nodes));
    vector<pair<double, int>> costs;
    costs.reserve(nodes);

    // Row reduction: for each row, keep the cheapest sqrt(N) arcs
    for (int i = depots; i < matrixSize; i++) {
        costs.clear();
        for (int j = depots; j < matrixSize; j++)
            if (costMatrix[i][j] != -1)
                costs.push_back({costMatrix[i][j], j});

        unsigned int n = min((unsigned int)costs.size(), keep);
        if (n < costs.size())
            nth_element(costs.begin(), costs.begin() + n, costs.end());

        for (unsigned int x = 0; x < n; x++)
            assignMatrix[i * matrixSize + costs[x].second] = 1;
    }

    // Column reduction: for each column, keep the cheapest sqrt(N) arcs
    for (int j = depots; j < matrixSize; j++) {
        costs.clear();
        for (int i = depots; i < matrixSize; i++)
            if (costMatrix[i][j] != -1)
                costs.push_back({costMatrix[i][j], i});

        unsigned int n = min((unsigned int)costs.size(), keep);
        if (n < costs.size())
            nth_element(costs.begin(), costs.begin() + n, costs.end());

        for (unsigned int x = 0; x < n; x++)
            assignMatrix[costs[x].second * matrixSize + j] = 1;
    }
}

#if HAS_JV

// R2: LAP-based reduction per depot using Jonker-Volgenant.
// Parallelized with OpenMP — each depot is independent.
// Each thread allocates its own LAP work arrays to avoid contention.
void reduction2(const IloArray<IloNumArray> &costMatrix, int depots, int nodes,
                int matrixSize, vector<char> &assignMatrix,
                vector<vector<vector<bool>>> &assignSImatrix)
{
    int enddim = nodes + nodes;

    #pragma omp parallel for schedule(dynamic)
    for (int depot = 0; depot < depots; depot++) {
        // Per-thread LAP work arrays
        vector<vector<cost>> origAssignCost(enddim, vector<cost>(enddim));
        vector<col> rowsol(enddim);
        vector<row> colsol(enddim);
        vector<cost> u(enddim);
        vector<cost> v(enddim);

        // Zero the bottom-right quadrant (same for all depots)
        for (int i = 0; i < nodes; i++)
            for (int j = 0; j < nodes; j++)
                origAssignCost[i + nodes][j + nodes] = 0;

        // Build pointer array for lap() which expects cost**
        vector<cost *> costPtrs(enddim);
        for (int i = 0; i < enddim; i++)
            costPtrs[i] = origAssignCost[i].data();

        // Fill cost matrix for this depot
        for (int i = 0; i < nodes; i++) {
            for (int j = 0; j < nodes; j++) {
                origAssignCost[i][j] = (costMatrix[i + depots][j + depots] == -1)
                    ? INFEASIBLE_VALUE : (cost)costMatrix[i + depots][j + depots];
                origAssignCost[i + nodes][j] = (costMatrix[depot][j + depots] == -1)
                    ? INFEASIBLE_VALUE : (cost)costMatrix[depot][j + depots];
                origAssignCost[j][i + nodes] = (costMatrix[j + depots][depot] == -1)
                    ? INFEASIBLE_VALUE : (cost)costMatrix[j + depots][depot];
            }
        }

        lap(enddim, costPtrs.data(), rowsol.data(), colsol.data(),
            u.data(), v.data());

        // Write results — assignMatrix uses char so distinct indices are safe,
        // but assignSImatrix[depot] is depot-exclusive so no contention.
        for (int dim = 0; dim < enddim; dim++) {
            int from = dim >= nodes ? depot : dim + depots;
            int to = rowsol[dim] >= nodes ? depot : rowsol[dim] + depots;
            assignSImatrix[depot][from][to] = true;
            assignMatrix[from * matrixSize + to] = 1;
        }
        assignSImatrix[depot][depot][depot] = false;
        assignMatrix[depot * matrixSize + depot] = 0;
    }
}

#endif // HAS_JV

// R3: LP relaxation-based reduction.
// Only creates variables for feasible arcs (costMatrix != -1), skipping
// depot-to-depot arcs. This shrinks the LP and removes big-M numerics.
void reduction3(const IloArray<IloNumArray> &costMatrix, int depots, int nodes,
                int matrixSize, vector<char> &assignMatrix,
                IloNumArray maxVehiclesPerDepot)
{
    IloEnv env;
    IloModel relaxMDVSP(env);
    IloCplex solver(relaxMDVSP);
    solver.setOut(env.getNullStream());

    IloRangeArray pi = IloAdd(relaxMDVSP, IloRangeArray(env, nodes, 1, 1));
    IloRangeArray sigma1 = IloAdd(relaxMDVSP, IloRangeArray(env, 0, maxVehiclesPerDepot));
    IloRangeArray sigma2 = IloAdd(relaxMDVSP, IloRangeArray(env, 0, maxVehiclesPerDepot));

    // Track which (i,j) have variables so we can skip them in constraints
    vector<char> hasVar(matrixSize * matrixSize, 0);
    IloArray<IloNumVarArray> x(env, matrixSize);
    for (int k = 0; k < matrixSize; k++)
        x[k] = IloNumVarArray(env, matrixSize);

    // Create variables only for feasible arcs
    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {
            // Skip infeasible arcs and depot-to-depot arcs
            if (costMatrix[i][j] == -1 || (i < depots && j < depots)) {
                x[i][j] = IloNumVar(env, 0, 0);
                continue;
            }

            hasVar[i * matrixSize + j] = 1;
            if (i < depots && j >= depots)
                x[i][j] = IloNumVar(pi[j - depots](1) + sigma1[i](1));
            else if (i >= depots && j < depots)
                x[i][j] = IloNumVar(sigma2[j](1));
            else // i >= depots && j >= depots
                x[i][j] = IloNumVar(pi[j - depots](1));
        }
    }

    // Flow conservation — only feasible arcs
    for (int i = depots; i < matrixSize; i++) {
        IloExpr flowIn(env);
        IloExpr flowOut(env);
        for (int j = 0; j < matrixSize; j++) {
            if (hasVar[j * matrixSize + i]) flowIn += x[j][i];
            if (hasVar[i * matrixSize + j]) flowOut += x[i][j];
        }
        relaxMDVSP.add(flowIn - flowOut == 0);
        flowIn.end();
        flowOut.end();
    }

    // Objective — only feasible arcs, no big-M costs
    IloExpr objective(env);
    for (int i = 0; i < matrixSize; i++)
        for (int j = 0; j < matrixSize; j++)
            if (hasVar[i * matrixSize + j])
                objective += costMatrix[i][j] * x[i][j];

    relaxMDVSP.add(IloMinimize(env, objective));

    if (solver.solve()) {
        for (int i = 0; i < matrixSize; i++)
            for (int j = 0; j < matrixSize; j++)
                if (hasVar[i * matrixSize + j] && solver.getValue(x[i][j]) > 0)
                    assignMatrix[i * matrixSize + j] = 1;
    }

    solver.end();
    relaxMDVSP.end();
    env.end();
}

bool applyReductions(IloArray<IloNumArray> &costMatrix, int depots, int nodes,
                     int matrixSize, IloNumArray maxVehiclesPerDepot,
                     vector<vector<int>> &predSI,
                     vector<vector<vector<bool>>> &assignSImatrix,
                     bool doR1, bool doR2, bool doR3, bool doSI)
{
    if (!doR1 && !doR2 && !doR3 && !doSI)
        return true;

    vector<char> assignMatrix(matrixSize * matrixSize, 0);

#if HAS_JV
    // Initial solution uses the LAP-based reduction2
    if (doSI) {
        reduction2(costMatrix, depots, nodes, matrixSize, assignMatrix, assignSImatrix);
        for (int k = 0; k < depots; k++)
            for (int j = 0; j < matrixSize; j++)
                for (int i = 0; i < matrixSize; i++)
                    if (assignSImatrix[k][j][i])
                        predSI[k][i] = j;
    }
#endif

    // Reductions are cumulative: each adds arcs to assignMatrix
    if (doR1)
        reduction1(costMatrix, depots, nodes, matrixSize, assignMatrix);

#if HAS_JV
    if (doR2 && !doSI)
        reduction2(costMatrix, depots, nodes, matrixSize, assignMatrix, assignSImatrix);
#endif

    if (doR3)
        reduction3(costMatrix, depots, nodes, matrixSize, assignMatrix, maxVehiclesPerDepot);

    // Prune arcs not selected by any reduction
    if (doR1 || doR2 || doR3) {
        for (int i = depots; i < matrixSize; i++)
            for (int j = depots; j < matrixSize; j++)
                if (!assignMatrix[i * matrixSize + j])
                    costMatrix[i][j] = -1;
    }

    return true;
}
