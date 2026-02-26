#include "reduction.h"
#include <cmath>

using namespace std;

void reduction1(IloArray<IloNumArray> costMatrix, int depots, int nodes,
                int matrixSize, vector<vector<bool>> &assignMatrix)
{
    unsigned int i, j, x;
    vector<vector<int>> reduceItems(matrixSize, vector<int>(matrixSize, 0));
    unsigned int reduction = ceil(sqrt((float)nodes));
    vector<unsigned int> reduceRow(reduction);
    vector<unsigned int> reduceCol(reduction);
    vector<int> indexRow(reduction);
    vector<int> indexCol(reduction);
    int posChangeRow, posChangeCol;

    for (i = depots; i < (unsigned int)matrixSize; i++) {
        for (x = 0; x < reduction; x++) {
            reduceRow[x] = 1e6;
            reduceCol[x] = 1e6;
            indexRow[x] = -1;
            indexCol[x] = -1;
            posChangeRow = -1;
            posChangeCol = -1;
        }
        for (j = depots; j < (unsigned int)matrixSize; j++) {
            posChangeRow = -1;
            posChangeCol = -1;

            // Row reduction: keep cheapest arcs by row
            if (costMatrix[i][j] != -1) {
                for (x = 0; x < reduction; x++) {
                    if (costMatrix[i][j] <= reduceRow[x]) {
                        posChangeRow = x;
                        break;
                    }
                }
                if (posChangeRow != -1) {
                    for (x = reduction - 1; x > (unsigned int)posChangeRow; x--) {
                        reduceRow[x] = reduceRow[x - 1];
                        indexRow[x] = indexRow[x - 1];
                    }
                    reduceRow[posChangeRow] = costMatrix[i][j];
                    indexRow[posChangeRow] = j;
                }
            }

            // Column reduction: keep cheapest arcs by column
            if (costMatrix[j][i] != -1) {
                for (x = 0; x < reduction; x++) {
                    if (costMatrix[j][i] <= reduceCol[x]) {
                        posChangeCol = x;
                        break;
                    }
                }
                if (posChangeCol != -1) {
                    for (x = reduction - 1; x > (unsigned int)posChangeCol; x--) {
                        reduceCol[x] = reduceCol[x - 1];
                        indexCol[x] = indexCol[x - 1];
                    }
                    reduceCol[posChangeCol] = costMatrix[j][i];
                    indexCol[posChangeCol] = i;
                }
            }
        }
        for (x = 0; x < reduction; x++) {
            if (indexRow[x] != -1)
                reduceItems[i][indexRow[x]] = 1;
            if (indexCol[x] != -1)
                reduceItems[indexCol[x]][i] = 1;
        }
    }
    for (i = depots; i < (unsigned int)matrixSize; i++)
        for (j = depots; j < (unsigned int)matrixSize; j++)
            if (reduceItems[i][j] == 1)
                assignMatrix[i][j] = true;
}

#if HAS_JV

void reduction2(IloArray<IloNumArray> costMatrix, int depots, int nodes,
                int matrixSize, vector<vector<bool>> &assignMatrix,
                vector<vector<vector<bool>>> &assignSImatrix)
{
    int enddim = nodes + nodes;
    cost **origAssignCost = new cost *[enddim];
    for (int i = 0; i < enddim; i++)
        origAssignCost[i] = new cost[enddim];

    col *rowsol = new col[enddim];
    row *colsol = new row[enddim];
    cost *u = new cost[enddim];
    cost *v = new cost[enddim];

    for (int depot = 0; depot < depots; depot++) {
        for (int i = 0; i < nodes; i++) {
            for (int j = 0; j < nodes; j++) {
                origAssignCost[i][j] = (cost)(costMatrix[i + depots][j + depots] == -1
                    ? INFEASIBLE_VALUE : costMatrix[i + depots][j + depots]);
                origAssignCost[i + nodes][j] = (cost)(costMatrix[depot][j + depots] == -1
                    ? INFEASIBLE_VALUE : costMatrix[depot][j + depots]);
                origAssignCost[j][i + nodes] = (cost)(costMatrix[j + depots][depot] == -1
                    ? INFEASIBLE_VALUE : costMatrix[j + depots][depot]);
                origAssignCost[i + nodes][j + nodes] = 0;
            }
        }

        lap(enddim, origAssignCost, rowsol, colsol, u, v);

        for (int dim = 0; dim < enddim; dim++) {
            int from = dim >= nodes ? depot : dim + depots;
            int to = rowsol[dim] >= nodes ? depot : rowsol[dim] + depots;
            assignSImatrix[depot][from][to] = true;
            assignMatrix[from][to] = true;
        }
        assignSImatrix[depot][depot][depot] = false;
        assignMatrix[depot][depot] = false;
    }

    for (int i = 0; i < enddim; i++)
        delete[] origAssignCost[i];
    delete[] origAssignCost;
    delete[] rowsol;
    delete[] colsol;
    delete[] u;
    delete[] v;
}

#endif // HAS_JV

void reduction3(IloArray<IloNumArray> costMatrix, int depots, int nodes,
                int matrixSize, vector<vector<bool>> &assignMatrix,
                IloNumArray maxVehiclesPerDepot)
{
    IloEnv env;
    IloModel relaxMDVSP(env);
    IloCplex solver(relaxMDVSP);
    IloRangeArray pi = IloAdd(relaxMDVSP, IloRangeArray(env, nodes, 1, 1));
    IloRangeArray sigma1 = IloAdd(relaxMDVSP, IloRangeArray(env, 0, maxVehiclesPerDepot));
    IloRangeArray sigma2 = IloAdd(relaxMDVSP, IloRangeArray(env, 0, maxVehiclesPerDepot));

    IloArray<IloNumVarArray> x(env, matrixSize);
    for (unsigned int k = 0; k < (unsigned int)matrixSize; k++)
        x[k] = IloNumVarArray(env, matrixSize);

    for (unsigned int i = 0; i < (unsigned int)matrixSize; i++) {
        for (unsigned int j = 0; j < (unsigned int)matrixSize; j++) {
            if ((int)i < depots && (int)j >= depots)
                x[i][j] = IloNumVar(pi[j - depots](1) + sigma1[i](1));
            else if ((int)i >= depots && (int)j < depots)
                x[i][j] = IloNumVar(sigma2[j](1));
            else if ((int)i >= depots && (int)j >= depots)
                x[i][j] = IloNumVar(pi[j - depots](1));
            else
                x[i][j] = IloNumVar(env);
        }
    }

    for (int i = depots; i < matrixSize; i++) {
        IloExpr flowIn(env);
        IloExpr flowOut(env);
        for (int j = 0; j < matrixSize; j++) {
            flowIn += x[j][i];
            flowOut += x[i][j];
        }
        relaxMDVSP.add(flowIn - flowOut == 0);
        flowIn.end();
        flowOut.end();
    }

    IloExpr objective(env);
    for (unsigned int i = 0; i < (unsigned int)matrixSize; i++)
        for (unsigned int j = 0; j < (unsigned int)matrixSize; j++)
            if (costMatrix[i][j] == -1)
                objective += INFEASIBLE_VALUE * x[i][j];
            else
                objective += costMatrix[i][j] * x[i][j];

    relaxMDVSP.add(IloMinimize(env, objective));

    if (solver.solve()) {
        for (unsigned int i = 0; i < (unsigned int)matrixSize; i++)
            for (unsigned int j = 0; j < (unsigned int)matrixSize; j++)
                if (solver.getValue(x[i][j]) > 0)
                    assignMatrix[i][j] = true;
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

    vector<vector<bool>> assignMatrix(matrixSize, vector<bool>(matrixSize, false));

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
                if (!assignMatrix[i][j])
                    costMatrix[i][j] = -1;
    }

    return true;
}
