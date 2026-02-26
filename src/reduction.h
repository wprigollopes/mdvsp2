#ifndef REDUCTION_H
#define REDUCTION_H

#include "types.h"
#include <ilcplex/ilocplex.h>
#include <vector>

// R1: keeps sqrt(N) cheapest arcs per row/column
void reduction1(IloArray<IloNumArray> costMatrix, int depots, int nodes,
                int matrixSize, std::vector<std::vector<bool>> &assignMatrix);

#if HAS_JV
#include "JV/lap.h"
// R2: Jonker-Volgenant LAP-based reduction (per depot)
void reduction2(IloArray<IloNumArray> costMatrix, int depots, int nodes,
                int matrixSize, std::vector<std::vector<bool>> &assignMatrix,
                std::vector<std::vector<std::vector<bool>>> &assignSImatrix);
#endif

// R3: LP relaxation-based reduction
void reduction3(IloArray<IloNumArray> costMatrix, int depots, int nodes,
                int matrixSize, std::vector<std::vector<bool>> &assignMatrix,
                IloNumArray maxVehiclesPerDepot);

// Apply selected reductions and prune costMatrix.
// Simplified from the original 8-case if-else chain: reductions are cumulative.
bool applyReductions(IloArray<IloNumArray> &costMatrix, int depots, int nodes,
                     int matrixSize, IloNumArray maxVehiclesPerDepot,
                     std::vector<std::vector<int>> &predSI,
                     std::vector<std::vector<std::vector<bool>>> &assignSImatrix,
                     bool doR1 = true, bool doR2 = true, bool doR3 = true,
                     bool doSI = true);

#endif
