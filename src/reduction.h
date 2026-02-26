#ifndef REDUCTION_H
#define REDUCTION_H

#include "types.h"
#include <ilcplex/ilocplex.h>
#include <vector>

// assignMatrix is a flat matrixSize x matrixSize array of char (0/1).
// Index with: assignMatrix[i * matrixSize + j]
// Using char instead of vector<bool> avoids bit-packing overhead and
// enables safe parallel writes to distinct indices.

// R1: keeps sqrt(N) cheapest arcs per row/column using nth_element — O(N^2)
void reduction1(const IloArray<IloNumArray> &costMatrix, int depots, int nodes,
                int matrixSize, std::vector<char> &assignMatrix);

#if HAS_JV
#include "JV/lap.h"
// R2: Jonker-Volgenant LAP-based reduction (per depot, parallelized with OpenMP)
void reduction2(const IloArray<IloNumArray> &costMatrix, int depots, int nodes,
                int matrixSize, std::vector<char> &assignMatrix,
                std::vector<std::vector<std::vector<bool>>> &assignSImatrix);
#endif

// R3: LP relaxation-based reduction (skips infeasible arcs from LP model)
void reduction3(const IloArray<IloNumArray> &costMatrix, int depots, int nodes,
                int matrixSize, std::vector<char> &assignMatrix,
                IloNumArray maxVehiclesPerDepot);

// Apply selected reductions and prune costMatrix.
bool applyReductions(IloArray<IloNumArray> &costMatrix, int depots, int nodes,
                     int matrixSize, IloNumArray maxVehiclesPerDepot,
                     std::vector<std::vector<int>> &predSI,
                     std::vector<std::vector<std::vector<bool>>> &assignSImatrix,
                     bool doR1 = true, bool doR2 = true, bool doR3 = true,
                     bool doSI = true);

#endif
