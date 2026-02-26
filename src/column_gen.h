#ifndef COLUMN_GEN_H
#define COLUMN_GEN_H

#include "types.h"
#include "shortest_path.h"
#include <ilcplex/ilocplex.h>
#include <vector>

class ColumnGenerator {
public:
    ColumnGenerator(IloEnv env, IloArray<IloNumArray> &costMatrix,
                    int depots, int nodes, IloNumArray &maxVehiclesPerDepot);

    void initialize();
    void buildSubproblemGraphs();

#if HAS_JV
    void addInitialColumns(const std::vector<std::vector<int>> &predSI,
                           const std::vector<std::vector<char>> &assignSImatrix);
#endif

    // Returns 0=success, 1=infeasible, 2=iteration limit
    int solve();

private:
    void addColumn(const std::vector<int> &pred, int depot);
    bool checkIntegrality(unsigned int upTo);
    void roundColumns(unsigned int upTo);
    bool feasibility2Test();
    bool feasibility6Test();

    IloEnv env_;
    IloArray<IloNumArray> &costMatrix_;
    int depots_, nodes_, matrixSize_;
    unsigned int capacity_; // sized from instance: nodes * depots * 4
    IloNumArray &maxVehiclesPerDepot_;
    bool useInitialSolution_;

    // Master model
    IloModel masterModel_;
    IloCplex master_;
    IloObjective pathsSelect_;
    IloRangeArray pi_;
    IloRangeArray sigma_;
    IloNumVarArray deltaPi_;
    IloNumVarArray deltaSigma_;

    // Columns
    unsigned int p_;
    IloNumArray c_;
    IloNumVarArray y_;
    IloArray<IloNumArray> a_;
    IloIntArray omega_;
    std::vector<bool> omegasSelected_;

    // Iteration state
    std::vector<int> masterValues_;
    unsigned int iterations_;
    int stabilized_;
    unsigned int lastRoundingIter_;
    bool lbSet_;

    // Duals and feasibility data
    IloNumArray piDuals_;
    IloNumArray sigmaDuals_;
    IloNumArray masterData_;

    // Subproblem graphs (shared trip arcs + per-depot graphs)
    TripArcs tripArcs_;
    std::vector<SubproblemGraph> graphs_;
    std::vector<int> pred_;
    std::vector<long long> dist_;
};

#endif
