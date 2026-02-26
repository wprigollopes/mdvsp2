#include "types.h"
#include "reduction.h"
#include "column_gen.h"
#include <ilcplex/ilocplex.h>
#include <cstdlib>
#include <string>
#include <iostream>

ILOSTLBEGIN

using namespace std;

static void readdata(const char *filename,
                     IloArray<IloNumArray> &costMatrix,
                     unsigned short int &depots,
                     unsigned int &nodes,
                     IloNumArray &maxVehiclesPerDepot)
{
    ifstream in(filename);
    if (!in) {
        cerr << "File not found: " << filename << endl;
        exit(1);
    }
    in >> depots;
    in >> nodes;
    in >> maxVehiclesPerDepot;
    in >> costMatrix;
}

int main(int argc, char **argv)
{
    IloEnv env;
    try {
        unsigned short int depots;
        unsigned int nodes;
        IloNumArray maxVehiclesPerDepot(env);
        IloArray<IloNumArray> costMatrix(env);

        const char *instanceFile = (argc > 1)
            ? argv[1]
            : "instancias/pepin/m4n500s0.inp";

        readdata(instanceFile, costMatrix, depots, nodes, maxVehiclesPerDepot);
        cout << instanceFile << endl;

        unsigned int matrixSize = depots + nodes;

        cout << "\n\n ========================================= " << endl;
        cout << instanceFile << endl
             << "Depots: " << depots << endl
             << "Nodes: " << nodes << endl
             << "OmegaMin: " << OMEGA_MIN << endl
             << "IterMin: " << ITER_MIN << endl
             << "Z min: " << ZMIN << endl
             << "Generation mode: "
             << ((CG_TYPE == 1) ? "With relaxation" : "Without relaxation")
             << endl;

        // Reduction data structures
        vector<vector<char>> assignSImatrix(depots,
            vector<char>(matrixSize * matrixSize, 0));

        vector<vector<int>> predSI(depots);
        for (int k = 0; k < depots; k++)
            predSI[k].assign(matrixSize, 0);

#if HAS_JV
        bool useInitialSolution = true;
        cout << "With Initial Solution: Yes" << endl;
#else
        bool useInitialSolution = false;
        cout << "With Initial Solution: No (JV not available)" << endl;
#endif
        cout << "==============================================" << endl;

        // Run reductions and optionally compute initial solution
        applyReductions(costMatrix, depots, nodes, matrixSize,
                        maxVehiclesPerDepot, predSI, assignSImatrix,
                        true, true, true, useInitialSolution);

        // Set up and run column generation
        ColumnGenerator cg(env, costMatrix, depots, nodes, maxVehiclesPerDepot);
        cg.initialize();
        cg.buildSubproblemGraphs();

#if HAS_JV
        if (useInitialSolution)
            cg.addInitialColumns(predSI, assignSImatrix);
#endif

        int result = cg.solve();
        env.end();
        return result;
    }
    catch (IloException &ex) {
        cerr << "Error: " << ex << endl;
    }
    catch (...) {
        cerr << "Error" << endl;
    }
    env.end();
    return 1;
}
