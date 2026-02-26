#include "column_gen.h"
#include <iostream>

using namespace std;

ColumnGenerator::ColumnGenerator(IloEnv env, IloArray<IloNumArray> &costMatrix,
                                 int depots, int nodes,
                                 IloNumArray &maxVehiclesPerDepot)
    : env_(env), costMatrix_(costMatrix), depots_(depots), nodes_(nodes),
      matrixSize_(depots + nodes), maxVehiclesPerDepot_(maxVehiclesPerDepot),
      useInitialSolution_(false), masterModel_(env), master_(masterModel_),
      p_(0), c_(env, ITERLIM), y_(env, ITERLIM), a_(env, ITERLIM),
      omega_(env, ITERLIM), omegasSelected_(ITERLIM, false),
      masterValues_(ITERLIM, 0), iterations_(1), stabilized_(0),
      lastRoundingIter_(0), lbSet_(CG_TYPE != 1), piDuals_(env),
      sigmaDuals_(env), masterData_(env), graphs_(depots),
      pred_(matrixSize)
{
}

void ColumnGenerator::initialize()
{
    pathsSelect_ = IloAdd(masterModel_, IloMinimize(env_));
    pi_ = IloAdd(masterModel_,
                 IloRangeArray(env_, nodes_, 1, (CG_TYPE == 1) ? IloInfinity : 1));
    sigma_ = IloAdd(masterModel_, IloRangeArray(env_, 0, maxVehiclesPerDepot_));

    deltaPi_ = IloNumVarArray(env_, nodes_, 0, IloInfinity, ILOFLOAT);
    deltaSigma_ = IloNumVarArray(env_, depots_, 0, IloInfinity, ILOFLOAT);

    for (int j = 0; j < nodes_; j++)
        deltaPi_[j] = IloNumVar(pathsSelect_((IloNum)BIGM) + pi_[j](1));

    for (int k = 0; k < depots_; k++)
        deltaSigma_[k] = IloNumVar(sigma_[k](1));

    master_.setOut(env_.getNullStream());
    masterValues_[0] = nodes_ * BIGM;
}

void ColumnGenerator::buildSubproblemGraphs()
{
    for (int k = 0; k < depots_; k++)
        graphs_[k].build(costMatrix_, depots_, k);
}

#if HAS_JV
void ColumnGenerator::addInitialColumns(
    const vector<vector<int>> &predSI,
    const vector<vector<vector<bool>>> &assignSImatrix)
{
    useInitialSolution_ = true;

    for (int k = 0; k < depots_; k++) {
        for (int i = 0; i < matrixSize_; i++) {
            if (!assignSImatrix[k][i][k])
                continue;

            a_[p_] = IloNumArray(env_, nodes_);
            int pathCost = 0;

            // Reconstruct path from JV predecessor data
            pathCost += costMatrix_[i][k];
            int pos = i;
            bool valid = true;

            while (pos != -1) {
                if (predSI[k][pos] < depots_ && predSI[k][pos] != -1) {
                    if (costMatrix_[k][pos] == -1) {
                        valid = false;
                        break;
                    }
                    pathCost += costMatrix_[k][pos];
                    a_[p_][pos - depots_] = 1;
                    break;
                } else if (predSI[k][pos] == -1) {
                    cerr << "Error in initial solution path reconstruction" << endl;
                    valid = false;
                    break;
                } else {
                    if (costMatrix_[predSI[k][pos]][pos] == -1) {
                        valid = false;
                        break;
                    }
                    pathCost += costMatrix_[predSI[k][pos]][pos];
                    a_[p_][pos - depots_] = 1;
                    pos = predSI[k][pos];
                }
            }

            if (valid) {
                y_[p_] = IloNumVar(pathsSelect_(pathCost) + pi_(a_[p_]) + sigma_[k](1));
                y_[p_].setBounds(0, 1);
                omegasSelected_[p_] = false;
                p_++;
            }
        }
    }
}
#endif

void ColumnGenerator::addColumn(const vector<int> &pred, int depot)
{
    int N = nodes_;
    a_[p_] = IloNumArray(env_, N);
    c_[p_] = 0;
    omega_[p_] = depot;

    // Walk predecessor chain from sink back to source
    int pos = pred[N + 1]; // predecessor of sink
    c_[p_] += costMatrix_[pos + depots_][depot]; // return-to-depot arc

    while (pos != N) { // N = source node
        if (pred[pos] == N) {
            // First trip: depot -> trip
            c_[p_] += costMatrix_[depot][pos + depots_];
        } else {
            // Trip-to-trip arc
            c_[p_] += costMatrix_[pred[pos] + depots_][pos + depots_];
        }
        a_[p_][pos] = 1;
        pos = pred[pos];
    }

    y_[p_] = IloNumVar(pathsSelect_(c_[p_]) + pi_(a_[p_]) + sigma_[depot](1));
    y_[p_].setBounds(0, 1);
    omegasSelected_[p_] = false;
    p_++;
}

bool ColumnGenerator::feasibility2Test()
{
    for (int j = 0; j < nodes_; j++)
        if (masterData_[j] != 0)
            return false;
    return true;
}

bool ColumnGenerator::feasibility6Test()
{
    for (int j = 0; j < nodes_; j++) {
        if (masterData_[j] != 0) {
            cout << "Infeasible solution, trying again" << endl;
            return false;
        }
    }
    return true;
}

bool ColumnGenerator::checkIntegrality(unsigned int upTo)
{
    for (unsigned int i = 0; i < upTo; i++)
        if (master_.getValue(y_[i]) != 0 && master_.getValue(y_[i]) != 1)
            return false;
    return true;
}

void ColumnGenerator::roundColumns(unsigned int upTo)
{
    bool checkOmega = false;
    float lastOmega = 0;
    unsigned int lastp = 0;

    for (unsigned int i = 0; i < upTo; i++) {
        if (omegasSelected_[i])
            continue;
        double val = master_.getValue(y_[i]);
        if (val >= OMEGA_MIN && val != 1) {
            y_[i].setLB(1);
            omegasSelected_[i] = true;
            master_.solveFixed();
            checkOmega = true;
        }
        if (lastOmega < val && val != 1) {
            lastOmega = val;
            lastp = i;
        }
    }
    if (!checkOmega) {
        omegasSelected_[lastp] = true;
        masterModel_.add(y_[lastp] == 1);
    }
}

int ColumnGenerator::solve()
{
    bool continueCG = true;

    while (continueCG) {
        unsigned int pPrev = p_;
        bool doFeasibility6 = false;
        continueCG = false;

        master_.solve();

        // Check stabilization
        if (masterValues_[iterations_ - 1] == (int)master_.getObjValue())
            stabilized_++;
        else
            stabilized_ = 0;

        masterValues_[iterations_] = master_.getObjValue();
        iterations_++;

        if (stabilized_ > NUM_STABILIZED) {
            continueCG = false;
        } else {
            // Get delta variable values (only needed without initial solution)
            if (!useInitialSolution_)
                master_.getValues(deltaPi_, masterData_);

            // Check if we should transition to feasibility/rounding phase
            if (useInitialSolution_) {
                if ((iterations_ - lastRoundingIter_ > ITER_MIN) &&
                    (masterValues_[iterations_ - 2] - masterValues_[iterations_ - 1] <= (int)ZMIN))
                    doFeasibility6 = true;
            } else if (feasibility2Test()) {
                if ((iterations_ - lastRoundingIter_ > ITER_MIN) &&
                    (masterValues_[iterations_ - 2] - masterValues_[iterations_ - 1] <= (int)ZMIN))
                    doFeasibility6 = true;
            }

            // Column generation: solve subproblems and add improving columns
            if (!doFeasibility6) {
                master_.getDuals(piDuals_, pi_);
                master_.getDuals(sigmaDuals_, sigma_);

                for (int k = 0; k < depots_; k++) {
                    long long spDist = dag_shortest_path(
                        graphs_[k], pred_, piDuals_, sigmaDuals_, k);

                    if (spDist < 0) {
                        addColumn(pred_, k);
                        continueCG = true;
                        master_.solve();
                        master_.getDuals(piDuals_, pi_);
                        master_.getDuals(sigmaDuals_, sigma_);
                    }
                }
            }
        }

        if (!continueCG)
            doFeasibility6 = true;

        if (doFeasibility6) {
            // Check feasibility of delta variables
            if (!useInitialSolution_) {
                if (!feasibility6Test()) {
                    cout << "No solution found, sorry :/" << endl;
                    return 1;
                }
            }

            if (!lbSet_) {
                // First transition: tighten covering constraints to equality
                for (int j = 0; j < nodes_; j++)
                    pi_[j].setBounds(1, 1);
                lastRoundingIter_ = iterations_;
                lbSet_ = true;
            } else {
                if (iterations_ > ITER_MIN) {
                    // Check if we have an integer solution
                    if (checkIntegrality(pPrev)) {
                        int nVehicles = 0;
                        for (unsigned int i = 0; i < p_; i++)
                            if (master_.getValue(y_[i]) != 0)
                                nVehicles++;
                        cout << "Number of vehicles: " << nVehicles << endl;
                        cout << "Objective value: " << master_.getObjValue() << endl;
                        return 0;
                    }
                }
                // Omega-rounding: fix fractional columns
                roundColumns(pPrev);
                lastRoundingIter_ = iterations_;
            }
        }
        continueCG = true;
    }

    if (p_ >= ITERLIM) {
        cout << "Iteration limit exceeded, result: "
             << masterValues_[iterations_ - 1] << endl;
    }
    return 2;
}
