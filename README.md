# MDVSP Solver

A fast heuristic solver for the **Multiple-Depot Vehicle Scheduling Problem** (MDVSP), implementing the algorithm described in:

> P.C. Guedes, W.P. Lopes, L.R. Rohde, D. Borenstein.
> **Simple and efficient heuristic approach for the multiple-depot vehicle scheduling problem.**
> *Optimization Letters*, 10:1449-1461, 2016.
> DOI: [10.1007/s11590-015-0944-x](https://doi.org/10.1007/s11590-015-0944-x)

## Problem

Given a set of *N* service trips with fixed start/end times and locations, and *K* depots each with a limited fleet of vehicles, the MDVSP finds a minimum-cost assignment of trips to vehicle schedules such that:

- Every trip is covered by exactly one vehicle
- Each vehicle starts and ends at the same depot
- Depot fleet capacities are respected
- Deadheading (empty travel between trips) cost is minimized

The MDVSP is NP-hard (Bertossi et al., 1987).

## Method

The solver is a two-phase heuristic:

### Phase 1: State Space Reduction

Three selection procedures prune the arc set before solving, removing arcs unlikely to appear in good solutions:

| Procedure | Strategy | Complexity |
|-----------|----------|------------|
| **R1** (k-SDVSP) | Solve each single-depot VSP via Jonker-Volgenant LAPJV; keep arcs appearing in any solution | O(\|K\| n^2 log n) |
| **R2** (Relaxed-MDVSP) | Solve the LP relaxation of the multicommodity flow formulation; keep arcs with positive flow | LP solve |
| **R1+R2** | Union of both selections | Combined |

R1 uses the LAPJV algorithm (Jonker & Volgenant, 1987) to solve each depot's single-depot VSP as an assignment problem. R2 solves a relaxed MDVSP where vehicles may end at a different depot than where they started.

### Phase 2: Modified Truncated Column Generation (MTCG)

The reduced problem is solved via Dantzig-Wolfe decomposition:

- **Master problem (RMP):** Set partitioning/covering formulation solved as an LP by CPLEX. Uses a relaxed formulation (covering constraints >= 1) initially, switching to equality constraints when stabilized.
- **Subproblems:** One shortest-path problem per depot on the MDVSP network. Trips are sorted by start time, so the subproblem graph is a DAG -- solved in O(V+E) by topological relaxation.
- **Column insertion:** Unlike the standard approach of adding up to |K| columns per iteration, each improving column is inserted immediately and the master is re-solved. This updates duals more frequently, stabilizing convergence.
- **Integrality:** The CG terminates when the master stagnates (Z_min, I parameters). Fractional columns are rounded using an Omega_min threshold. The LAPJV initial solution provides a warm start.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Z_min | 0 | Minimum improvement to continue CG |
| I (ITER_MIN) | 30 | Minimum iterations between rounding phases |
| Omega_min | 0.9 | Threshold for rounding fractional columns to 1 |

## Results

Tested on Huisman's benchmark instances (4/8 depots, 500/1000/1500 trips):

- **Solution quality:** Average gap below 0.7% from best-known solutions across all variants
- **R1+R2 variant:** Best trade-off -- 0.51% average gap, 16x faster than previous CG results
- **R1 variant:** Fastest -- 0.62% average gap, 34x faster than previous results
- All variants find the optimal number of vehicles for every instance

## Project Structure

```
src/
  main.cc              Entry point, data loading
  types.h              Constants, Arc struct
  shortest_path.h/cc   DAG shortest path (O(V+E)), SPFA fallback
  reduction.h/cc       R1 (LAPJV), R2 (LP relaxation), R3 (sqrt(N) pruning)
  column_gen.h/cc      MTCG solver (master LP, subproblems, rounding)
JV/
  lap.h, lap.cpp       Jonker-Volgenant LAPJV algorithm
  gnrl.h               JV type definitions
Makefile               Build with CPLEX and OpenMP
```

## Dependencies

- **IBM CPLEX** (Concert Technology C++ API) -- LP/MIP solver for the master problem and R2 reduction
- **GCC** with C++17 and OpenMP support

## Building

```bash
# Set CPLEX version and system architecture
make CPLEXVERSION=2211 SYSTEM=x86-64_linux

# Or edit the defaults in Makefile and run
make
```

CPLEX is expected at `/opt/ibm/ILOG/CPLEX_Studio<VERSION>/`.

## Usage

```bash
./mdvsp <instance_file>
```

Instance files follow CPLEX Concert text format:
```
<depots>
<nodes>
[max_vehicles_per_depot]
[[cost_matrix]]
```

Where the cost matrix is (depots+nodes) x (depots+nodes), with -1 indicating infeasible arcs.

Benchmark instances are available at [Huisman's website](http://people.few.eur.nl/huisman/instances.htm).

## License

Academic use. The LAPJV implementation is by Roy Jonker (MagicLogic Optimization Inc.), modified by Yong Yang.
