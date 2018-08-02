using DegeneracyHunter

# Create an optimization model 'm'
include("zero_flowrates.jl")

# Load MILP solver
using Gurobi

# Search for irreducible degenerate sets
sets, dd, ds = DegeneracyHunter.degeneracyHunter(m, GurobiSolver());

println("Done.")
