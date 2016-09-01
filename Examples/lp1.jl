######################################################
# Created by Alex Dowling (alexdowling.net)
# while at the University of Wisconsin-Madison
# ####################################################

# Define problem
include("../DegeneracyHunter.jl")
using JuMP

# Define problem
function lp1()

	m = Model()

	@variable(m, 0 <= x[1:3] <= 5)
	@variable(m, 0 <= y <= 0)
	@constraint(m, x[1] + x[2] >= 1)
	@constraint(m, x[1] + x[2] + x[3] == 1)
	@constraint(m, x[2] - 2*x[3] <= 1)
	@constraint(m, x[1] + x[3] >= 1)
	@constraint(m, x[1] + x[2] + x[3] == 1) # Redundant constraint - makes problem degenerate
	@objective(m, Min, sum{x[i],i=1:3})
	
	return m
end

function initialize!(m::Model)
	x = getvariable(m,:x)
	setvalue(x[1], 1.0)
	setvalue(x[2], 5.0)
	setvalue(x[3], -1.0)
	
	# Intentionally do not initialize y
	
	return nothing
end

# Create and initialize model
m = lp1()
initialize!(m) 

# Check initialization
DegeneracyHunter.printVariableDiagnostics(m)

# Print equations with absolute residuals greater than 0.001
DegeneracyHunter.printInfeasibleEquations(m, 0.001)

# Check for degeneracy constraints
DegeneracyHunter.degeneracyHunter(m, includeBounds=true)

# Solve model
tic()
status = solve(m)
tm = toq()

# Print inactive constraints at solution
DegeneracyHunter.printInactiveEquations(m)

# Package and print problem size information
ps = assembleProblemStats(m,status,tm)
printProblemStats(ps)