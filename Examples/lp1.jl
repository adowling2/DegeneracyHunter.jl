######################################################
# Created by Alex Dowling (alexdowling.net)
# while at the University of Wisconsin-Madison
# ####################################################

include("../DegeneracyHunter.jl")
using JuMP

m = Model()

@variable(m, 0 <= x[1:3] <= 5)
@variable(m, 0 <= y <= 0)
@constraint(m, x[1] + x[2] >= 1)
@constraint(m, x[1] + x[2] + x[3] == 1)
@constraint(m, x[2] - 2*x[3] <= 1)
@constraint(m, x[1] + x[3] >= 1)
@constraint(m, x[1] + x[2] + x[3] == 1) # Redundant constraint - makes problem degenerate
@objective(m, Min, sum{x[i],i=1:3})

DegeneracyHunter.printVariableDiagnostics(m)

#=

tic()
status = solve(m)
tm = toq()

for i = 1:3
	println("x[$i] = ",string(getValue(x[i])))
end
println("y = ",string(getValue(y)))

println(" ")

=#


(ids, dd, ds) = DegeneracyHunter.degeneracyHunter(m, includeBounds=true);

# ps = assembleProblemStats(m,status,tm)