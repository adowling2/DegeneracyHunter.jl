include("../DegeneracyHunter.jl")
using JuMP
using Ipopt

m = Model(solver=IpoptSolver())

@variable(m, -10 <= a[1:5] <= 10, start=1.0)
@variable(m, -100 <= z <= 100, start=5.0)

# Fix this variable by setting upper bound = lower bound
@variable(m, 0 <= y <= 0)


@constraint(m, a[1] + a[2] - a[3] >= 10)
@NLconstraint(m, a[1]*a[4] - a[2]*y + a[2] >= 0)
@NLconstraint(m, a[3]*a[2] + a[1]*a[4] - a[5] == 0)
@NLconstraint(m, z == sum{a[i], i=1:5} - 0.001*sum{a[i]^2, i=1:5})
@constraint(m, a[3] + a[4] == 10*y)

@objective(m, Min, z)

solve(m)

DegeneracyHunter.printInfeasibleEquations(m, 1E-10)