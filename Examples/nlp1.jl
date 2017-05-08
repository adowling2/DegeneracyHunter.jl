include("../DegeneracyHunter.jl")
using JuMP
#using Ipopt
using AmplNLWriter

#s = IpoptSolver(bound_relax_factor=0.0)
s = AmplNLSolver("ipopt")

m = Model(solver=s)

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

println("Manually check constrain violations")

r = zeros(5)
a_ = zeros(5)

for i = 1:5
	a_[i] = getvalue(a[i])
end
	
z_ = getvalue(z)
y_ = getvalue(y)

r[1] = max(10 - a_[1] - a_[2] + a_[3], 0)
r[2] = max(-a_[2] + a_[2]*y_ - a_[1]*a_[4], 0)
r[3] = a_[3]*a_[2] + a_[1]*a_[4] - a_[5]
r[4] = z_ - sum(a_) + 0.001*sum(a_.^2)
r[5] = a_[3] + a_[4] - 10*y_

for i = 1:5
	println("r[",i,"] = ",r[i])
end