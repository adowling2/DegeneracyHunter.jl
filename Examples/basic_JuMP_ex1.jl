using JuMP

##### Simple LP
# Example taken directly from https://jump.readthedocs.org/en/latest/quickstart.html

println("***** Problem 1 *****")

m = Model()
@variable(m, 0 <= x <= 2 )
@variable(m, 0 <= y <= 30 )

@objective(m, Max, 5x + 3*y )
@constraint(m, 1x + 5y <= 3.0 )

print(m)

status = solve(m)

println("Objective value: ", getobjectivevalue(m))
println("x = ", getvalue(x))
println("y = ", getvalue(y))

##### Another simple LP with logical conditions in constraint definition
# Example from https://groups.google.com/forum/#!topic/julia-opt/pMHz-9YHN2o
# Note: This example only works with Julia-v0.4 and later

println(" ")
println("***** Problem 2 *****")

m = Model()

I = 1:5

@variable(m, 0 <= x[I] <= 10)
@variable(m, 0 <= y[I] <= 10)

# This constraint is only considered when i + j are less than or equal to 3
# Note: Only supported by Julia-v0.4 and later
@constraint(m, constr[i=I,j=I; i+j <= 3], x[i] - y[j] == 1)

@objective(m, Min, sum{x[i] + y[i], i=I})

print(m)

status = solve(m)

println(" ")
println("Objective value: ", getobjectivevalue(m))
println("x = ", getvalue(x))
println("y = ", getvalue(y))