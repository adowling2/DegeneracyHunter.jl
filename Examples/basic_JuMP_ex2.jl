using JuMP

### Place optimization problem definition and solve in a function
function solveOptProb1(a,b)
	m = Model()
	@variable(m, 0 <= x <= 2 )
	@variable(m, 0 <= y <= 30 )

	@objective(m, Max, a*x + 3*y )
	@constraint(m, 1x + 5y <= b )
	
	status = solve(m)
	
	return getobjectivevalue(m), [getvalue(x); getvalue(y)], status
	
end

### Test function

(obj, z, status) = solveOptProb1(1,1)


### Iterate over possible values of a and b
aOpt = 0:5
bOpt = 1:3

aN = length(aOpt)
bN = length(bOpt)

objValues = zeros(aN,bN)
zValues = zeros(2,aN,bN)

for i in 1:length(aOpt)
	for j in 1:length(bOpt)
		(objValues[i,j], zValues[:,i,j], status) = solveOptProb1(aOpt[i],bOpt[j])
	end
end

println("objValues = ")
println(objValues)
println(" ")

println("zValues = ")
println(zValues)
println(" ")