using JuMP

##### This file looks at solving an optimization in parallel using pmap()
# More info: http://julia.readthedocs.org/en/latest/manual/parallel-computing/
# Use "julia -p n" to start Julia with n 

### Import optimization problem definition/solve in a function

push!(LOAD_PATH,pwd())
using DummyModule

### Test function

my = solveOptProb2((1,1))


### Iterate over possible values of a and b
aOpt = 0:5
bOpt = 1:3

aN = length(aOpt)
bN = length(bOpt)

# Assemble array of input parameters for the optimization problem
P = Array(Tuple,aN*bN)

for i in 1:length(aOpt)
	for j in 1:length(bOpt)
		P[(i-1)*bN + j] = (aOpt[i],bOpt[j])
	end
end

# Solve the optimization problem in parallel
# Note: results is an array of tuples, with the same dimensions as P
results = pmap(solveOptProb2, P)

for i = 1:length(results)
	DummyModule.print(results[i])
end