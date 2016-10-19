using JuMP
using Ipopt

# Define sets
# Streams
S = 2:9

# Components (chemical species)
C = ["A", "B"]

# Units
U = ["U1", "U2", "U3"]

# Define stream/unit connectivity using dictionaries
inlets = Dict{ASCIIString,Integer}("U1"=>2, "U2"=>5, "U3"=>4)
outletV = Dict{ASCIIString,Integer}("U1"=>3, "U2"=>6, "U3"=>8)
outletL = Dict{ASCIIString,Integer}("U1"=>4, "U2"=>7, "U3"=>9)

# steams
feeds = [2, 5]

# Equilibrium coefficients

K = Dict{Tuple{ASCIIString, ASCIIString}, Float64}(
		("U1","A")=>1.008,		("U1","B")=>0.9,
		("U2","A")=>1.099,		("U2","B")=>0.9,
		("U3","A")=>1.093,		("U3","B")=>0.9	)	 

		
# Feed flowrate specification
# (for components, same for all feed streams)
feedflow = Dict{AbstractString, Float64}("A"=>0.55, "B"=>0.45)

# Define model
m = Model(solver=IpoptSolver())

##### General flowsheet unit model #####

# Total molar flowrate
@variable(m, f[S] >= 0)

# Component molar flowrate
@variable(m, fc[S,C] >= 0)

# Component mole fraction
@variable(m, 0.001 <= x[S,C] <= 1)

# Unit model
for u in U
	
	i = inlets[u]
	v = outletV[u]
	l = outletL[u]
	
	# Overall mass balance
	@constraint(m, f[i] == f[v] + f[l])

	for c in C
		# Component mass balance
		@constraint(m, fc[i,c] == fc[v,c] + fc[l,c])
	
		# vapor-liquid equilibrium
		@constraint(m, x[v, c] == K[(u,c)]*x[l,c])
	
	end
	
	# Unit summation (Rachford-Rice equation)
	# This constraint is redundant!!!
	# @constraint(m, sum{x[v,c] - x[l,c], c in C} == 0)
	
end

# Stream model: mole fraction specification
for s in S
	for c in C
		@constraint(m, fc[s,c] == f[s]*x[s,c])
	end
end

##### Economics and recovery/purity specifications

@variable(m, 0.551 <= purityA <= 1.0)
@variable(m, 0.9 <= recoveryA <= 1.0)

# Equipment cost
ecost = Dict{AbstractString,Float64}("U1"=>1.5, "U2"=>1.0, "U3"=>0.5)

# Set total feed to one
@constraint(m, sum{f[s], s in feeds} == 1)

# Calculate purity of A
@constraint(m, purityA*sum{sum{f[s], s in outletV[u]}, u in U} == sum{sum{fc[s,"A"], s in outletV[u]}, u in U})

# Calculate recovery of A
@constraint(m, recoveryA == sum{sum{fc[s,"A"], s in outletV[u]}, u in U} / feedflow["A"])

# Objective
@objective(m, Min, sum{sum{f[s], s in inlets[u]}*ecost[u], u in U} - 100*purityA)

for s in feeds
	for c in C
		setupperbound(x[s,c], feedflow[c])
		setlowerbound(x[s,c], feedflow[c])
	end
end

##### Initialization

# Guess all molar flowrates are one
# And molar fractions match feed
for s in S
	
	setvalue(f[s], 1.0)
	
	for c in C
		setvalue(fc[s,c], feedflow[c])
		setvalue(x[s,c], feedflow[c])
	end
end

print(m)

solve(m)