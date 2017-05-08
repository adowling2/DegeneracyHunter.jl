module DegeneracyHunter

export IrreducibleDegenerateSet, DegenSettings, DegenData, printBound, printVariablesInEquation,
	printIDS, degeneracyHunter, printInfeasibleEquations, printInactiveEquations,
	checkVarBounds, printRows, ProblemStats, assembleProblemStats, stringProblemStats, printProblemStats,
	printVariableDiagnostics, round_to_bounds!, checkEquationScaling

importall JuMP
import MathProgBase
# using Gurobi

##############
##### Functions and Types for setup

type IrreducibleDegenerateSet
	c::Int64
	nElements::Float64
	lambda::Array{Float64,1}
	y::Array{Float64,1}
end

type DegenSettings
	includeBounds::Bool
	includeWeaklyActive::Bool
	removeFixedVar::Bool
	onlyCandidateSearch::Bool
	epsiActive::Float64
	epsiLambda::Float64
	lambdaM::Float64
end

function DegenSettings()
	return DegenSettings(false, true, true, true, 1E-6, 1E-6, 1.0E5)
end

type DegenData

# Point to analyze
	x::Array{Float64,1}

# Jacobian expressed as a sparse matrix
	J_active

#	NLP evaluator
	d

# Full Jacobian (pull entries for bounds, if appropriate per settings)
	J::Array{Float64,1}
	iR::Array{Int64,1}
	jC::Array{Int64,1}

# Constraint information
	gLB::Array{Float64,1}
	gUB::Array{Float64,1}
	g::Array{Float64,1}

# Maps from active Jacobian to full Jacobian
	gMap::Array{Int64,1}
	bMap::Array{Int64,1}
	vMap::Array{Int64,1}

# Size info
	nLambda::Int64
	nVarActive::Int64
	nVar::Int64
	nConstr::Int64

end

function DegenData()
	return DegenData(Float64[],
						Void,
						Void,
						Float64[],
						Int64[],
						Int64[],
						Float64[],
						Float64[],
						Float64[],
						Int64[],
						Int64[],
						Int64[],
						0,
						0,
						0,
						0)
end

function DegenData(m::Model, f=STDOUT, status::Symbol=:Unknown)

	dd = DegenData()

# Select point to analyze
	dd.x = get_x(m, status)

	dd.nVar = length(dd.x)

# Evaluate constraint bounds
	(dd.gLB, dd.gUB) = JuMP.constraintbounds(m)

	if(!m.internalModelLoaded)
		JuMP.build(m)
	end
	
	if(typeof(m.internalModel) == MathProgBase.SolverInterface.NonlinearToLPQPBridge)
		processNonlinearModel!(m, dd, f)
	
	# TO DO: Update processLinearModel! to work with quadratic constraints
	elseif(typeof(m.internalModel) <: MathProgBase.SolverInterface.AbstractLinearQuadraticModel)
		# processLinearModel!(m, dd, f)
		processNonlinearModel!(m, dd, f)
		
	elseif(typeof(m.internalModel) <: MathProgBase.SolverInterface.AbstractNonlinearModel)
		processNonlinearModel!(m, dd, f)
	
	else
		println("Warning: Conic models are not yet supported!")
	end

	dd.nConstr = MathProgBase.numconstr(m)

	return dd

end

function get_x(m::Model, status::Symbol = :Unknown)
	if(typeof(m.internalModel) == Void || status == :InfeasibleOrUnbounded || status == :Unsolved)
		x = m.colVal
	else
	# Note: getsolution will throw an error if model is infeasible or unbounded with Gurobi
		x = MathProgBase.getsolution(m.internalModel)
	end
	
	return x
end

function processNonlinearModel!(m::Model, dd::DegenData, f=STDOUT)

# Initialize NLP evaluatorSetup

	print(f,"Setting up NLP evaluator... ")
	tic()
	dd.d = JuMP.NLPEvaluator(m)
	tm = toq()
	println(f, string(tm, " seconds"))

	print(f,"Initializing... ")
	tic()
	MathProgBase.initialize(dd.d, [:ExprGraph, :Jac])
	tm = toq()
	println(f, string(tm," seconds"))

# Evaluate constraints
	n = MathProgBase.numconstr(m)

	dd.g = zeros(n)
	MathProgBase.eval_g(dd.d,dd.g,dd.x)

# Query Jacobian structure

	# Rows: Constraints
	# Columns: Variables

	print(f,"Determining Jacobian structure... ")
	tic()
	(dd.iR, dd.jC) = MathProgBase.jac_structure(dd.d)
	tm = toq()
	println(f,string(tm," seconds"))

	dd.J = zeros(length(dd.iR))

# Query Jacobian at point
	print(f,"Evaluating Jacobian... ")
	tic()
	MathProgBase.eval_jac_g(dd.d,dd.J,dd.x)
	tm = toq()
	println(f,string(tm, " seconds"))

	return nothing

end

# TODO: Update to work with quadratic models.
function processLinearModel!(m::Model, dd::DegenData, f=STDOUT)

	# Jacobian as a SparseMatrixCSC
	if(!m.internalModelLoaded)
		JuMP.build(m)
	end
	A = MathProgBase.getconstrmatrix(internalmodel(m))

	# Translate to row, column, element form
	(dd.iR, dd.jC, dd.J) = findnz(A)

	# Get RHS of constraints
	# dd.g = MathProgBase.getconstrsolution(internalmodel(m))
	dd.g = A*dd.x

	return nothing

end


###############
##### Exported functions for model diagnostics

function round_to_bounds!(m::Model, epsilon::Float64=1.0E-6, f=STDOUT)
	# This function rounds each element of the solution vector to the closest bound if
	# it is within epsilon.
	# Why is this necessary? Some solvers (e.g., IPOPT) give solutions that are slightly
	# infeasible or violate bounds. This is problematic in a meta-algorithm with calls
	# to other solvers.
	
	x = get_x(m, :Unsolved)
	
	up = m.colUpper
	lo = m.colLower
	
	for i = 1:length(x)
	
		v = Variable(m, i)
	
		if(up[i] < lo[i])
			println(f,"Warning: Upper bound (",up[i],") is smaller than lower bound (",lo[i],") for variable ",v)
		end
	
		x_ = NaN
	
		if(x[i] > up[i])
			x_ = up[i]
		elseif(x[i] < lo[i])
			x_ = lo[i]
		elseif(up[i] - x[i] < epsilon)
			x_ = up[i]
		elseif(x[i] - lo[i] < epsilon)
			x_ = lo[i]
		end
		
		if(!isnan(x_))
			setvalue(v, x_)
		end
	
	end
	
	return nothing

end

function printVariableDiagnostics(m::Model, epsilon::Float64=1.0E-6, f=STDOUT, status::Symbol=:Unknown)
	x = get_x(m, status)

	println(f,"Uninitialized Variables: ")
	for i = 1:length(x)
		if(isnan(x[i]))
			v = Variable(m, i)
			#println("x[",i,"] = ",v)
			println(f,v)
		end
	end
	println(f," ")

	println(f,"Variable Violating Bounds: ")
	up = m.colUpper
	lo = m.colLower

	for i = 1:length(x)
		viol_up = x[i] - up[i] > epsilon
		viol_lo = x[i] - lo[i] < -epsilon
		if(viol_up || viol_lo)

			v = Variable(m, i)
			print(f,v, " = ",x[i], ", bounds: [ ")

			if(viol_lo)
				print_with_color(:red, f, string(lo[i]))
			else
				print(f, string(lo[i]))
			end

			print(f, " , ")

			if(viol_up)
				print_with_color(:green, f, string(up[i]))
			else
				print(f, string(up[i]))
			end

			println(f, " ]")
		end
	end

	return nothing

end

function printInfeasibleEquations(m2::Model, eqns, f=STDOUT, status::Symbol=:Unknown)

	dd = DegenData(m2, status)
	return printInfeasibleEquations(m2, dd, eqns, f)

end

function printInfeasibleEquations(m2::Model, dd::DegenData, eqns, f=STDOUT)

	println(f,"Checking infeasibility for specified equations: ")

	r = reshape(min(dd.g - dd.gLB, 0) + max(dd.g - dd.gUB, 0),length(dd.g))

	for i = eqns
		print(f, "r[",string(i),"] = ")
		print_with_color(:green,f,string(r[i]))
		print(f,"\n")

		printEquation(m2, dd, i, f)

		println(f, " ")

	end

	return nothing

end

function printInfeasibleEquations(m2::Model, epsilon::Float64, f=STDOUT, status::Symbol=:Unknown)

	dd = DegenData(m2, f, status)
	return printInfeasibleEquations(m2, dd, epsilon, f)
end

function printInfeasibleEquations(m2::Model, dd::DegenData, epsilon::Float64, f=STDOUT)

	println(f,"Infeasible equations: ")

	#=
	println("size(dd.g) = ",size(dd.g))
	println("size(dd.gLB) = ",size(dd.gLB))
	println("size(dd.gUB) = ",size(dd.gUB))
	=#

	r = reshape(min(dd.g - dd.gLB, 0) + max(dd.g - dd.gUB, 0),length(dd.g))

	k = sortperm(abs(r))

	for j = 1:length(k)

		i = k[j]

		if(abs(r[i]) >= epsilon || isnan(r[i]))
			print(f,"r[",string(i),"] = ")
			print_with_color(:green,f,string(r[i]))
			print(f,"\n")

			printEquation(m2, dd, i,f)

			printVariablesInEquation(m2, dd, i, true, true, f)

			println(f," ")

		end

	end

	return nothing

end

function printInactiveEquations(m2::Model, epsilon::Float64=1E-6, f=STDOUT)

	dd = DegenData(m2, f)
	return printInactiveEquations(m2, dd, epsilon, f)

end

function printInactiveEquations(m2::Model, dd::DegenData, epsilon::Float64=1E-6, f=STDOUT)

	println(f,"Inactive Equations: ")

	inactive = ((dd.gUB - dd.g) .> epsilon) & ((dd.g - dd.gLB) .> epsilon)

	for i = 1:length(dd.g)
		if(((dd.gUB[i] - dd.g[i]) .> epsilon) & ((dd.g[i] - dd.gLB[i]) .> epsilon))
			printEquation(m2, dd, i, f)
		end
	end

	return nothing

end

function checkVarBounds(m::Model, f=STDOUT)

	low = m.colLower
	up = m.colUpper
	val = m.colVal

	n = length(val)

	for i = 1:n
		if(up[i] < low[i] || val[i] < low[i] || val[i] > up[i])

			print(f,getname(m,i),":\t\t","Lower: ")

			if(up[i] < low[i] || val[i] < low[i])
				print_with_color(:red,f,string(low[i]))
			else
				print(f,low[i])
			end

			print(f,"\t Upper: ")

			if(up[i] < low[i] || val[i] > up[i])
				print_with_color(:blue,f,string(up[i]))
			else
				print(f,up[i])
			end

			print(f,"\t Value: ")

			if(val[i] < low[i] || val[i] > up[i])
				print_with_color(:green,f,string(val[i]))
			else
				print(f,val[i])
			end

			print(f,"\n")
		end
	end

	return nothing

end

function printBound(m::Model, i::Int64, epsiActive::Float64, f=STDOUT)

	v = Variable(m, i)

	up = m.colUpper[i]
	lo = m.colLower[i]
	val = m.colVal[i]

	if(abs(val-lo) < epsiActive)
		print(f,string(lo)," <= ")
	end

	print(getname(v))

	if(abs(val-up) < epsiActive)
		print(f," <= ",string(up))
	end

	print(f,"\n")

	return nothing

end

function printRows(m::Model, rows::Int64, f=STDOUT; lite::Bool=false)
	return printRows(m, Array([rows]), f, lite)
end

function printRows(m::Model, rows::Array{Int64,1}, f=STDOUT; lite::Bool=false)

	dd = DegenData(m, f)
	return printRows(m, dd, rows, f, lite)

end

function printRows(m::Model, dd::DegenData, rows::Int64, f=STDOUT; lite::Bool=false)
	return printRows(m, dd, Array([rows]), f, lite=lite)
end

function printRows(m::Model, dd::DegenData, rows::Array{Int64,1}, f=STDOUT; lite::Bool=false)

	if(lite && length(rows) == 1)
	
		r = rows[1]
		
		if(r <= dd.nConstr)
			printEquation(m, dd, rows[1], f)
		else
			printBound(m, r - dd.nConstr, 1E-6, f)
		end
	else

		for i = 1:length(rows)
			r = rows[i]
			
			println(f,"**********************************")
			
			if(r <= dd.nConstr)
					
				println(f,"Constraint ",string(r),"...")
				printEquation(m, dd, r, f)
				println(f,string(dd.gLB[r])," <= ", string(dd.g[r]), " <= ",string(dd.gUB[r]))

				println(f," ")
				println(f,"Involved Variables:")
				k = find(dd.iR .== r)

				for j in unique(dd.jC[k])
					v = Variable(m,j)
					print(f,getname(v)," \t")
					println(f,getlowerbound(v)," <= ",getvalue(v)," <= ",getupperbound(v))
				end
			
			else
				v = r - dd.nConstr
				println(f,"Bound for variable ",v)
				printBounds(m, v, 1E-6, f)
				
			end

			println(f," ")

		end
	end


	return nothing

end

function checkEquationScaling(m::Model, status::Symbol=:Unknown, f=STDOUT)

	dd = DegenData(m, status)
	
	checkJacobianScaling(m, dd, f)

	return nothing

end

##############
##### Degeneracy Hunter features

function degeneracyHunter(m::Model;
	includeBounds::Bool = false,
	includeWeaklyActive::Bool = true,
	removeFixedVar::Bool = true,
	onlyCandidateSearch::Bool = true,
	epsiActive::Float64 = 1E-6,
	epsiLambda::Float64 = 1E-6,
	lambdaM::Float64 = 1E5,
	f=STDOUT)

	ds = DegenSettings(includeBounds, includeWeaklyActive, removeFixedVar,
		onlyCandidateSearch, epsiActive, epsiLambda, lambdaM)

	return degeneracyHunter(m, ds, f)

end

function degeneracyHunter(m::Model, ds::DegenSettings, f)

	println(f,"******************************************")
	println(f,"Welcome to Degeneracy Hunter!")

	println(f," ")
	println(f,"Options: ")
	println(f,"includeBounds = ",string(ds.includeBounds))
	println(f,"includeWeaklyActive = ",string(ds.includeWeaklyActive))
	println(f,"removeFixedVar = ",string(ds.removeFixedVar))
	println(f,"epsiActive = ",string(ds.epsiActive))
	println(f,"epsiLambda = ",string(ds.epsiLambda))
	println(f,"lambdaM = ",string(ds.lambdaM))
	println(f," ")

### Step 0: Setup

# Setup main data structure, evaluate full Jacobian

	dd = DegenData(m, f)

	checkJacobianScaling(m, dd, f)


### Step 1: Assemble "active" Jacobian


# Determine which constraints and bounds are active

	# Grab all equality constraints and inequality constraints that are active
	gLow = (abs(dd.gLB - dd.g) .< ds.epsiActive)
	gUp = (abs(dd.gUB - dd.g) .< ds.epsiActive)
	gEql = (abs(dd.gUB - dd.gLB) .< ds.epsiActive)

	gActive = gLow | gUp | gEql

	# Tip: "$" is the bitwise xor operator
	actIneql = sum(gLow $ gUp)
	eql = sum(gEql)
	totalIneql = length(dd.g) - eql
	inactiveIneql = totalIneql - actIneql

	println(f," ")
	println(f,"(Weakly) Active Inequality Constraints: ",string(sum(gLow $ gUp)))
	println(f,"Inactive Inequality Constraints: ",string(inactiveIneql))
	println(f,"Equality Constraints: ",string(eql))
	println(f," ")


	if(ds.includeWeaklyActive)
		# Do nothing
	else
		# Query multipliers

		# Need to implement

	end

	vFixed = abs(m.colLower - m.colUpper) .< ds.epsiActive

	bActive = (abs(m.colLower - dd.x) .< ds.epsiActive) | (abs(m.colUpper - dd.x) .< ds.epsiActive)

	println(f,"Variables: ",string(length(dd.x)))
	println(f,"(Weakly) Active Variable Bounds: ",string(sum(bActive)))
	println(f,"Fixed Variables: ",string(sum(vFixed)))
	println(f," ")

	if(ds.includeBounds)

		if(ds.includeWeaklyActive)
			# Do nothing
		else
			# Query multiplier

			# Need to implement
		end
	else
		bActive[:] = false
	end


	# Remove fixed variables from bActive
	if(ds.removeFixedVar)
		bActive = bActive & !vFixed
	end

# Add only active bounds to Jacobian
	print(f,"Adding Jacobian elements for bounds... ")
	tic()

	# Number of constraints
	n = length(dd.g)

	dd.bMap = find(bActive)
	for i = 1:length(dd.bMap)
		j = dd.bMap[i]
		push!(dd.iR, i + n)
		push!(dd.jC, j)
		push!(dd.J, 1.0)
	end

	tm = toq()
	println(f,string(tm," seconds"))

	print(f,"Assembling J_sparse... ")
	tic()
	J_sparse = sparse(dd.iR, dd.jC, dd.J, n + length(dd.bMap), dd.nVar)
	tm = toq()
	println(f,string(tm," seconds"))

# Remove inactive constraints and variables from the Jacobian

	dd.gMap = find(gActive)

	if(ds.removeFixedVar)
		dd.vMap = find(!vFixed)
		cols = !vFixed
	else
		dd.vMap = 1:nVar
		cols = trues(nVar)
	end


	rows = [gActive; ones(length(dd.bMap)) .> 0]

	print(f,"Assembling J_active... ")
	tic()
	dd.J_active = J_sparse[rows, cols]
	tm = toq()
	println(f, string(tm," seconds"))

	dd.nLambda = sum(rows)
	dd.nVarActive = sum(cols)

# Step 3: Check if any equations/bounds are degenerate themselves

	checkActiveJacobianRows(m, dd, f)
	checkActiveJacobianEntries(m, dd, f)

	flush(f)

# Step 4: Identify candidate equations

	cand, ids = findCandidates(dd, ds, true)
	println(f,"Identified ",length(cand)," candidate constraints/bounds.")

	if(length(cand) > 0)

# Check for each candidate
		println(f,"Degenerate set from candidate search:")
		printIDS(ids, m, dd, ds, f)

		println(f,"")
		println(f,"nVar = ",string(dd.nVarActive))
		println(f,"nLambda = ",string(dd.nLambda))
		println(f,"\t Constraints = ",string(length(dd.gMap)))
		println(f,"\t Variable Bounds = ",string(length(dd.bMap)))
		println(f,"")

		# This is probably not the most efficient approach.
		# To Do: Extract row and column index vectors from J_active
		# and try working with those
		# Where is this used?
		# sparsityPattern = J_active .!= 0

		if(!ds.onlyCandidateSearch)

			sets = Array(IrreducibleDegenerateSet,length(cand))

			print(f,"Setting up MILP... ")
			tic()
			dh = setupMILP(dd, ds)
			tm = toq()
			println(f,string(tm," seconds"))

			for c = 1:length(cand)
				sets[c] = solveMILP(dh, dd, ds, c)
				printIDS(sets[c], m, dd, ds, f);
			end
		
		else
			sets = ids
		end

	else
		sets = Array(IrreducibleDegenerateSet,0)
	end

	println(f," ")
	println(f,"Leaving Degeneracy Hunter.")
	println(f,"***************************")
	println(f," ")

	return sets, dd, ds

end

function setupMILP(dd::DegenData, ds::DegenSettings)

	M = ds.lambdaM

	#m2 = Model(solver=GurobiSolver(IntFeasTol=1E-6, NumericFocus=3))
	m2 = Model()
	L = 1:dd.nLambda
	@variable(m2,y[L],Bin)
	@variable(m2, -M <= lambda[L] <= M)

	@constraint(m2, degen[j=1:dd.nVarActive], sum{dd.J_active[i,j]*lambda[i], i=L} == 0)


	@constraint(m2, lower[i=L], -M*y[i] <= lambda[i])
	@constraint(m2, upper[i=L],  lambda[i] <= M*y[i])

	@objective(m2, Min, sum{y[i],i=L})

	return m2

end

function solveMILP(m2::Model, dd::DegenData, ds::DegenSettings, c::Int64)

	M = ds.lambdaM

	lambda = getvariable(m2, :lambda)
	L = 1:dd.nLambda

	for i = L
		setlowerbound(lambda[i],-M)
		setupperbound(lambda[i],M)
	end

	setupperbound(lambda[c],1.0)
	setlowerbound(lambda[c],1.0)

	status = solve(m2)

	if(status == :Optimal)

		lmbd_ = getvalue(lambda)
		lmbd = zeros(dd.nLambda)
		for i = L
			lmbd[i] = lmbd_[i]
		end

		y = getvariable(m2, :y)

		y_ = getvalue(y)
		ySln = zeros(dd.nLambda)
		for i = L
			ySln[i] = y_[i]
		end

		return IrreducibleDegenerateSet(c, sum(ySln), lmbd, ySln)

	else

		return IrreducibleDegenerateSet(c, 0, [], [])

	end

end

function findCandidates(dd::DegenData, ds::DegenSettings, explicit::Bool)

	if(explicit)
		mSmall = 1E-4
	else
		mSmall = 1.0
	end

	M = ds.lambdaM

	#m2 = Model(solver=GurobiSolver(IntFeasTol=1E-6, NumericFocus=3))
	m2 = Model()
	L = 1:dd.nLambda
	@variable(m2, -M - mSmall <= lambda[L] <= M + mSmall)
	@variable(m2, yPos[i=L], Bin)
	@variable(m2, yNeg[i=L], Bin)
	@constraint(m2, logic[i=L], yPos[i] + yNeg[i] <= 1)

	if(explicit)

		target = 1

		@variable(m2, ySelect[i=L], Bin)
		@constraint(m2, selectLow[i=L], -(M+target)*(1-ySelect[i]) + target <= lambda[i])
		@constraint(m2, selectHigh[i=L], lambda[i] <= (M-target)*(1-ySelect[i]) + target)

		@constraint(m2, positiveUp[i=L], lambda[i] <= M*yPos[i])
		@constraint(m2, positiveLow[i=L], mSmall*yPos[i] - M*(1-yPos[i]) <= lambda[i])

		@constraint(m2, negativeLow[i=L], -M*yNeg[i] <= lambda[i])
		@constraint(m2, negativeUp[i=L], lambda[i] <= -mSmall*yNeg[i] + M*(1-yNeg[i]))

		@constraint(m2, sum{ySelect[i], i=L} >= 1)

		@objective(m2, Min, sum{yPos[i] + yNeg[i], i=L})

	else
		@variable(m2, 0 <= bound[i=L] <= M + mSmall)


		@constraint(m2, positiveUp[i=L], lambda[i] <= M*yPos[i] + mSmall)
		@constraint(m2, negativeLow[i=L], -M*yNeg[i] - mSmall <= lambda[i])

		# When yPos = 1, lambda >= mSmall
		# When yPos = 0, lambda >= -mSmall
		@constraint(m2, positiveLow[i=L], -mSmall + 2*mSmall*yPos[i] <= lambda[i])

		# When yNeg = 1, lambda <= -mSmall
		# When yNeg = 0, lambda <= mSmall
		@constraint(m2, negativeUp[i=L], lambda[i] <= mSmall - 2*mSmall*yNeg[i])


		@constraint(m2, lower[i=L], -bound[i] <= lambda[i])
		@constraint(m2, upper[i=L], lambda[i] <= bound[i])

		# At least one constraint must be in the degenerate set
		@constraint(m2, sum{yPos[i] + yNeg[i], i=L} >= 1)

		@objective(m2, Min, sum{bound[i], i=L})
	end

	@constraint(m2, degen[j=1:dd.nVarActive], sum{dd.J_active[i,j]*lambda[i], i=L} == 0)

	status = solve(m2)

	c = Array(Int64,0)

	if(status == :Optimal)

		lmbd_ = getvalue(lambda)
		lmbd = zeros(dd.nLambda)
		for i = L
			lmbd[i] = lmbd_[i]
		end

		abs_lmbd = abs(lmbd)

		lmax = 0
		for i = L
			if(lmax < abs_lmbd[i])
				lmax = abs_lmbd[i]
			end
		end

		for i = L
			if(abs_lmbd[i] > 0.1*lmax)
				push!(c,i)
			end
		end

		y = 1.0*(abs_lmbd .> 1/M)

		ids = IrreducibleDegenerateSet(0, sum(y), lmbd, y)

	else
		ids = IrreducibleDegenerateSet(0, 0.0, zeros(0), zeros(0))
	end

	return c, ids

end

############################
# ProblemStats

type ProblemStats
	nCont::Int64
	nContFixed::Int64
	nContFree::Int64
	nContLowerOnly::Int64
	nContUpperOnly::Int64
	nContTwoBounds::Int64
	nContUnbounded::Int64
	nContActiveBounds::Int64
	nBin::Int64
	nBinFixed::Int64
	nBinFree::Int64

	linear::Int64
	linEql::Int64
	linIneqlAct::Int64
	linIneqlInact::Int64
	quad::Int64
	quadEql::Int64
	quadIneqlAct::Int64
	quadIneqlInact::Int64
	nonlinear::Int64
	nonlinEql::Int64
	nonlinIneqlAct::Int64
	nonlinIneqlInact::Int64

	status::Symbol
	solveTime::Float64
end

function ProblemStats(m::Model, status::Symbol=:Unsolved, solveTime::Float64=0.0)
	return assembleProblemStats(m, status, solveTime)
end

function assembleProblemStats(m::Model, status=:Unsolved, solveTime=0.0)

	epsilon = 1E-6

	fixed = m.colUpper - m.colLower .< epsilon
	cont = m.colCat .== :Cont
	lwr = m.colLower .> -Inf
	upr = m.colUpper .< Inf
	bnry = m.colCat .== :Bin

	n = MathProgBase.numconstr(m)

	nLin = MathProgBase.numlinconstr(m)
	nQuad = MathProgBase.numquadconstr(m)
	nNonLin = n - nLin - nQuad

	(gLB, gUB) = JuMP.constraintbounds(m)

	# This portion of code converts m into an NLP
	d = JuMP.NLPEvaluator(m)
	MathProgBase.initialize(d, [:ExprGraph, :Jac])
	
	x = get_x(m, status)

	g = zeros(n,1)
	MathProgBase.eval_g(d,g,x)

	eql = (gUB - gLB) .< epsilon
	ineql = !eql
	active = ((gUB - g) .< epsilon) | ((g - gLB) .< epsilon)

	N = 1:n

	lin = N .<= nLin
	quad = (N .> nLin) & (N .<= nLin + nQuad)
	nonlin = N .> nLin + nQuad

	activeBnds = ((x - m.colLower) .< epsilon) | ((m.colUpper - x) .< epsilon)

	return ProblemStats( sum(cont), # Continuous
						sum(cont & fixed), # Fixed
						sum(cont & !fixed), # Free
						sum(cont & lwr & !upr), # Only lower bounded
						sum(cont & !lwr & upr), # Only upper bounded
						sum(cont & lwr & upr & !fixed), # Lower and upper bounded
						sum(cont & !lwr & !upr), # Unbounded
						sum(activeBnds & !fixed), # Number of active bounds
						sum(bnry), # Binary
						sum(bnry & fixed), # Binary fixed
						sum(bnry & !fixed), # Binary free
						nLin, # Linear constraints
						sum(lin & eql),
						sum(lin & ineql & active),
						sum(lin & ineql & !active),
						nQuad, # Quadratic constraints
						sum(quad & eql),
						sum(quad & ineql & active),
						sum(quad & ineql & !active),
						nNonLin, # Nonlinear constraints
						sum(nonlin & eql),
						sum(nonlin & ineql & active),
						sum(nonlin & ineql & !active),
						status, # Solver status
						solveTime # Solver time
						)

end

function stringProblemStats(ps::ProblemStats)

	return join([ps.nCont, ps.nContFixed, ps.nContFree, ps.nContLowerOnly, ps.nContUpperOnly,ps.nContTwoBounds,ps.nContUnbounded,ps.nContActiveBounds,
					ps.nBin,ps.nBinFixed,ps.nBinFree,
					ps.linear,ps.linEql,ps.linIneqlAct,ps.linIneqlInact,
					ps.quad,ps.quadEql,ps.quadIneqlAct,ps.quadIneqlInact,
					ps.nonlinear,ps.nonlinEql,ps.nonlinIneqlAct,ps.nonlinIneqlInact,
					string(ps.status),ps.solveTime],",")

end

function stringProblemStats()

	return string("nContinous,nContFixed,nContFree,nContLowrOnly,nContUpprOnly,nContBndBoth,nContUnbnd,nContActiveBnds,nBinary,nBinFixed,nBinFree,nLinear,nLinEql,nLinIneqlAct,nLinIneqlInact,nQuadratic,nQuadEql,nQuadIneqlAct,nQuadIneqlInact,nNonlinear,nNLEql,nNLIneqlAct,nNLIneqlInact,status,solveTime")

end

function printProblemStats(ps::ProblemStats, f=STDOUT)

	println(f," ")
	println(f,"Continous Variables (total): ",ps.nCont)
	println(f,"\t Fixed: ",ps.nContFixed)
	println(f,"\t Free: ",ps.nContFree)
	println(f,"\t\t Lower Bounded Only: ",ps.nContLowerOnly)
	println(f,"\t\t Upper Bounded Only: ",ps.nContUpperOnly)
	println(f,"\t\t Lower \& Upper Bounded: ",ps.nContTwoBounds)
	println(f,"\t\t Unbounded: ",ps.nContUnbounded)
	println(f,"\t\t Active Bounds: ",ps.nContActiveBounds)
	println(f," ")
	println(f,"Binary Variable (total): ",ps.nBin)
	println(f,"\t Fixed: ",ps.nBinFixed)
	println(f,"\t Free: ",ps.nBinFree)
	println(f," ")

	println(f,"Linear Constraints: ",ps.linear)
	println(f,"\t Equality: ",ps.linEql)
	println(f,"\t Ineq. (Active): ",ps.linIneqlAct)
	println(f,"\t Ineq. (Inactive): ",ps.linIneqlInact)
	println(f," ")

	println(f,"Quadratic Constraints: ",ps.quad)
	println(f,"\t Equality: ",ps.quadEql)
	println(f,"\t Ineq. (Active): ",ps.quadIneqlAct)
	println(f,"\t Ineq. (Inactive): ",ps.quadIneqlInact)
	println(f," ")

	println(f,"Nonlinear Constraints: ",ps.nonlinear)
	println(f,"\t Equality: ",ps.nonlinEql)
	println(f,"\t Ineq. (Active): ",ps.nonlinIneqlAct)
	println(f,"\t Ineq. (Inactive): ",ps.nonlinIneqlInact)
	println(f," ")

	println(f,"Solver Exit Status: ",string(ps.status))
	println(f,"Solve Time: ",string(ps.solveTime), "s")
	println(f," ")

	return nothing

end

################
# Hidden Auxiliary Functions

function checkJacobianScaling(m::Model, dd::DegenData, f=STDOUT)

	# Determine location of non-zeros. This is important as fixed variables
	# are still in the Jacobian
	nz_ = find(dd.J .!= 0)

	# Determine location of smallest and largest elements
	s_ = nz_[indmin(abs(dd.J[nz_]))]
	l_ = indmax(abs(dd.J))


	println(f," ")
	println(f,"Smallest non-zero element in Jacobian = ",dd.J[s_])
	println(f,"Variable: ")
	printVariable(m, dd.jC[s_], f)
	println(f,"Equation: ")
	printEquation(m, dd, dd.iR[s_], f)

	println(f," ")
	println(f,"Largest non-zero element in Jacobian = ",dd.J[l_])
	println(f,"Variable: ")
	printVariable(m, dd.jC[l_], f)
	println(f,"Equation: ")
	printEquation(m, dd, dd.iR[l_], f)
	println(f," ")

end

function checkActiveJacobianRows(m::Model, dd::DegenData, f=STDOUT)
	for i = 1:dd.nLambda
		if(sum(abs(dd.J_active[i,:])) < 1E-4)
			println(f,"Warning!!! (Near) singular Jacobian row for ")
			printRows(m, dd, i, f)
			println(f," ")
		end
	end

	return nothing
end

function checkActiveJacobianEntries(m::Model, dd::DegenData, f=STDOUT)
	for i = 1:size(dd.J_active,1)
		nan_ = sum(isnan(dd.J_active[i,:])) > 0
		inf_ = sum(isinf(dd.J_active[i,:])) > 0
		if(nan_ || inf_)
			if(nan_)
				print(f,"NaN")
			end
			if(nan_ && inf_)
				print(f," and ")
			end
			if(inf_)
				print(f,"Inf")
			end
			println(f," in Jacobian")
			printJacobianForEquation(m, dd, f)
			println(f," ")
		end
	end

	return nothing
end

function printEquation(m::Model, dd::DegenData, i::Int64, f=STDOUT)

	nLin = MathProgBase.numlinconstr(m)
	nQuad = MathProgBase.numquadconstr(m)

	if(i <= nLin)
		#println(f,ConstraintRef{LinearConstraint}(m, i))
		println(f, m.linconstr[i])
	elseif(i <= nLin + nQuad)
		#println(f,ConstraintRef{QuadConstraint}(m, i - nLin))
		println(f, m.quadconstr[i - nLin])
	else
		println(f,MathProgBase.constr_expr(dd.d,i))
	end

	return nothing

end

#=
function printRow(m::Model, dd::DegenData, j::Int64, epsiActive::Float64, f=STDOUT)

	if(j <= length(dd.gMap))
		i = dd.gMap[j]

		printEquation(m, dd, i, f)

	else
		i = dd.bMap[j - length(dd.gMap)]

		printBound(m, i, epsiActive, f)

	end

	return nothing

end
=#

function printVariable(m::Model, dd::DegenData, i::Int64, f=STDOUT)

	v = Variable(m, dd.vMap[i])
	println(f,"x[",string(i),"] = ",getname(v))

	return nothing

end

function printVariable(m::Model, i::Int64, f=STDOUT)
	v = Variable(m, i)
	println(f,"x[",string(i),"] = ",getname(v))
	return nothing
end

function printJacobianForEquation(m::Model, dd::DegenData, r::Int64, f = STDOUT)

	printEquation(m, dd, f)

	if(dd.Jactive == Void)
		# Use sparse tuple representation of full Jacobian

		k = find(r == dd.iR)
		for i=1:length(k)
			J_ = dd.J[k[i]]
			if(isnan(J_) || J_ != 0)
				c = jC[k]
				v = Variable(m,c)
				println(f,"J = ",J_," for x[",c,"] = ",getname(v)," = ",dd.x[c])
			end
		end

	else
		# Sparse matrix representation of active Jacobian
		for c = 1:size(dd.Jactive,2)
			if(isnan(dd.Jactive[r,c]) || dd.Jactive[r,c] != 0)
				i = vMap[c]
				v = Variable(m, i)
				println(f,"J = ",Jactive[r,c]," for x[",i,"] = ",getname(v)," = ",dd.x[i])
			end
		end

	end

	return nothing

end

function printVariablesInEquation(m::Model, r::Int64)

	dd = DegenData(m)
	printVariablesInEquation(m, dd, r, true)

end

function printVariablesInEquation(m::Model, dd::DegenData, r::Int64, printVarNumber::Bool, printVarBounds::Bool=true, f=STDOUT)
	printVariablesInEquation(m, dd, Array([r]), printVarNumber, printVarBounds, f)
end

function printVariablesInEquation(m::Model, dd::DegenData, r::Array{Int64,1}, printVarNumber::Bool, printVarBounds::Bool=true, f=STDOUT)

	if(length(r) > 1)
		k = Array(Int64,0)

		for i = 1:length(r)
			t = find(dd.iR .== r[i])
			for j = 1:length(t)
				push!(k, t[j])
			end
		end

	else
		k = find(dd.iR .== r)
	end

	if(length(unique(dd.jC[k])) < 50)
		for j in unique(dd.jC[k])
			v = Variable(m,j)

			if(printVarNumber)
				print(f,"x[",string(j),"] = ")
			end

			print(f,getname(v)," = ",dd.x[j])

			if(printVarBounds)
				print(f,"  \t [")
				lo = getlowerbound(v)
				up = getupperbound(v)

				if(abs(dd.x[j] - lo) < 1E-4)
					print_with_color(:red,f, string(lo))
				else
					print(f,string(lo))
				end

				print(f,",")

				if(abs(dd.x[j] - up) < 1E-4)
					print_with_color(:blue,f, string(up))
				else
					print(f,string(up))
				end

				print(f,"]")
			end

			print(f,"\n")

		end
	else
		println(f,"Over 50 variables in this equation... not printing variable values")
	end

	return nothing

end

function printIDS(ids::IrreducibleDegenerateSet, m::Model, dd::DegenData, ds::DegenSettings, f=STDOUT)

	epsiLambda = ds.epsiLambda
	epsiActive = ds.epsiActive

	if(ids.nElements > 0)

		println(f,"IDS ",string(ids.c),"...")

		nLambda = length(ids.lambda)

		r = find(abs(ids.lambda) .> epsiLambda)

		for j = r

			s = @sprintf "%.4f" ids.y[j]
			print(f,"y = ", s,", ")

			s = @sprintf "%.4f" ids.lambda[j]
			print(f,"l = ", s ,"  \t ")

			printRows(m, dd, j, f, lite=true)
		end

		println(f," ")

		println(f,"Involved variables: ")

		r2 = intersect(r, 1:length(dd.gMap))
		printVariablesInEquation(m, dd, dd.gMap[r2], true, true, f)

		println(f," ")


	else

		#println("The following "," does not contribute significantly to a IDS")
	end

	return nothing

end

end
