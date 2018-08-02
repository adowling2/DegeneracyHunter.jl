##### classic_degeneracy_hunter.jl
# Created by Alexander Dowling

##############
##### Degeneracy Hunter features

function degeneracyHunter(m::Model, mySolver;
	includeBounds::Bool = false,
	includeWeaklyActive::Bool = true,
	removeFixedVar::Bool = true,
	onlyCandidateSearch::Bool = true,
	epsiActive::Float64 = 1E-6,
	epsiLambda::Float64 = 1E-6,
	lambdaM::Float64 = 1E5,
	f=STDOUT)

	ds = DegenSettings(includeBounds, includeWeaklyActive, removeFixedVar,
		onlyCandidateSearch, epsiActive, epsiLambda, lambdaM, mySolver)

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
	gLow = (abs.(dd.gLB - dd.g) .< ds.epsiActive)
	gUp = (abs.(dd.gUB - dd.g) .< ds.epsiActive)
	gEql = (abs.(dd.gUB - dd.gLB) .< ds.epsiActive)

	gActive = gLow .| gUp .| gEql

	# Tip: "$" was the bitwise xor operator; ⊻ is the new operator.
	actIneql = sum(xor.(gLow, gUp))
	eql = sum(gEql)
	totalIneql = length(dd.g) - eql
	inactiveIneql = totalIneql - actIneql

	println(f," ")
	println(f,"(Weakly) Active Inequality Constraints: ",string(actIneql))
	println(f,"Inactive Inequality Constraints: ",string(inactiveIneql))
	println(f,"Equality Constraints: ",string(eql))
	println(f," ")


	if(ds.includeWeaklyActive)
		# Do nothing
	else
		# Query multipliers

		# Need to implement

	end

	vFixed = abs.(m.colLower - m.colUpper) .< ds.epsiActive

	bActive = (abs.(m.colLower - dd.x) .< ds.epsiActive) .| (abs.(m.colUpper - dd.x) .< ds.epsiActive)

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
		bActive = bActive .& .!vFixed
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
		dd.vMap = find(.!vFixed)
		cols = .!vFixed
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
# 			of if any Jacobian entries are NaN or Inf

	checkActiveJacobianRows(m, dd, f)
	checkActiveJacobianEntries(m, dd, f)

	flush(f)

# Step 3.5: Perform SVD.

	print("Computing SVD... ")
	tic()
	# This is ugly, as it converts dd.J_active from a sparse to dense matrix for SVD.
	# TODO: Use a sparse SVD implementation.
	F = LinAlg.svdfact(Matrix(dd.J_active))
	tm = toq()
	println(f, string(tm," seconds"))
	
	println(" ")
	println(f,"Largest Singular Value: ",F[:S][1])
	println(f,"Smallest Singular Value: ",F[:S][end])
	
	println(" ")

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
		sets = Array{IrreducibleDegenerateSet}(0)
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
	m2 = Model(solver=mySolver)
	L = 1:dd.nLambda
	@variable(m2,y[L],Bin)
	@variable(m2, -M <= lambda[L] <= M)

	@constraint(m2, degen[j=1:dd.nVarActive], sum(dd.J_active[i,j]*lambda[i] for i in L) == 0)


	@constraint(m2, lower[i=L], -M*y[i] <= lambda[i])
	@constraint(m2, upper[i=L],  lambda[i] <= M*y[i])

	@objective(m2, Min, sum(y[i] for i in L))

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

	#m2 = Model(solver=Solver(IntFeasTol=1E-6, NumericFocus=3))
	m2 = Model(solver=ds.mySolver)
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

		@constraint(m2, sum(ySelect[i] for i in L) >= 1)

		@objective(m2, Min, sum(yPos[i] + yNeg[i] for i in L))

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
		@constraint(m2, sum(yPos[i] + yNeg[i] for i in L) >= 1)

		@objective(m2, Min, sum(bound[i] for i in L))
	end

	@constraint(m2, degen[j=1:dd.nVarActive], sum(dd.J_active[i,j]*lambda[i] for i in L) == 0)

	status = solve(m2)

	c = Array{Int64}(0)

	if(status == :Optimal)

		lmbd_ = getvalue(lambda)
		lmbd = zeros(dd.nLambda)
		for i = L
			lmbd[i] = lmbd_[i]
		end

		abs_lmbd = abs.(lmbd)

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
