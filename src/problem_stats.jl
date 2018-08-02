##### problem_stats.jl
# Created by Alexander Dowling

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
	ineql = .!eql
	active = ((gUB - g) .< epsilon) .| ((g - gLB) .< epsilon)

	N = 1:n

	lin = N .<= nLin
	quad = (N .> nLin) .& (N .<= nLin + nQuad)
	nonlin = N .> nLin + nQuad

	activeBnds = ((x - m.colLower) .< epsilon) .| ((m.colUpper - x) .< epsilon)

	return ProblemStats( sum(cont), # Continuous
						sum(cont .& fixed), # Fixed
						sum(cont .& .!fixed), # Free
						sum(cont .& lwr .& .!upr), # Only lower bounded
						sum(cont .& .!lwr .& upr), # Only upper bounded
						sum(cont .& lwr .& upr .& .!fixed), # Lower and upper bounded
						sum(cont .& .!lwr .& .!upr), # Unbounded
						sum(activeBnds .& .!fixed), # Number of active bounds
						sum(bnry), # Binary
						sum(bnry .& fixed), # Binary fixed
						sum(bnry .& .!fixed), # Binary free
						nLin, # Linear constraints
						sum(lin .& eql),
						sum(lin .& ineql .& active),
						sum(lin .& ineql .& .!active),
						nQuad, # Quadratic constraints
						sum(quad .& eql),
						sum(quad .& ineql .& active),
						sum(quad .& ineql .& .!active),
						nNonLin, # Nonlinear constraints
						sum(nonlin .& eql),
						sum(nonlin .& ineql .& active),
						sum(nonlin .& ineql .& .!active),
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