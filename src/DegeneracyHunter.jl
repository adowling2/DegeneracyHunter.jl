##### DegeneracyHunter.jl
# Created by Alexander Dowling

module DegeneracyHunter

export IrreducibleDegenerateSet, DegenSettings, DegenData, printBound, printVariablesInEquation,
	printIDS, degeneracyHunter, printInfeasibleEquations, printInactiveEquations,
	checkVarBounds, printRows, ProblemStats, assembleProblemStats, stringProblemStats, printProblemStats,
	printVariableDiagnostics, round_to_bounds!, checkEquationScaling

importall JuMP
import MathProgBase

include("types.jl")
include("model_diagnostics.jl")
include("classic_degeneracy_hunter.jl")
include("problem_stats.jl")
include("auxillary.jl")

end # module
