# DegeneracyHunter.jl

Creating good optimization models is difficult. Degeneracy Hunter is a collection of utilities to help JuMP/Julia model development. Key features include:
* Initial point analysis: display uninitialized variables and values violating bounds 
* Print infeasible constraint residuals
* Print inactive inequality constraints
* Determine Irreducible Set of Infeasible Equations (original *Degeneracy Hunter*)

## Degeneracy Hunter Algorithms

Degeneracy Hunter seeks to find Irreducible Degenerate Sets (IDS) by analyzing the Jacobian of the constraints at a specified point. This helps modelers identify problems, such as redundant constraints, overspecifications, etc., which all impact solver reliability. Additional details are available in the following conference paper:

Alexander W. Dowling and Lorenz T. Biegler (2015). "Degeneracy Hunter: An Algorithm for Determining Irreducible Sets of Degenerate Constraints in Mathematical Programs". 12th International Symposium on Process Systems Engineering and 25th European Symposium on Computer Aided Process Engineering. Ed. by Krist V. Gernaey, Jakob K. Huusom, and Raqul Gani. Computer-Aided Chemical Engineering, Vol. 37, p. 809 - 814. [link](http://www.sciencedirect.com/science/article/pii/B9780444635785501304) [pdf](http://sites.google.wisc.edu/alex-dowling/uploads/ESCAPE25_Dowling.pdf?attredirects=0)

## Example

```
using JuMP

function lp1()

	m = Model()

	@variable(m, 0 <= x[1:3] <= 5)
	@variable(m, 0 <= y <= 0)
	@constraint(m, x[1] + x[2] >= 1)
	@constraint(m, x[1] + x[2] + x[3] == 1)
	@constraint(m, x[2] - 2*x[3] <= 1)
	@constraint(m, x[1] + x[3] >= 1)
	@constraint(m, x[1] + x[2] + x[3] == 1) # Redundant constraint - makes problem degenerate
	@objective(m, Min, sum{x[i],i=1:3})
	
	return m
end

function initialize!(m::Model)
	x = getvariable(m,:x)
	setvalue(x[1], 1.0)
	setvalue(x[2], 5.0)
	setvalue(x[3], -1.0)
	
	# Intentionally do not initialize y
	
	return nothing
end
```

### Create and initialize model
```
m = lp1()
initialize!(m)
print(m)
```

**Output:**
```
Min x[1] + x[2] + x[3]
Subject to
 x[1] + x[2] ≥ 1
 x[1] + x[2] + x[3] = 1
 x[2] - 2 x[3] ≤ 1
 x[1] + x[3] ≥ 1
 x[1] + x[2] + x[3] = 1
 0 ≤ x[i] ≤ 5 ∀ i ∈ {1,2,3}
 0 ≤ y ≤ 0
```

### Check initial point
```
DegeneracyHunter.printVariableDiagnostics(m)
```

**Output:**
```
Uninitialized Variables: 
y
 
Variable Violating Bounds: 
x[3] = -1.0, bounds: [ 0.0 , 5.0 ]
```

### Check for infeasible equations
Print equations with absolute residuals greater than 0.001
```
DegeneracyHunter.printInfeasibleEquations(m, 0.001)
```

**Output:**
```
Infeasible equations: 
r[4] = -1.0
x[1] + x[3] ≥ 1
x[1] = x[1] = 1.0  	 [0.0,5.0]
x[3] = x[3] = -1.0  	 [0.0,5.0]
 
r[2] = 4.0
x[1] + x[2] + x[3] = 1
x[1] = x[1] = 1.0  	 [0.0,5.0]
x[2] = x[2] = 5.0  	 [0.0,5.0]
x[3] = x[3] = -1.0  	 [0.0,5.0]
 
r[5] = 4.0
x[1] + x[2] + x[3] = 1
x[1] = x[1] = 1.0  	 [0.0,5.0]
x[2] = x[2] = 5.0  	 [0.0,5.0]
x[3] = x[3] = -1.0  	 [0.0,5.0]
 
r[3] = 6.0
x[2] - 2 x[3] ≤ 1
x[2] = x[2] = 5.0  	 [0.0,5.0]
x[3] = x[3] = -1.0  	 [0.0,5.0]
```

### Check for degenerate constraints
Use default settings
```
DegeneracyHunter.degeneracyHunter(m)
```

**Output:**
```
******************************************
Welcome to Degeneracy Hunter!
```


Option values used:
```
Options: 
includeBounds = true
includeWeaklyActive = true
removeFixedVar = true
epsiActive = 1.0e-6
epsiLambda = 1.0e-6
lambdaM = 100000.0
```

Smallest and largest (in magnitude) non-zero Jacobian elements.
```
Smallest non-zero element in Jacobian = 1.0
Variable: 
x[1] = x[1]
Equation: 
x[1] + x[2] ≥ 1
 
Largest non-zero element in Jacobian = -2.0
Variable: 
x[3] = x[3]
Equation: 
x[2] - 2 x[3] ≤ 1
```

Information about problem size:
``` 
(Weakly) Active Inequality Constraints: 0
Inactive Inequality Constraints: 3
Equality Constraints: 2
 
Variables: 4
(Weakly) Active Variable Bounds: 4
Fixed Variables: 1
```

Time to complete some behind-the-scenes steps:
```
Adding Jacobian elements for bounds... 4.589e-6 seconds
Assembling J_sparse... 2.844e-6 seconds
Assembling J_active... 4.61e-6 seconds
```

Output from Gurobi. IDS calculation is posed as a MILP.
```
Optimize a model with 39 rows, 20 columns and 84 nonzeros
Coefficient statistics:
  Matrix range    [1e+00, 1e+05]
  Objective range [1e+00, 1e+00]
  Bounds range    [1e+00, 1e+05]
  RHS range       [1e+00, 1e+05]
Found heuristic solution: objective 5
Presolve time: 0.00s
Presolved: 39 rows, 20 columns, 134 nonzeros
Variable types: 5 continuous, 15 integer (15 binary)

Root relaxation: objective 1.000003e+00, 17 iterations, 0.00 seconds

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0    1.00000    0    7    5.00000    1.00000  80.0%     -    0s
H    0     0                       2.0000000    1.00000  50.0%     -    0s

Explored 0 nodes (21 simplex iterations) in 0.00 seconds
Thread count was 2 (of 2 available processors)

Optimal solution found (tolerance 1.00e-04)
Best objective 2.000000000000e+00, best bound 2.000000000000e+00, gap 0.0%
```

```
Identified 2 candidate constraints/bounds.
```

IDS information from initial serach. The values for ``l`` show how the equations are added together to create rank deficiency.
```
Degenerate set from candidate search:
IDS 0...
y = 1.0000, l = -1.0000  	 **********************************
Constraint 1...
x[1] + x[2] ≥ 1
1.0 <= 0.0 <= Inf
 
Involved Variables:
x[1] 	0.0 <= 1.0 <= 5.0
x[2] 	0.0 <= 5.0 <= 5.0
 
y = 1.0000, l = 1.0000  	 **********************************
Constraint 2...
x[1] + x[2] + x[3] = 1
1.0 <= 0.0 <= 1.0
 
Involved Variables:
x[1] 	0.0 <= 1.0 <= 5.0
x[2] 	0.0 <= 5.0 <= 5.0
x[3] 	0.0 <= -1.0 <= 5.0
 
 
Involved variables: 
x[1] = x[1] = 0.0  	 [0.0,5.0]
x[2] = x[2] = 0.0  	 [0.0,5.0]
x[3] = x[3] = 0.0  	 [0.0,5.0]
```

Information on the MILP problem size for additional searches:
```
nVar = 3
nLambda = 5
	 Constraints = 2
	 Variable Bounds = 3

Setting up MILP... 9.6227e-5 seconds
```

For each candidate equation, attempt to compute an IDS where l = 1 for said equation. If the MILP is feasible, an IDS exists containing this equation.

Gurobi results for candidate 1:
```
Optimize a model with 13 rows, 10 columns and 29 nonzeros
Coefficient statistics:
  Matrix range    [1e+00, 1e+05]
  Objective range [1e+00, 1e+00]
  Bounds range    [1e+00, 1e+05]
  RHS range       [0e+00, 0e+00]
Found heuristic solution: objective 5
Presolve removed 13 rows and 10 columns
Presolve time: 0.00s

Explored 0 nodes (0 simplex iterations) in 0.00 seconds
Thread count was 1 (of 2 available processors)

Optimal solution found (tolerance 1.00e-04)
```

IDS containing candidate 1:
```
Best objective 2.000000000000e+00, best bound 2.000000000000e+00, gap 0.0%
IDS 1...
y = 1.0000, l = 1.0000  	 **********************************
Constraint 1...
x[1] + x[2] ≥ 1
1.0 <= 0.0 <= Inf
 
Involved Variables:
x[1] 	0.0 <= 1.0 <= 5.0
x[2] 	0.0 <= 5.0 <= 5.0
 
y = 1.0000, l = -1.0000  	 **********************************
Constraint 2...
x[1] + x[2] + x[3] = 1
1.0 <= 0.0 <= 1.0
 
Involved Variables:
x[1] 	0.0 <= 1.0 <= 5.0
x[2] 	0.0 <= 5.0 <= 5.0
x[3] 	0.0 <= -1.0 <= 5.0
 
 
Involved variables: 
x[1] = x[1] = 0.0  	 [0.0,5.0]
x[2] = x[2] = 0.0  	 [0.0,5.0]
x[3] = x[3] = 0.0  	 [0.0,5.0]
```

Repeat for candidate 2:
```
Optimize a model with 13 rows, 10 columns and 29 nonzeros
Coefficient statistics:
  Matrix range    [1e+00, 1e+05]
  Objective range [1e+00, 1e+00]
  Bounds range    [1e+00, 1e+05]
  RHS range       [0e+00, 0e+00]
Found heuristic solution: objective 5
Presolve removed 13 rows and 10 columns
Presolve time: 0.00s

Explored 0 nodes (0 simplex iterations) in 0.00 seconds
Thread count was 1 (of 2 available processors)

Optimal solution found (tolerance 1.00e-04)
Best objective 2.000000000000e+00, best bound 2.000000000000e+00, gap 0.0%
```

Print IDS containing candidate 2:
```
IDS 2...
y = 1.0000, l = -1.0000  	 **********************************
Constraint 1...
x[1] + x[2] ≥ 1
1.0 <= 0.0 <= Inf
 
Involved Variables:
x[1] 	0.0 <= 1.0 <= 5.0
x[2] 	0.0 <= 5.0 <= 5.0
 
y = 1.0000, l = 1.0000  	 **********************************
Constraint 2...
x[1] + x[2] + x[3] = 1
1.0 <= 0.0 <= 1.0
 
Involved Variables:
x[1] 	0.0 <= 1.0 <= 5.0
x[2] 	0.0 <= 5.0 <= 5.0
x[3] 	0.0 <= -1.0 <= 5.0
 
 
Involved variables: 
x[1] = x[1] = 0.0  	 [0.0,5.0]
x[2] = x[2] = 0.0  	 [0.0,5.0]
x[3] = x[3] = 0.0  	 [0.0,5.0]
```

Finished iterative over candidate equations. Done.
```
Leaving Degeneracy Hunter.
***************************
```

### Solve model
solve(m)

### Print inactive constraints
```
DegeneracyHunter.degeneracyHunter(m)
```

**Output:**
```
Inactive Equations: 
x[2] - 2 x[3] ≤ 1
```

### Print problem statistics
```
ps = DegeneracyHunter.assembleProblemStats(m,status,tm)
DegeneracyHunter.printProblemStats(ps)
```

**Output:**
```
Continous Variables (total): 4
	 Fixed: 1
	 Free: 3
		 Lower Bounded Only: 0
		 Upper Bounded Only: 0
		 Lower & Upper Bounded: 3
		 Unbounded: 0
		 Active Bounds: 2
 
Binary Variable (total): 0
	 Fixed: 0
	 Free: 0
 
Linear Constraints: 5
	 Equality: 2
	 Ineq. (Active): 2
	 Ineq. (Inactive): 1
 
Quadratic Constraints: 0
	 Equality: 0
	 Ineq. (Active): 0
	 Ineq. (Inactive): 0
 
Nonlinear Constraints: 0
	 Equality: 0
	 Ineq. (Active): 0
	 Ineq. (Inactive): 0
 
Solver Exit Status: Optimal
Solve Time: 0.000768846s
```
