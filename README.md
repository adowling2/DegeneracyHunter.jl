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

function example()

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
```

###