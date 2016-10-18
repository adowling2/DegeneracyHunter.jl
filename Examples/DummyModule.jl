module DummyModule

using JuMP

export solveOptProb2, MyResults

function solveOptProb2(p)
	a = p[1]
	b = p[2]
	m = Model()
	@variable(m, 0 <= x <= 2 )
	@variable(m, 0 <= y <= 30 )

	@objective(m, Max, a*x + 3*y )
	@constraint(m, 1x + 5y <= b )
	
	status = solve(m)
	
	return MyResults(getobjectivevalue(m), getvalue(x), getvalue(y), status)
	
end

type MyResults
	obj::Float64
	x::Float64
	y::Float64
	status::Symbol
end

function print(r::MyResults)

	println(" ")
	println("objective = ",r.obj)
	println("x = ",r.x,"\t y = ",r.y)
	println(r.status)

	return nothing

end

end