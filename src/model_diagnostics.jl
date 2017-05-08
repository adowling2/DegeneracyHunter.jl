##### model_diagnostics.jl
# Created by Alexander Dowling

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