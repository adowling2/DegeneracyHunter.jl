##### auxillary.jl
# Created by Alexander Dowling

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
	printVariablesInEquation(m, dd, dd.iR[s_], true, true, f)

	println(f," ")
	println(f,"Largest non-zero element in Jacobian = ",dd.J[l_])
	println(f,"Variable: ")
	printVariable(m, dd.jC[l_], f)
	println(f,"Equation: ")
	printEquation(m, dd, dd.iR[l_], f)
	printVariablesInEquation(m, dd, dd.iR[l_], true, true,  f)
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
