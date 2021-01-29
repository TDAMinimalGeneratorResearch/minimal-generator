"""
Functions for finding the optimal uniform weighted generators and
	optimal length weighted minimal generators. Implemented based
	on paper: "Optimal Cycles for Persistent Homology Via Linear Programming"
	by Emerson G. Escolar and Yasuaki Hiraoka.
"""

function findUniformWeightedOptimalCycles(C, d, l, intSol = false, pivot = true)
	"""
	input:
	C --> Homology Object
	d --> dimension of the generator we want to minimize
	l --> index of the generator we want to minimize
	intSol --> whether require integral solutions or not
	pivot --> using only pivot columns in the optimization problem. (recommended to improve computation time)
	output:
	1. minimal generator vector
	2. the length weighted cost of the optimal generator
	3. the uniform weighted cost of the optimal generator
	"""
	allgens = C.generators[d]
	othergens = C_d_l__minimal(C, 1, l, false, intSol, C.generators[1][1:l-1], pivot)
	a = zeros(length(othergens), 1)
	nonzero = findall(x->x!=0, othergens)
	a[nonzero, 1] .= 1
	len = transpose(C.distVec[C.grainVec[2]]) * a
	return othergens, len[1], length(nonzero)
end

function findLengthWeightedOptimalCycles(C, d, l, intSol = false, pivot = true)
	"""
	input:
	C --> Homology Object
	d --> dimension of the generator we want to minimize
	l --> index of the generator we want to minimize
	intSol --> whether require integral solutions or not
	pivot --> using only pivot columns in the optimization problem.
	output:
	1. minimal generator vector
	2. the length weighted cost of the optimal generator
	3. the uniform weighted cost of the optimal generator
	"""
	allgens = C.generators[d]
	othergens = C_d_l__minimal(C, 1, l, true, intSol, C.generators[1][1:l-1], pivot)
	a = zeros(length(othergens), 1)
	nonzero = findall(x->x!=0, othergens)
	a[nonzero, 1] .= 1
	len = transpose(C.distVec[C.grainVec[2]]) * a
	return othergens, len[1], length(nonzero)
end


function C_d_l__minimal(C,d,l, with_length, int_sol, existing_gens = [], pivot = true)
	"""
	* (C,d,l) 	  --> vector xx which is a minimal generator homologous to the
						lth d-dimensional generator (c_l) of C, unless existing gens
						are given, in which case it is homologous to a sum of
						c_l and existing_gens.
						with_length indicates whether to
						minimize euclidean distance (if false, will instead
						minimize number of simplices) *don't set with_length to
						true for dimension not equal to 1
						setting int_sol = false may be faster, but may produce
						bad solutions sometimes
	"""
	bdrMatrices = C.bdrMatrices
	ddi1 = C_d_l__filt(C,d,l, pivot)
	gen = C_d_l__gen(C,d,l)

    #if given existing generators, include them as cols of bdr matrix
	if length(existing_gens) != 0
		for i in 1 : length(existing_gens)
			ddi1 = hcat(ddi1, existing_gens[i])
		end
	end

	## establish model and varible lengths
	xlen = ddi1.m
	ylen = length(ddi1.colptr) - 1
	m2 = Model(Gurobi.Optimizer) # Clp.Optimizer
	## create variable x which will be our new generator, and y for the boundary i+1
	if int_sol  # --- later, int.  vs. {0, \pm 1}
		println("int sol")
		@variable(m2, x_pos[1:xlen], Int)
		@variable(m2, x_neg[1:xlen], Int)
	else
		println("not int sol")
		@variable(m2, x_pos[1:xlen])
		@variable(m2, x_neg[1:xlen])
	end
	@variable(m2, y[1:ylen], Int)
	## we want to minimize sum abs x such that x differs from our generator by a boundary
	if with_length
		@objective(m2, Min, transpose(C.distVec[C.grainVec[2]])*(x_pos[1:xlen] .+ x_neg[1:xlen]) ) #shortest length
	else # sum(x_pos[1:xlen].+x_neg[1:xlen])) #
		@objective(m2, Min, sum(x_pos[1:xlen].+x_neg[1:xlen]))# length(findall( x-> x!=0 , x_pos[1:xlen] .+ x_neg[1:xlen]))) # for smallest sum of coeff.
	end
	@constraints(m2, begin x_pos[1:xlen] .>= 0; x_neg[1:xlen] .>= 0 end)
	@constraint(m2, x_pos[1:xlen] .- x_neg[1:xlen] .- ddi1 * y .==  Array(gen))
	# print([2].==[2])
	## Optimize the function and show the nonzero values of our vector x
	status = optimize!(m2)
	x = x_pos .- x_neg
	xx = spzeros(length(x))
	for i in 1:xlen
		if abs(JuMP.value(x[i])) >= .001
			xx[i] = JuMP.value(x[i])
		end
	end
	yy = zeros(ylen)
	for i in 1:ylen
		if abs(JuMP.value(y[i])) >= .001
			yy[i] = JuMP.value(y[i])
		end
	end
	return xx
end
