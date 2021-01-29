"""
This script implements the algorithm in Escolar's paper on finding uniform-weighted optimal cycles.
"""

"""
Common Input Variable Definitions

	C				 =>		 a Homology Object returned by the computeHomology function.
	d  				 => 	 an integer representing the dimension of the original cycle, e.g. 1 implies H1 homology.
	l 				 =>	   	 an integer representing the index of the original cycle.
	bdrMat 			 =>  	 a sparse matrix representing the boundary matrix with rows indexed by 1-simplices(edges),
					 		 and columns indexed by 2-simplices(triangles).
	integer			 => 	 a boolean, if true, forces the optimal cycle to have integral entries.

Common Output Variable Definitions
	volGens		   	 => 	 a sparse vector representing the area optimal cycle.
	numOfTri	 	 => 	 an integer representing the volume(number of triangles) of the optimal cycle.
	len				 => 	 a float representing the length of the optimal cycle.
	numEdge			 => 	 an integer representing the number of the edges of the optimal cycle.
	volVec			 => 	 a vector representing the coefficients for the triangles making up of the optimal cycle.
	firstTimeSuccess =>   	 The original paper mentions one performance improving trick,
							 this boolean represents whether the optimization indeed saved us time.
"""
function findZeroNormOptimalCycles(C, d, l, intSol = false, verbose = false, pivot = true)
	allgens = C.generators[d]
	othergens, yy, ddi1 = C_d_l__minimal(C, 1, l, false, intSol, C.generators[1][1:l-1], pivot)
	a = zeros(length(othergens), 1)
	nonzero = findall(x->x!=0, othergens)
	a[nonzero, 1] .= 1
	if verbose
		len = transpose(C.distVec[C.grainVec[2]]) * a
		return othergens, len, length(nonzero), yy, ddi1
	else
		return othergens, yy, ddi1
	end
end

function findOneNormOptimalCycles(C, d, l, intSol = false, verbose = false, pivot = true)
	allgens = C.generators[d]
	othergens, yy, ddi1 = C_d_l__minimal(C, 1, l, true, intSol, C.generators[1][1:l-1], pivot)
	a = zeros(length(othergens), 1)
	nonzero = findall(x->x!=0, othergens)
	a[nonzero, 1] .= 1
	if verbose
		len = transpose(C.distVec[C.grainVec[2]]) * a
		return othergens, len, length(nonzero), yy, ddi1
	else
		return othergens, yy, ddi1
	end
end


function C_d_l__minimal(C,d,l, with_length, int_sol, existing_gens = [], pivot = true)
	"""
	* (C,d,l) 	  -->   vector xx which is a minimal generator homologous to the
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
	if int_sol
		println("int sol")
		@variable(m2, x_pos[1:xlen], Int)
		@variable(m2, x_neg[1:xlen], Int)
		@variable(m2, y[1:ylen], Int)
	else
		println("not int sol")
		@variable(m2, x_pos[1:xlen])
		@variable(m2, x_neg[1:xlen])
		@variable(m2, y[1:ylen])
	end
	## we want to minimize sum abs x such that x differs from our generator by a boundary
	if with_length
		@objective(m2, Min, transpose(C.distVec[C.grainVec[2]])*(x_pos[1:xlen] .+ x_neg[1:xlen]) ) #shortest length
	else
		@objective(m2, Min, sum(x_pos[1:xlen].+x_neg[1:xlen]))
	end
	@constraints(m2, begin x_pos[1:xlen] .>= 0; x_neg[1:xlen] .>= 0 end)
	@constraint(m2, x_pos[1:xlen] .- x_neg[1:xlen] .- ddi1 * y .==  Array(gen))

	## Optimize the function and show the nonzero values of our vector x
	status = optimize!(m2)
	x = x_pos .- x_neg
	xx = spzeros(length(x))
	for i in 1:xlen
		if abs(JuMP.value(x[i])) >= .001
			xx[i] = JuMP.value(x[i])
		end
	end
	yy = spzeros(ylen)
	for i in 1:ylen
		if abs(JuMP.value(y[i])) >= .001
			yy[i] = JuMP.value(y[i])
		end
	end
	return xx, yy, ddi1
end
