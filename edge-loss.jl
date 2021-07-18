"""
Compute uniform/length weighted edge-loss optimal cycles.
prs: persistent homological cycle basis
fcb: filtered complex basis
"""


function findUnifEdgeOptimalCycles_prs(C, d, l, optimized_gens, intSol = false, verbose = false, pivot = true)
	bi,di = C.barCode[1][l]
	bcVec = hcat(C.barCode[1]...)
	biVec = bcVec[1,:] # all the birth distances as a vector
	diVec = bcVec[2,:] # all the death distances as a vector
	genIdx = intersect(findall(x -> x <= bi, biVec), findall(x-> x <= di, diVec)) # all the indices with a death time smaller than di
	for i in filter(i -> !isassigned(optimized_gens, i), 1:l)
		if (i in genIdx)
			deleteat!(genIdx, findall(x->x==i, genIdx)[1])
		end
	end
	prev = vcat(optimized_gens[filter(i -> isassigned(optimized_gens, i), 1:l-1)], C.generators[1][l:length(C.generators[1])])
	result = C_d_l__minimal(C, 1, l, false, intSol, prev[genIdx], pivot)
	if length(result) < 4
		return false
	else
		othergens, yy, ddi1, len = result
		nonzero = findall(x->x!=0, othergens)
		if verbose
			return othergens, len, length(nonzero), yy, ddi1
		else
			return othergens, yy, ddi1
		end
	end
end

function findUnifEdgeOptimalCycles_fcb(C, d, l, intSol = false, verbose = false, pivot = true)
	result = C_d_l__minimal(C, 1, l, false, intSol, C.generators[1][1:l-1], pivot)
	if length(result) < 4
		return false
	else
		othergens, yy, ddi1, len = result
		nonzero = findall(x->x!=0, othergens)
		if verbose
			return othergens, len, length(nonzero), yy, ddi1
		else
			return othergens, yy, ddi1
		end
	end
end

function findLenEdgeOptimalCycles_fcb(C, d, l, intSol = false, verbose = false, pivot = true)
	result = C_d_l__minimal(C, 1, l, true, intSol, C.generators[1][1:l-1], pivot)
	if length(result) < 4
		return false
	else
		othergens, yy, ddi1, len = result
		nonzero = findall(x->x!=0, othergens)
		if verbose
			return othergens, len, length(nonzero), yy, ddi1
		else
			return othergens, yy, ddi1
		end
	end
end

function findLenEdgeOptimalCycles_prs(C, d, l,optimized_gens, intSol = false, verbose = false, pivot = true)
	bi,di = C.barCode[1][l]
	bcVec = hcat(C.barCode[1]...)
	biVec = bcVec[1,:] # all the birth distances as a vector
	diVec = bcVec[2,:] # all the death distances as a vector
	genIdx = intersect(findall(x -> x <= bi, biVec), findall(x-> x <= di, diVec)) # all the indices with a death time smaller than di
	for i in filter(i -> !isassigned(optimized_gens, i), 1:l)
		if (i in genIdx)
			deleteat!(genIdx, findall(x->x==i, genIdx)[1])
		end
	end
	prev = vcat(optimized_gens[filter(i -> isassigned(optimized_gens, i), 1:l-1)], C.generators[1][l:length(C.generators[1])])
	result = C_d_l__minimal(C, 1, l, true, intSol, prev[genIdx], pivot)
	if length(result) < 4
		return false
	else
		othergens, yy, ddi1, len = result
		nonzero = findall(x->x!=0, othergens)
		if verbose
			return othergens, len, length(nonzero), yy, ddi1
		else
			return othergens, yy, ddi1
		end
	end
end


function C_d_l__minimal(C,d,l, with_length_or_custom, int_sol, existing_gens = [], pivot = true)
	"""
	* (C,d,l) 	  --> vector xx which is a minimal generator homologous to the
						lth d-dimensional generator (c_l) of C, unless existing gens
						are given, in which case it is homologous to a sum of
						c_l and existing_gens.
						with_length_or_custom indicates whether to
						minimize euclidean distance (if false, will instead
						minimize number of simplices) *don't set with_length_or_custom to
						true for dimension not equal to 1
						setting int_sol = false may be faster, but may produce
						bad solutions sometimes
	"""
	bdrMatrices = C.bdrMatrices
	ddi1 = C_d_l__filt(C,d,l, pivot)
	gen = C_d_l__gen(C,d,l)
	bi, di = C.barCode[1][l]
	birthGrain = findall(x->x==bi, C.distVec)[1]
	lastcol  = maximum(findall(x -> x == birthGrain, C.grainVec[2]))
	ddi1 = ddi1[collect(1:lastcol), :]
    #if given existing generators, include them as cols of bdr matrix
	if length(existing_gens) != 0
		for i in 1 : length(existing_gens)
			ddi1 = hcat(ddi1, existing_gens[i][1:lastcol])
		end
	end
	## establish model and varible lengths
	xlen = ddi1.m
	ylen = length(ddi1.colptr) - 1
	m2 = Model(GLPK.Optimizer) # Clp.Optimizer
	## create variable x which will be our new generator, and y for the boundary i+1
	if int_sol  # --- later, int.  vs. {0, \pm 1}
		println("int sol")
		@variable(m2, x_pos[1:xlen], Int)
		@variable(m2, x_neg[1:xlen], Int)
		# @variable(m2, y[1:ylen], Int)
	else
		println("not int sol")
		@variable(m2, x_pos[1:xlen])
		@variable(m2, x_neg[1:xlen])

	end
	@variable(m2, y[1:ylen])
	## we want to minimize sum abs x such that x differs from our generator by a boundary
	if typeof(with_length_or_custom) == Bool
		if with_length_or_custom
			@objective(m2, Min, transpose(C.distVec[C.grainVec[2]][1:xlen]) * (x_pos[1:xlen] .+ x_neg[1:xlen]) ) #shortest length
		else
			@objective(m2, Min, sum(x_pos[1:xlen].+x_neg[1:xlen]))
		end
	else
		@objective(m2, Min, with_length_or_custom * (x_pos[1:xlen] .+ x_neg[1:xlen]) ) #shortest length
	end
	@constraints(m2, begin x_pos[1:xlen] .>= 0; x_neg[1:xlen] .>= 0 end)
	@constraint(m2, x_pos[1:xlen] .- x_neg[1:xlen] .- ddi1 * y .==  convert(SparseMatrixCSC{Rational{Int64}, Int64}, gen[1:lastcol]))
	# Optimize the function and show the nonzero values of our vector x
	status = optimize!(m2)
	x = x_pos .- x_neg
	xx = spzeros(xlen,1)
	if has_values(m2)
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
		len = sum(C.distVec[C.grainVec[2][1:lastcol]][xx.rowval])
		return xx, yy, ddi1, len
	else
		return false
	end
end
