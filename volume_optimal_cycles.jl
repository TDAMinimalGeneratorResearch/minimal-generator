"""
This script implements the algorithm in Obayashi's paper on finding volume optimal cycles.
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

function find_area_optimal_cycle(C::hObject, d, l,  bdrMat, integer = false)
	"""
	This function finds the area optimal cycle.
	"""
	bi, di = C.barCode[1][l]
	bdrMat, grainstwo, distances, hverts = findbasis(C, C.pointCloud, bi, di+exp(-16))
	return find_optimal_cycle_helper(C, d, l, bdrMat, true, integer)
end

function find_volume_optimal_cycle(C::hObject, d, l,integer = false)
	"""
	This function finds the volume optimal cycle.
	"""
	bi, di = C.barCode[1][l]
	bdrMat, grainstwo, distances, hverts = findbasis(C, C.pointCloud, bi, di+exp(-16))
	return find_optimal_cycle_helper(C, d, l, bdrMat, distances, hverts, false, integer)
end

function find_optimal_cycle_helper(C::hObject, d::Int, l::Int, bdrMat, distances, triVerts, Area = false, integer = false)
	"""
	This function
	"""
	volGens, success, numOfTri, numEdge, len, volVec = volume_optimize(C, d, l, bdrMat, distances, triVerts, true, false, Area, integer)
	firstTimeSuccess = false
	if success
		firstTimeSuccess = true
		return volGens, numOfTri, len, numEdge, volVec, firstTimeSuccess
	else
		println()
		println("running second time!")
		println()
		volGens, success, numOfTri, numEdge, len, volVec = volume_optimize(C, d, l, bdrMat, distances, triVerts, true, true, Area, integer)
		if success
			return volGens, numOfTri, len, numEdge, volVec, firstTimeSuccess
		else
			print("Unable to find volume optimal cycles. ")
		end
	end
end

function volume_optimize(C, d, l, bdrMat, distances, triVerts, secondTime, positive, Area, integer = false)
	"""
	This is the main algorithm of the volume optimize algorithm.
	INPUT:
	"""
	# find all 2-simplicies born at di.
    bi, di = C.barCode[d][l] # find birth and death filtration distance, bi, di
    num_triangles_born_at_di = findall(x-> abs(x-di)<exp(-16), distances)
    num_edges = size(bdrMat, 1)
	num_triangles = size(bdrMat, 2)
    xlen = num_triangles # number of q+1 simplices born IN (bi, di)
	ylen = num_triangles
    zlen = length(num_triangles_born_at_di)

    m2 = Model(Gurobi.Optimizer)
    ## create variable x which is the optimal volume.
    if integer
        @variable(m2, x_pos[1:xlen], Int)
        @variable(m2, x_neg[1:xlen], Int)
    else
        @variable(m2, x_pos[1:xlen])
        @variable(m2, x_neg[1:xlen])
    end

    @variable(m2, y[1:ylen]) # a vector of the coefficients for the 2 simplices born in (bi, di)
    @variable(m2, z_pos[1:zlen]) # variable z is a vector of the coefficients for the 2 simplices born at di.
    @variable(m2, z_neg[1:zlen])

    ### minimize totalarea or not
    @objective(m2, Min, sum(x_pos[1:xlen] .+ x_neg[1:xlen])) # for smallest sum of coeff.
    if !Area
        @objective(m2, Min, sum(x_pos[1:xlen] .+ x_neg[1:xlen])) # for smallest sum of coeff.
    else
	    totalareaVec = Array{Float64}(undef, length(C.grainVec[d + 2]))
	    @time for i in 1:size(triVerts,2)
	        points = C.pointCloud[:,triVerts[:,i]]
	        totalarea = findArea(points)
	        totalareaVec[i] = totalarea
	    end
        @objective(m2, Min, transpose(totalareaVec) * (x_pos[1:xlen] .+ x_neg[1:xlen]))
    end

    # constraint 1: the generator must contain the edge born AT the birth time of the cycle.

	# extract edges born at bi
	birthGrain = findall(x -> x == bi, C.distVec)[1] # find the corresponding grain numbers of bi and di
    deathGrain = findall(x -> x == di, C.distVec)[1]
	firstcol = minimum(findall(x -> x == birthGrain, C.grainVec[2]))
    lastcol  = maximum(findall(x -> x == birthGrain, C.grainVec[2]))
    edges_born_at_bi = spzeros(num_edges, lastcol - firstcol + 1)
    for i in firstcol:lastcol
        edges_born_at_bi[i, i + 1 - firstcol] = 1
    end

    if secondTime
        if positive
            @constraint(m2, (transpose((bdrMat * (x_pos[1:xlen] .- x_neg[1:xlen])))) * edges_born_at_bi .>= 1.0)
        else
            @constraint(m2, (transpose((bdrMat * (x_pos[1:xlen] .- x_neg[1:xlen])))) * edges_born_at_bi .<= -1.0)
        end
    end

    # constraint2: the optimal volume contains the triangles born at di and triangles born in (bi, di)

	triangles_born_at_di = spzeros(xlen, zlen)
	for i in xlen-zlen+1: xlen
		triangles_born_at_di[i, i - ylen+zlen] = 1
	end
	triangles_born_in_bidi = sparse(I, xlen, xlen)

    @constraint(m2,  (x_pos[1:xlen] .- x_neg[1:xlen]) .== triangles_born_at_di * (z_pos[1:zlen] .- z_neg[1:zlen]) .+ triangles_born_in_bidi * y)

    # constraint3: the generator does not contain any edge born in (bi, di)
	firstcol = minimum(findall(x -> x < birthGrain, C.grainVec[d + 1])) # index of the first simplex born after the birth of the cycles
    lastcol  = maximum(findall(x -> x >= deathGrain, C.grainVec[d + 1])) # index of the last simplex born before the death of the cycles
    edges_born_in_bidi = spzeros(num_edges, lastcol - firstcol + 1) # a matrix whose ij th entry is 1 if there is a q simplex born at step i.
    for i in firstcol:lastcol
        edges_born_in_bidi[i, i + 1 -  firstcol] = 1
    end
    @constraint(m2, transpose(bdrMat * (x_pos[1:xlen] .- x_neg[1:xlen])) * edges_born_in_bidi .== 0)

    @constraints(m2, begin x_pos[1:xlen] .>= 0; x_neg[1:xlen] .>= 0 end)
    @constraint(m2, sum(z_pos .+ z_neg) .>= 0.00000001) # coefficients for the triangles born at di cannot be 0.
    @time status = optimize!(m2)
    x = x_pos .- x_neg
    volVec = spzeros(Rational,xlen)
    a = spzeros(xlen)
    z = z_pos .- z_neg

    for i in 1:xlen
        if abs(JuMP.value(x[i])) >= .001
            volVec[i] = JuMP.value(x[i])
            a[i] = 1
        end
    end

    volGens = dropzeros(convert(SparseMatrixCSC{Rational,Int64}, bdrMat * volVec[:,:]))
    success = abs((transpose(bdrMat * volVec[:,:]) * sparse(edges_born_at_bi[:,:]))[1]) >= 0.8

	if Area
    	totalarea = transpose(totalareaVec) * a
	end
    numOfTri = length(findall(x -> x!= 0, volVec))
    numEdge = length(findall(x->x!=0, volGens.nzval))
    a = spzeros(length(volGens),1)
    a[findall(x->x!=0, volGens),1] .= 1
    len = (transpose(C.distVec[C.grainVec[2]][1:size(edges_born_in_bidi,1)]) * a)[1]
	if !Area
		return volGens, success, numOfTri, numEdge, len, volVec
	else
		return volGens, success, numOfTri, numEdge, len, volVec, totalarea
	end
end
