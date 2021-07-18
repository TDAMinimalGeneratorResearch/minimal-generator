########## volume optimal cycles
"""
This method implements the algorithm in Obayashi's paper on volume optimal cycles.
The pseudo algorithm can be found on page 520 of the original paper.
"""

function findAreaOptimalCycle(C, d, l, dd, triVerts, integer = false)
	return findVolOptimalCycle(C, d, l, dd, triVerts, true, integer)
end

function findVolumeOptimalCycle(C, d, l, dd, triVerts, integer = false)
	return findVolOptimalCycle(C, d, l,  dd, triVerts,false, integer)
end

function findVolOptimalCycle(C, d, l, dd, triVerts, Area = false, integer = false)
	bi, di = C.barCode[1][l]
	volGens, success, numOfTri, numEdge, len, area, xx ,zz = volumeOptimize(C, d, l,dd, triVerts, false, true, Area, integer)
	firstTimeSuccess = false
	if success
		firstTimeSuccess = true
		return volGens, numOfTri, len, numEdge,area, xx,zz,  firstTimeSuccess
	else
		println()
		println("running second time!")
		println()
		volGens, success, numOfTri, numEdge, len, area, xx ,zz = volumeOptimize(C, d, l, dd, triVerts, true, false, Area, integer)
		volGens1, success1, numOfTri1, numEdge1, len1, area1, xx1 ,zz1 = volumeOptimize(C, d, l, dd, triVerts,true, true, Area, integer)
		if area < area1
			return volGens, numOfTri, len, numEdge, area, xx,zz,  firstTimeSuccess
		else
			return volGens1, numOfTri1, len1, numEdge1, area1, xx1,zz1, firstTimeSuccess
			print("Unable to find volume optimal cycles. ")
		end
	end
end

function constructInput(D, d, l)
		bi, di = D.barCode[d][l] # find birth and death filtration distance, bi, di
		ddi1 = findbasis(D, D.pointCloud, bi, di)
	    deathGrain = findall(x -> x == di, D.distVec)[1]
		dis = findall(x->x == deathGrain, D.grainVec[3])
		# pivotdis = dis[findall(in(broadcast(abs,D.pcola[2])), dis)] # triangles born at the death time that are pivot columns

		sigmadi = Int(D.simplexcode[d][l][2]) # the trianle born at the death time of the cycle
		divec = pivotdis[findall(x->x<= sigmadi, pivotdis)] # all pivot triangles born before the death triangle
		dd = ddi1[1]

		didi = ddi1[1][:, findall(x->x==di, ddi1[3])]
		dd = dd[:, 1:(size(dd,2) - size(didi,2))] # all triangles born before the death time of the cycle
		oribdr = SparseMatrixCSC(length(D.bdrMatrices["cp"][1])-1, length(D.bdrMatrices["cp"][2])-1, D.bdrMatrices["cp"][2], D.bdrMatrices["rv"][2],  D.bdrMatrices["vl"][2])
		dd = hcat(dd, oribdr[1:size(dd,1), dis]) # all triangles born before the death time, and triangles born at the death time are pivots
		triVerts = hcat(ddi1[4][2][:, 1:(size(ddi1[1],2) - size(didi,2))], D.permutedhverts[2][:,dis])
		return dd, triVerts
end

function constructInputAll(D, d, l, a)
		bi, di = D.barCode[d][l] # find birth and death filtration distance, bi, di
	    deathGrain = findall(x -> x == di, D.distVec)[1]
		dis = findall(x->x == deathGrain, D.grainVec[3])
		pivotdis = dis[findall(in(broadcast(abs,D.pcola[2])), dis)] # triangles born at the death time that are pivot columns
		sigmadi = Int(D.simplexcode[d][l][2]) # the trianle born at the death time of the cycle
		divec = pivotdis[findall(x->x<= sigmadi, pivotdis)] # all pivot triangles born before the death triangle
		x = maximum(findall(x->  x==di, a[5]))
		keep = findall(x-> bi<= x<di, a[3])
		dd = a[1][ 1: x,keep]
		oribdr = SparseMatrixCSC(length(D.bdrMatrices["cp"][1])-1, length(D.bdrMatrices["cp"][2])-1, D.bdrMatrices["cp"][2], D.bdrMatrices["rv"][2],  D.bdrMatrices["vl"][2])
		dd = hcat(dd, oribdr[1:size(dd,1), divec]) # all triangles born before the death time, and triangles born at the death time are pivots
		triVerts = hcat(a[4][2][:, keep], D.permutedhverts[2][:,divec])
		return dd, triVerts
end


function volumeOptimize(D, d, l, dd, triVerts, secondTime, positive, Area, integer = false)

	numEdges = size(dd, 1)
    xlen = size(dd,2)
    zlen = 1
	ylen = xlen  -zlen
    m2 = Model(GLPK.Optimizer) # Clp.Optimizer
    ## create variable x which is the optimal volumD.
	println("variable")
	println(xlen)

    if integer
        @variable(m2, x_pos[1:xlen], Int)
        @variable(m2, x_neg[1:xlen], Int)
    else
        @variable(m2, x_pos[1:xlen] )
        @variable(m2, x_neg[1:xlen])
    end

	println(ylen)
    # @variable(m2, y[1:ylen]) # a vector of the coefficients for the q+1 simplices born in (bi, di)
    # @variable(m2, z[1:zlen]) # variable z is a vector of the coefficients for the q+1 simplices born at di.
	println("done")
	# totalareaVec = Array{Float64}(undef, size(triVerts,2))
	# # dm = pairwise(Euclidean(), pc)
	# for i in 1:size(triVerts,2)
	# 	points = triVerts[:,i]
	# 	totalarea = findAreaDM(points, D.pointCloud)
	# 	global totalareaVec[i] = totalarea
	# end
	if !Area
		@objective(m2, Min, sum(x_pos[1:xlen] .+ x_neg[1:xlen])) # for smallest sum of coeff.
	else
		@objective(m2, Min, transpose(totalareaVec) * (x_pos[1:xlen] .+ x_neg[1:xlen]))
	end
	# using simplexwise filtration
	firstcol = Int(D.simplexcode[d][l][1])

    # sigmabi = spzeros(numEdges,  1)
	# sigmabi[firstcol,1] = 1
    # constraint1: the generator must contain the edge born at the birth time of the cycle.
    if secondTime
		if positive
            @constraint(m2, (transpose((dd[:,1:xlen] * (x_pos[1:xlen] .- x_neg[1:xlen]))))[firstcol] .>= 1.0)
        else
            @constraint(m2, (transpose((dd[:,1:xlen] * (x_pos[1:xlen] .- x_neg[1:xlen]))))[firstcol] .<= -1.0)
        end
    end

    # constraint3: the generator does not contain any edge born in (bi, di)
    @constraint(m2, transpose(dd[firstcol+1:numEdges,1:xlen] * (x_pos[1:xlen] .- x_neg[1:xlen])) .== 0)
    # @constraints(m2, begin x_pos[1:xlen] .<= 1; x_neg[1:xlen] .<= 1 end)
	@constraints(m2, begin x_pos[1:xlen] .>= 0; x_neg[1:xlen] .>= 0 end)
	# constraint2: the optimal volume contains the triangles born at di and triangles born in the range (bi, di)
    @constraint(m2, (x_pos[1:xlen] .- x_neg[1:xlen])[xlen] .== 1)
    @time status = optimize!(m2)
    x = x_pos .- x_neg
    xx = spzeros(Float64,xlen)
    a = spzeros(xlen)
    for i in 1:xlen
        if abs(JuMP.value(x[i])) >= .001
            xx[i] = JuMP.value(x[i])
            a[i] = 1
        end
    end
	zz = spzeros(Float64,zlen)
	# for i in 1:zlen
    #     if abs(JuMP.value(z[i])) >= .00000001
    #         zz[i] = JuMP.value(z[i])
    #     end
    # end
    volGens = dropzeros(convert(SparseMatrixCSC{Float64,Int64}, dd[:,1:xlen] * xx[:,:]))


    success = (firstcol in (dropzeros(dd[:,1:xlen] * xx[:,:]).rowval))
    totalarea = -1# transpose(totalareaVec) * a
    numOfTri = length(findall(x -> x!= 0, xx))
    numEdge = length(findall(x->x!=0, volGens.nzval))
    a = spzeros(length(volGens),1)
    a[findall(x->x!=0, volGens),1] .= 1
    len = (transpose(D.distVec[D.grainVec[2]][1:numEdges]) * a)[1]
	return volGens, success, numOfTri, numEdge, len, totalarea, xx, zz
end
