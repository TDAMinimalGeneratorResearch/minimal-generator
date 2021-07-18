using LinearAlgebra
using Distances
using Distributions

"""

"""

function distance_matrix__bdrMatrix(distance_matrix, verts, lower, upper)
	"""
	Given a distance matrix and distance threshold, change all distances outside
	the input range to Inf
	"""
    t = 1
    m = size(distance_matrix)[1]
    n = size(distance_matrix)[2]
    distance_matrix_new = deepcopy(distance_matrix)
    for i in 1:m
        for j in i+1:m
            if  distance_matrix[i, j] <= upper
                distance_matrix_new[i, j] = t
                distance_matrix_new[j, i] = t
                t = t + 1
            else
                distance_matrix_new[i, j] = Inf
				distance_matrix_new[j, i] = Inf
            end
        end
    end
    return distance_matrix_new
end

function build_lowfaces(distance_matrix, lower, upper, verts)
	n = length(verts)
    lowfaces = Array{Int}(undef, 2, Int(n * (n-1)/2))
	mapped = Array{Int}(undef, 2, Int(n * (n-1)/2))
    t = 1
    distances = Array{Float64}(undef,0)
    for i in 1 : n
        for j in i + 1 : n
            if  distance_matrix[i,j] <= upper
                lowfaces[t : t + 1,] = [i,j]
				mapped[t : t + 1,] = sort([verts[i], verts[j]])
                append!(distances, distance_matrix[i, j])
                t = t + 2
            end
        end
    end
    return lowfaces[1:2, 1:Int((t+1)/2) - 1], mapped[1:2, 1:Int((t+1)/2) - 1], distances, length(distances) .- sortperm(sortperm(distances)) .+ 1
end


function build_highfaces(lower, upper, distance_matrix, edgesLength, verts)
    t = 1
    distances = Array{Float64}(undef,0)
	n = length(verts)
    highfaces = Array{Int}(undef, 3, Int(n * (n-1) * (n-2)/6))
	mapped = Array{Int}(undef, 3, Int(n * (n-1) * (n-2)/6))
	grain =  Array{Int}(undef, 0)
	previousGrain = sortperm(sortperm(edgesLength))
    for i in 1: n
        for j in i + 1: n
            for k in j + 1: n
				if !isinf(distance_matrix[i,j]) && distance_matrix[i,j]> 0 && !isinf(distance_matrix[i,k]) && distance_matrix[i,k]> 0 && !isinf(distance_matrix[j,k]) && distance_matrix[j,k]> 0
					ij = Int(distance_matrix[i, j])
					jk = Int(distance_matrix[j, k])
					ik = Int(distance_matrix[i, k])
	                max_length = maximum([edgesLength[ij], edgesLength[ik], edgesLength[jk]])
					curgrain = maximum([previousGrain[ij], previousGrain[ik], previousGrain[jk]])
	                if lower <= max_length <= upper
	                    highfaces[t : t + 2,] = [Int(ij), Int(ik), Int(jk)]
						mapped[t : t + 2,] = sort([verts[i],verts[k],verts[j]])
	                    append!(distances, max_length)
						append!(grain, curgrain)
	                    t = t + 3
	                end
				end
            end
        end
    end
    return highfaces[1:3, 1:Int((t+2)/3) - 1], mapped[1:3, 1:Int((t+2)/3) - 1], distances, length(edgesLength) .- grain .+ 1
end

function dim__signedboundarymatrix(ddim, simplices, lowverts, higverts)
	brv = vec(reshape(simplices, (:,1)))
	bcp = [1:size(simplices)[1]:length(brv)+1;]
	bvl 					=	zeros(Rational{Int64},length(brv))
	# 	COUNT THE NUMBER OF CELLS IN DIMENSION dim AND dim-1
	nhig 					=	size(higverts)[2]
	nlow 					=	size(lowverts)[2]


	# 	OBTAIN SIGNS FOR THE NONZERO ENTRIES
	bvlind 					=	0
	for celhig = 1:nhig  								# celhig := a cell of the higher dimension
		vtahig 				=	higverts[:,celhig] 		# vertices of the higher cell
		cellowa 			=	crows(bcp,brv,celhig) 	# nonzero rows of the column indexed by celhig
		for cellow = cellowa 							# cellow := a cell of the lower dimension (and a face of the higher cell)
			vtalow 			=	lowverts[:,cellow]		# vertices of the lower cell
			p 				=	0
			for slot = 1:ddim
				if vtalow[slot] != vtahig[slot]
					p 		=	slot
					break
				end
			end
			if p == 0
				p = ddim+1 # recall that the high dim cell has dim+1 vertices
			end
			if p % 2 == 0
				coefficient 	=	-1 		# if the missing vertex is in an even location, give a negative entry
			else
				coefficient 	=	1		# otherwise, give a positive
			end

			bvlind 			=	bvlind+1
			bvl[bvlind] 	=	coefficient

		end
	end

	# 	RECOVER THE SIZE OF THE MATRIX
	bm 						=	nlow
	bn 						=	nhig

	# 	RETURN
	return brv,bcp,bvl,bm,bn,lowverts,higverts
end

function signbdrMatrix(verts, lowverts_orig, higverts_orig, maxdim, triangles, edge, grainsone, grainstwo)
	# retrieves information from C
	bdrMatrices = Dict{String, Vector}()  # storing bdr matrices info in a dictionary
	bdrMatrices["rv"] 	 	= 	Array{Array{Int64, 1}}(undef, maxdim + 1)
	bdrMatrices["cp"] 	 	= 	Array{Array{Int64, 1}}(undef, maxdim + 1)
	bdrMatrices["vl"] 	 	= 	Array{Array{Rational, 1}}(undef, maxdim + 1) # the value entries should be rational.
	bdrMatrices["m"] 	 	= 	Array{Int64}(undef, maxdim + 1) # m is the number of rows of the bdr matrix.
	bdrMatrices["grain"] 	= 	Vector{Array{Int, 1}}(undef, 3) # grain vectors gives us information about the birth time of each simplex.
	bdrMatrices["lverts"] 	= 	Array{Array{Int64, 2}}(undef, maxdim + 1) # row labels
	bdrMatrices["hverts"] 	= 	Array{Array{Int64, 2}}(undef, maxdim + 1) # col labels
	bdrMatrices["grain"][1] =   repeat([maximum(grainsone) + 1], length(verts))
	# find bdr matrices from dimension 1 up to the max dimension specified by the users and store them in a dictionary.

	for i in 1: maxdim + 1
		if i == 2
			bdrMatrices["grain"][i+1]   	= 	grainstwo
			brv,bcp,bvl,bm,bn,lowverts,higverts = dim__signedboundarymatrix(2, triangles, lowverts_orig, higverts_orig)
		else
			bdrMatrices["grain"][i+1]   	= 	grainsone
			brv,bcp,bvl,bm,bn,lowverts,higverts = dim__signedboundarymatrix(1, edge, verts, lowverts_orig)
		end
		bdrMatrices["rv"][i]  		= 	brv
		bdrMatrices["cp"][i]  		= 	bcp
		bdrMatrices["vl"][i]  		= 	bvl
		bdrMatrices["m"][i]   		= 	bm

		bdrMatrices["lverts"][i] 	= 	deepcopy(lowverts)
		bdrMatrices["hverts"][i] 	= 	deepcopy(higverts)
	end
	return bdrMatrices
end

function findbasis(D, pc, lower, upper, distmat=false)
	vecverts = vec(D.bdrMatrices["lverts"][1])
	verts = transpose(vecverts)[:,:]
	pcreorder= pc[:,vecverts]
	if distmat
		distance_matrix= pcreorder[vecverts,:]
	else
		distance_matrix = Distances.pairwise(Euclidean(), pcreorder)
	end
	di=upper
	x = distance_matrix__bdrMatrix(distance_matrix, vecverts, 0, di)
	edge, lowverts, edgeLength, grainsone = build_lowfaces(distance_matrix, 0, di, vecverts)
	triangles, higverts, distances,grainstwo = build_highfaces(lower, upper, x, edgeLength, vecverts)
	a = signbdrMatrix(verts, lowverts, higverts, 1, triangles, edge, grainsone, grainstwo)
	filtrationorder = bdrMatrices__filtrationorder(a,1)
	bdr = SparseMatrixCSC(length(filtrationorder["cp"][1]) - 1, length(filtrationorder["cp"][2]) - 1, filtrationorder["cp"][2], filtrationorder["rv"][2], filtrationorder["vl"][2])
	return bdr, filtrationorder["hverts"] 
end
