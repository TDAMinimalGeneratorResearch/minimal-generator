include("EireneCode.jl")
function transpose_sparse(M)
	"""
	Takes in an array A, returns the transpose of it in sparse form.
	"""
	return sparse(transpose(M))
end


function flip(M)
	"""
	Takes in an array M, returns this matrix in flipped upside down form.
	(reverses the order of the rows).
	"""
	num_rows = size(M,1)
	flipPerm = num_rows + 1 .- Vector(1:num_rows)
	return M[flipPerm, :]
end

function flip_sparse(rv, cp, vl, m)
	"""
	takes in a sparse matrix and flips it (puts the rows in backwards order)
	"""
	if length(rv) === 0
		return rv, cp, vl
	else
		newrv = m + 1 .- rv
		for i in 1: length(cp)-1
			range = cran(cp, i)
			newrv[range] = reverse(newrv[range])
			vl[range] = reverse(vl[range])
		end
		return newrv, cp, vl
	end
end

function left_right_flip(rv, cp, vl)
	"""
	Takes a sparse matrix and flips it vertically instead of horizontally
	(reverses the order of the columns)
	"""
	reversedrv = reverse(deepcopy(rv))
	reversedvl = reverse(deepcopy(vl))
	newcp = Array{Int64}(undef, length(cp))
	newcp[1] = 1
	for i in 1: length(cp) - 1
		newcp[i+1] = newcp[i] + cp[length(cp) - i + 1] - cp[length(cp) - i]
	end
	for j in 1: length(newcp) - 1
		reversedrv[cran(newcp, j)] = reverse(reversedrv[cran(newcp, j)])
		reversedvl[cran(newcp, j)] = reverse(reversedvl[cran(newcp, j)])
	end
	newcp = newcp[findall(x-> isassigned(newcp,x), 1:length(newcp))]
	return reversedrv, newcp, reversedvl
end

function rv_cp_vl_m__L_R(rv, cp, vl, m, isLast = false)
	"""
	Takes in rv, cp, vl, and m for a boundary matrix and returns L and R
	as well as the reduced matrix (called P). This function assumes the boundary
	matrix's rows and columns are already in filtration order.
	"""
	m = m
	n = length(cp) - 1
	# Mrv_Mcp_Mvl_Mm__reducedform uses the bottom left entry as a pivot and clears
	# to the right. To get it to use the bottom left entry and clear rows above,
	# we reflect our boundary matrix over the off-diagonal line (i.e. flip-transpose-flip it)
	frv, fcp, fvl = flip_sparse(rv, cp, vl, m)
	tfrv, tfcp, tfvl = Mrv_Mcp_Mvl_Mm_rows_cols__slicetranspose(frv, fcp, fvl, m)
	ftfrv, ftfcp, ftfvl = flip_sparse(tfrv, tfcp, tfvl, n)
	L_reduced_dict = Mrv_Mcp_Mvl_Mm__reducedform(ftfrv, ftfcp, ftfvl, n,
												  returntransform=true,
												  returntempreduced=false,
												  returnpivotindices=true)
	# The transformation matrix is a reflected L, so we must flip-trans-flip it back
	Lrv = L_reduced_dict["Trv"]
	Lcp = L_reduced_dict["Tcp"]
	Lvl = L_reduced_dict["Tvl"]

	frv, fcp, fvl = flip_sparse(Lrv, Lcp, Lvl, m)
	tfrv, tfcp, tfvl = Mrv_Mcp_Mvl_Mm_rows_cols__slicetranspose(frv, fcp, fvl, m)
	ftfrv, ftfcp, ftfvl = flip_sparse(tfrv, tfcp, tfvl, m)

	Lrv = ftfrv
	Lcp = ftfcp
	Lvl = ftfvl
	#  We also flip-trans-flip our boundary matrix back to its orignal state
	reduced_transposed_flipped_sparserv = L_reduced_dict["Prv"]
	reduced_transposed_flipped_sparsecp = L_reduced_dict["Pcp"]
	reduced_transposed_flipped_sparsevl = L_reduced_dict["Pvl"]

	frv, fcp, fvl = flip_sparse(reduced_transposed_flipped_sparserv, reduced_transposed_flipped_sparsecp, reduced_transposed_flipped_sparsevl, n)
	tfrv, tfcp, tfvl = Mrv_Mcp_Mvl_Mm_rows_cols__slicetranspose(frv, fcp, fvl, n)
	ftfrv, ftfcp, ftfvl = flip_sparse(tfrv, tfcp, tfvl, m)

	# we have found L, now we find R column-reducing our row-reduced boundary matrix

	if isLast  # if this is the last bdr matrix, we only need to perform column reduce on
			   # the pivot columns
		pcola = sort(n + 1 .- broadcast(abs, L_reduced_dict["prowa"]))
		newrv = Array{Int64}(undef, length(ftfrv))
		newcp = Array{Int64}(undef, length(ftfcp))
		newvl = Array{Rational{Int64}, 1}(undef, length(ftfvl))
		curLen = 1
		newcp[1] = 1
		for i in 1: length(pcola)
			rvVals= cran(ftfcp, pcola[i])
			newrv[curLen: curLen + length(rvVals) - 1] = ftfrv[rvVals]
			newvl[curLen: curLen + length(rvVals) - 1] = ftfvl[rvVals]
			curLen = curLen + length(rvVals)
			if i != length(pcola)
				newcp[pcola[i] + 1 : pcola[i+1]] .= curLen
			else
				newcp[pcola[i] + 1 : length(newcp)] .= curLen
			end
		end
		ftfcp   = newcp
		ftfrv   = newrv
		ftfvl   = newvl
	end
	R_reduced_dict = Mrv_Mcp_Mvl_Mm__reducedform(ftfrv, ftfcp, ftfvl,
													   m,
													   returntransform=true,
													   returntempreduced=false,
													   returnpivotindices=true)

	Rrv = R_reduced_dict["Trv"]
	Rcp = R_reduced_dict["Tcp"]
	Rvl = R_reduced_dict["Tvl"]

	Pcp = R_reduced_dict["Pcp"]
	Pvl = R_reduced_dict["Pvl"]
	Prv = R_reduced_dict["Prv"]

	#save the pivot rows and columns
	pcola = R_reduced_dict["pcola"]
	prowa = R_reduced_dict["prowa"]

	# L = L, R = R, reduced boundary = P

	returndict 					=	Dict(	"Lrv"	=> 	Lrv,
											"Lcp"	=> 	Lcp,
											"Lvl" 	=> 	Lvl,
											"Rrv"   =>  Rrv,
											"Rcp"   =>  Rcp,
											"Rvl"   =>  Rvl,
											"Pcp"   => 	Pcp,
											"Pvl"   =>  Pvl,
											"Prv"   =>  Prv,
											"pcola" =>  pcola,
											"prowa" =>  prowa
											)
	return returndict
end

function bdrMatrices__filtrationorder(bdrMatrices)
	"""
	sorts the rows and columns of bdr matrices so that they are in filtration order
	"""
	D = deepcopy(bdrMatrices)
	for i in 1:length(D["rv"])
		# sort rows first
		rowperm = sortperm(sortperm(-D["grain"][i]))
		newrv = rowperm[D["rv"][i]]
		D["rv"][i] = newrv
		D["lverts"][i][:,rowperm] = D["lverts"][i]
		D["grain"][i] = reverse(D["grain"][i][sortperm(D["grain"][i])])
		# sort columns next, by transposing and sorting the rows
		transDi = transpose_sparse(SparseMatrixCSC(length(D["grain"][i]),
												   length(D["grain"][i+1]),
												   D["cp"][i],
												   D["rv"][i],
												   D["vl"][i]))

		colperm = sortperm(sortperm(-D["grain"][i+1]))
		D["hverts"][i][:,colperm] = D["hverts"][i]
		newrv = colperm[transDi.rowval]
		transDi = SparseMatrixCSC(length(D["grain"][i+1]),
								  length(D["grain"][i]),
								  transDi.colptr,
								  newrv,
								  transDi.nzval)
		Di = transpose_sparse(transDi)
		D["cp"][i] = Di.colptr
		D["vl"][i] = Di.nzval
		D["rv"][i] = Di.rowval
		if i === length(D["rv"]) # if we're on the last matrix, change the grain vector.
			# else wait until the next step and change it with the rows of Di+1
			D["grain"][i+1] = reverse(D["grain"][i+1][sortperm(D["grain"][i+1])])
		end
	end
		return D
end


function bdrMatrices__reducedbdr(bdrMatrices)
	"""
	Takes in a dictionary of boundary matrices in sparse form, and returns a dictionary of
	their reduced form using the Greg H process.
	"""
	bdrMatrices = deepcopy(bdrMatrices)
	numMatrices = length(bdrMatrices["rv"])

	rva = bdrMatrices["rv"]
	cpa = bdrMatrices["cp"]
	vla = bdrMatrices["vl"]
	ma  = bdrMatrices["m"]

	# ReducedDicta = Array{Dict}(undef, numMatrices)
	# ReducedDicta = Array{Dict}(undef, numMatrices)
	ReducedDicta = Dict{String, Vector}()
	ReducedDicta["pcola"] = Array{Array{Int64, 1}}(undef, numMatrices)
	ReducedDicta["prowa"] = Array{Array{Int64, 1}}(undef, numMatrices)
	ReducedDicta["Lrv"] = Array{Array{Int64, 1}}(undef, numMatrices)
	ReducedDicta["Lvl"] = Array{Array{Rational{Int64}, 1}}(undef, numMatrices)
	ReducedDicta["Lcp"] = Array{Array{Int64, 1}}(undef, numMatrices)
	ReducedDicta["Rrv"] = Array{Array{Int64, 1}}(undef, numMatrices)
	ReducedDicta["Rvl"] = Array{Array{Rational{Int64}, 1}}(undef, numMatrices)
	ReducedDicta["Rcp"] = Array{Array{Int64, 1}}(undef, numMatrices)
	ReducedDicta["Prv"] = Array{Array{Int64, 1}}(undef, numMatrices)
	ReducedDicta["Pvl"] = Array{Array{Rational{Int64}, 1}}(undef, numMatrices)
	ReducedDicta["Pcp"] = Array{Array{Int64, 1}}(undef, numMatrices)

	grainvec = bdrMatrices["grain"]

	for i = 1: numMatrices
		m = ma[i]
		n = length(bdrMatrices["cp"][i])-1
		if i > 1 # rather than multiply boundary i + 1 on the left by R inverse,
				# we can just delete the rows of boundary i + 1 which correspond
				# to pivot columns in boundary i
			indices = findall(x -> x in broadcast(abs, ReducedDicta["pcola"][i-1]), rva[i])
			deleteat!(rva[i], indices)
			deleteat!(vla[i], indices)
			original = deepcopy(cpa[i])
			for index in indices
				a = findall(x -> x > index, original)
				cpa[i][a] = cpa[i][a] .- 1
			end
		end
		if i === numMatrices
			dict = rv_cp_vl_m__L_R(rva[i],cpa[i],vla[i], ma[i], true) # compute L and R
		else
			dict = rv_cp_vl_m__L_R(rva[i],cpa[i],vla[i], ma[i], false) # compute L and R
		end
		ReducedDicta["pcola"][i] = dict["pcola"]
		ReducedDicta["prowa"][i] = dict["prowa"]
		ReducedDicta["Lrv"][i]   = dict["Lrv"]
		ReducedDicta["Lvl"][i]   = dict["Lvl"]
		ReducedDicta["Lcp"][i]   = dict["Lcp"]
		ReducedDicta["Rrv"][i]   = dict["Rrv"]
		ReducedDicta["Rvl"][i]	 = dict["Rvl"]
		ReducedDicta["Rcp"][i]	 = dict["Rcp"]
		ReducedDicta["Prv"][i]	 = dict["Prv"]
		ReducedDicta["Pvl"][i]	 = dict["Pvl"]
		ReducedDicta["Pcp"][i]	 = dict["Pcp"]

		if i != 1
			rv, cp, vl = left_right_flip(ReducedDicta["Lrv"][i],ReducedDicta["Lcp"][i], ReducedDicta["Lvl"][i])
			Ln = length(ReducedDicta["Pcp"][i - 1]) - 1
			reduced = rv_cp_vl_m__L_R(rv,cp,vl,Ln)
			Li1idinv = I + SparseMatrixCSC(Ln, Ln, reduced["Lcp"], reduced["Lrv"], reduced["Lvl"] )

			S = SparseMatrixCSC(ma[i-1],
								length(ReducedDicta["Pcp"][i-1]) - 1,
								ReducedDicta["Pcp"][i-1],
								ReducedDicta["Prv"][i-1],
								ReducedDicta["Pvl"][i-1]) * Li1idinv
			ReducedDicta["Prv"][i-1]	 = S.rowval
			ReducedDicta["Pvl"][i-1]	 = S.nzval
			ReducedDicta["Pcp"][i-1]	 = S.colptr
		end

		Ra = I + SparseMatrixCSC(n, n,
											   ReducedDicta["Rcp"][i],
											   ReducedDicta["Rrv"][i],
											   ReducedDicta["Rvl"][i])
		La = I + SparseMatrixCSC(m, m,
										    ReducedDicta["Lcp"][i],
											ReducedDicta["Lrv"][i],
											ReducedDicta["Lvl"][i])

	end
	return ReducedDicta,grainvec,bdrMatrices["lverts"], bdrMatrices["hverts"]
end


# Next: compute H_i by looking at the appropriate columns of R_i L_i+1^-1

# assume boundary i+1 exists

function reduceddicta_i__gens_bars(ReducedDicta, grainvec, i)
	if i === 0
		cols1 = ReducedDicta["Pcp"][1]
		rows = ReducedDicta["Prv"][1]
		generators = Array{SparseMatrixCSC}(undef, length(grainvec[1]))
		persistence = Array{Array{Rational{Int64},1}}(undef, length(grainvec[1]))
		for j in 1:length(grainvec[1])
			bi = maximum(grainvec[1]) + 1 - grainvec[1][j]
			jindex = findall(x->x===j, rows)
			if isempty(jindex)
				di = Inf
			else
				jval = jindex[1]
				jcol = findall(x->x===jval, cols1)[1]
				di = maximum(grainvec[1]) + 1 - grainvec[i + 2][jcol]
			end
			if bi != di
				gensj = Array{Rational{Int64},1}(undef,length(grainvec[1])).=0
				gensj[j] = 1
				generators[j] = gensj
				persistence[j] = [bi, di]
			end
		end
		defindices = filter(i -> isassigned(generators, i), 1:length(generators))
		generators = generators[defindices]
		persistence = persistence[defindices]
		return generators,persistence

	end
	if i >= length(grainvec)-1  # if it is greater max dim return 0 generators
		generators = Vector{Array{Rational{Int64},1}}(undef,0)
		persistence = Vector{Array{Rational{Int64},1}}(undef,0)
		return generators,persistence
	end
	# find homology through LR decomposition
	Rn = length(ReducedDicta["Rcp"][i]) - 1  	# dimension of Ri
	Ln = length(ReducedDicta["Lcp"][i+1]) - 1   # dimension of Li+1

	Ri = SparseMatrixCSC{Rational{Int64}, Int64}(Rn, Rn,
												 ReducedDicta["Rcp"][i],
												 ReducedDicta["Rrv"][i],
												 ReducedDicta["Rvl"][i])  # Ri matrix
	Riid = Ri + I  # add back the identity

	# find Li+1 inverse through row reduction. Need to first flip the matrix left to right.
	rv, cp, vl = left_right_flip(ReducedDicta["Lrv"][i+1],
								 ReducedDicta["Lcp"][i+1],
								 ReducedDicta["Lvl"][i+1])

	reduced = rv_cp_vl_m__L_R(rv,cp,vl,Ln) # inversed Li+1

	Li1idinv = I + SparseMatrixCSC{Rational{Int64}, Int64}(Ln, Ln,
														   reduced["Lcp"],
														   reduced["Lrv"],
														   reduced["Lvl"] ) # add back the identity

	goodmatrix = Riid * Li1idinv  # calculate Ri * (Li+1)^-1

	# find birth and death times for each pivot col. in reduced bdr i
	# find zero columns of bdr which are not nonzero rows of bdr1

	ntimesteps = maximum(grainvec[1])  # num of filtration steps is the maximum

	cols = ReducedDicta["Pcp"][i]
	cols1 = ReducedDicta["Pcp"][i+1]
	Pn = length(cols) - 1 # number of columns in the reduced matrix
	generators = Array{SparseMatrixCSC{Rational}}(undef, Pn)
	persistence = Array{Array{Rational{Int64},1}}(undef, Pn)
	for j in 1 : Pn
		if cols[j] === cols[j+1] # if this is a zero column
			bi = ntimesteps + 1 - grainvec[i + 1][j] # birth time comes from grain vector for boundary i
			# have to find the column with an entry in the jth row in order to find di
			if j in ReducedDicta["prowa"][i + 1] || j in ReducedDicta["prowa"][i + 1] .* -1
				jval = findall(x -> abs(x) === j, ReducedDicta["prowa"][i + 1])[1]
				jcol = abs(ReducedDicta["pcola"][i + 1][jval])
				di   = ntimesteps + 1 - grainvec[i + 2][jcol]
			else
				di = Inf
			end
			if bi != di && di != Inf
				col = zeros(Rational{Int64}, length(grainvec[i+1]))
				col[ReducedDicta["Prv"][i + 1][cran(ReducedDicta["Pcp"][i+1], jcol)]] = ReducedDicta["Pvl"][i + 1][cran(ReducedDicta["Pcp"][i + 1], jcol)]
				generators[j] = goodmatrix * sparse(col)
				persistence[j] = [bi, di]
			end
		end
	end

	defindices = filter(i -> isassigned(generators, i), 1:length(generators))
	generators = generators[defindices]
	persistence = persistence[defindices]
	return generators, persistence, goodmatrix
end


##############

function pointcloud_maxDim__bdrmatrices(pointcloud, maxdim = 1)
	"""
	takes in a pointcloud, returns bdr matrices dictionary, by default max dim = 1.
	"""
	C = eirene(pointcloud, maxdim = maxdim, model = "pc") # calls Eirene on the pointcloud and gets a dictionary C that contains useful information

	# retrieves information from C
	bdrMatrices = Dict{String, Vector}()  # storing bdr matrices info in a dictionary
	bdrMatrices["rv"] 	 	= 	Array{Array{Int64, 1}}(undef, maxdim + 1)
	bdrMatrices["cp"] 	 	= 	Array{Array{Int64, 1}}(undef, maxdim + 1)
	bdrMatrices["vl"] 	 	= 	Array{Array{Rational, 1}}(undef, maxdim + 1) # the value entries should be rational.
	bdrMatrices["m"] 	 	= 	Array{Int64}(undef, maxdim + 1) # m is the number of rows of the bdr matrix.
	bdrMatrices["grain"] 	= 	C["grain"][1:maxdim + 2] # grain vectors gives us information about the birth time of each simplex.
	bdrMatrices["lverts"] 	= 	Array{Array{Int64, 2}}(undef, maxdim + 1) # row labels
	bdrMatrices["hverts"] 	= 	Array{Array{Int64, 2}}(undef, maxdim + 1) # col labels

	# find bdr matrices from dimension 1 up to the max dimension specified by the users and store them in a dictionary.
	for i in 1: maxdim + 1
		brv,bcp,bvl,bm,bn,lowverts,higverts = eirened_dim__signedboundarymatrix(C, dim = i; useoriginalvertexlabels=true)
		bdrMatrices["rv"][i]  		= 	brv
		bdrMatrices["cp"][i]  		= 	bcp
		bdrMatrices["vl"][i]  		= 	bvl
		bdrMatrices["m"][i]   		= 	bm
		bdrMatrices["lverts"][i] 	= 	lowverts
		bdrMatrices["hverts"][i] 	= 	higverts
	end
	return bdrMatrices
end


struct hObject
		"""
		a custom object to store all the information we need for calculations below
		"""
		pointCloud :: Array{Float64, 2}
		bdrMatrices :: Dict{String, Vector}
		reducedbdr :: Dict{String, Vector}
		barCode :: Vector{Array{Array{Float64, 1}, 1}}
		grainVec :: Vector{Array{Int64, 1}}
		distVec :: Vector{Float64}
		pcola :: Vector{Array{Int64, 1}}
		prowa :: Vector{Array{Int64, 1}}
		generators :: Vector{Array{SparseMatrixCSC, 1}}
		permutedlverts :: Vector{Array{Int64, 2}}
		permutedhverts :: Vector{Array{Int64, 2}}
end

# function pc_lverts__distvec(pc,permutedlverts)
# 	"""
# 	creates the vector distvec (dv) which contains grain information in terms of
# 	actual distances (e.g. dv[33] represents the epsilon value of an edge with
# 	grain = 33)
# 	"""
# 	dv = zeros(size(permutedlverts[2])[2] + 1)
# 	dv[size(permutedlverts[2])[2] + 1] = 0
# 	for i in 1:size(permutedlverts[2])[2]
# 		dv[i] = euclidean(pc[:,permutedlverts[2][1,i]],pc[:,permutedlverts[2][2,i]])
# 	end
# 	dv = unique(dv)
# 	dv = dv[sortperm(-dv)]
# 	return dv
# end

function pc_lverts__distvec(pc,grainvec,permutedlverts)
	"""
	creates the vector distvec (dv) which contains grain information in terms of
	actual distances (e.g. dv[33] represents the epsilon value of an edge with
	grain = 33)
	"""
	dv = zeros(maximum(grainvec[2]) + 1)
	dv[maximum(grainvec[2]) + 1] = 0
	for i in 1:length(grainvec[2])
		if i == length(grainvec[2])
			dv[grainvec[2][i]] = euclidean(pc[:,permutedlverts[2][1,i]],pc[:,permutedlverts[2][2,i]])
		elseif grainvec[2][i] != grainvec[2][i+1]
			dv[grainvec[2][i]] = euclidean(pc[:,permutedlverts[2][1,i]],pc[:,permutedlverts[2][2,i]])
		end
	end
	return dv
end 

function computeHomology(pointcloud, maxDim = 1)
	"""
	performs the full gregH "pipeline" and creates a homologyObject (as above)
	which contains information about boundary matrices, generators, etc.
	"""
	bdrMatrices 	 = pointcloud_maxDim__bdrmatrices(pointcloud, maxDim)
	filtrationorder  = bdrMatrices__filtrationorder(bdrMatrices)
	@time reduced,grainvec,permutedlverts, permutedhverts = bdrMatrices__reducedbdr(filtrationorder)
	generators = Array{Array{SparseMatrixCSC}}(undef, maxDim)
	barCodes = Array{Array{Array{Rational{Int64},1}}}(undef, maxDim)
	for i in 1: maxDim
		@time gens,bars= reduceddicta_i__gens_bars(reduced,grainvec,i)
		generators[i] = gens
		barCodes[i] = bars
	end
	distVec = pc_lverts__distvec(pointcloud, grainvec, permutedlverts)
	# convert barcode to distances
	ldv1 = length(distVec) + 1
	distcode = Vector{Vector{Array{Float64,1}}}(undef, maxDim)
	for d in 1:maxDim
		distcode[d] = Vector{Array{Float64,1}}(undef,length(barCodes[d]))
		for i in 1:length(barCodes[d])
			distcode[d][i] = [distVec[ldv1 - Int(barCodes[d][i][1])],distVec[ldv1 - Int(barCodes[d][i][2])]]
		end
	end

	C = hObject(
		pointcloud,
		filtrationorder,
		reduced,
		distcode,
		grainvec,
		distVec,
		reduced["pcola"],
		reduced["prowa"],
		generators,
		permutedlverts,
		permutedhverts)
	return C
end
