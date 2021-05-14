include("EireneCode.jl")
using Distances
function transpose_sparse(M)
	"""
	Takes in an array A, returns the transpose of it in sparse form.
	"""
	return sparse(transpose(M))
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

function rv_cp_vl_m__L_R(rv, cp, vl, m, isLast = false)
	"""
	Takes in rv, cp, vl, and m for a boundary matrix and returns L and R
	as well as the reduced matrix (called P). This function assumes the boundary
	matrix's rows and columns are already in filtration order.
	"""
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
	Prv = L_reduced_dict["Prv"]
	Pcp = L_reduced_dict["Pcp"]
	Pvl = L_reduced_dict["Pvl"]

	frv, fcp, fvl = flip_sparse(Prv, Pcp, Pvl, n)
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

function bdrMatrices__filtrationorder(bdrMatrices, maxDim)
	"""
	sorts the rows and columns of bdr matrices so that they are in filtration order
	"""
	for i in 1: (maxDim+1)
		# sort rows first
		rowperm = sortperm(sortperm(-bdrMatrices["grain"][i]))
		newrv = rowperm[bdrMatrices["rv"][i]]
		bdrMatrices["rv"][i] = newrv
		bdrMatrices["lverts"][i][:,rowperm] = bdrMatrices["lverts"][i]
		bdrMatrices["grain"][i] = bdrMatrices["grain"][i][sortperm(-bdrMatrices["grain"][i])]
		# sort columns next, by transposing and sorting the rows
		transDi = transpose_sparse(SparseMatrixCSC(length(bdrMatrices["grain"][i]),
												   length(bdrMatrices["grain"][i+1]),
												   bdrMatrices["cp"][i],
												   bdrMatrices["rv"][i],
												   bdrMatrices["vl"][i]))

		colperm = sortperm(sortperm(-bdrMatrices["grain"][i+1]))
		bdrMatrices["hverts"][i][:,colperm] = bdrMatrices["hverts"][i]
		newrv = colperm[transDi.rowval]
		transDi = SparseMatrixCSC(length(bdrMatrices["grain"][i+1]),
								  length(bdrMatrices["grain"][i]),
								  transDi.colptr,
								  newrv,
								  transDi.nzval)
		Di = transpose_sparse(transDi)
		bdrMatrices["cp"][i] = Di.colptr
		bdrMatrices["vl"][i] = Di.nzval
		bdrMatrices["rv"][i] = Di.rowval
		if i === length(bdrMatrices["rv"]) # if we're on the last matrix, change the grain vector.
			# else wait until the next step and change it with the rows of Di+1
			bdrMatrices["grain"][i+1] = reverse(bdrMatrices["grain"][i+1][sortperm(bdrMatrices["grain"][i+1])])
		end
	end
		return bdrMatrices
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
			indices = findall(x -> x in broadcast(abs, ReducedDicta["pcola"][i-1]), rva[i]) #
			vla[i][indices] .= 0
		    newOne = dropzeros!(SparseMatrixCSC(m, n, cpa[i], rva[i], vla[i]))
			rva[i] = newOne.rowval
			cpa[i] = newOne.colptr
			vla[i] = newOne.nzval
			ma[i] = newOne.m
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
			Ln = length(ReducedDicta["Pcp"][i - 1]) - 1
			Li1idinv = rv_cp_vl_m__inverse( deepcopy(ReducedDicta["Lrv"][i]),
											deepcopy(ReducedDicta["Lcp"][i]),
											deepcopy(ReducedDicta["Lvl"][i]), Ln)
			S = dropzeros(SparseMatrixCSC(ma[i-1],
								Ln,
								ReducedDicta["Pcp"][i-1],
								ReducedDicta["Prv"][i-1],
								ReducedDicta["Pvl"][i-1]) * Li1idinv)
			ReducedDicta["Prv"][i-1]	 = S.rowval
			ReducedDicta["Pvl"][i-1]	 = S.nzval
			ReducedDicta["Pcp"][i-1]	 = S.colptr
		end
	end
	return ReducedDicta,grainvec,bdrMatrices["lverts"], bdrMatrices["hverts"]
end

function rv_cp_vl_m__inverse(rv,cp,vl,m)
	mat = dropzeros(SparseMatrixCSC(m,m,cp,rv,vl) + I)
	flipped = mat[m:-1:1, :]
	reduced = Mrv_Mcp_Mvl_Mm__reducedform(flipped.rowval,flipped.colptr,flipped.nzval, m,returntransform=true,returntempreduced=false,returnpivotindices=true)
	inverse = I + dropzeros(SparseMatrixCSC(m, m, reduced["Tcp"], reduced["Trv"], reduced["Tvl"]))
	return inverse
end

function reduceddicta_i__gens_bars(ReducedDicta, grainvec, i)
	if i === 0
		cols = ReducedDicta["Pcp"][1]
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
				jcol = findall(x->x===jval, cols)[1]
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
	Ri = dropzeros(SparseMatrixCSC{Rational{Int64}, Int64}(Rn, Rn,
												 ReducedDicta["Rcp"][i],
												 ReducedDicta["Rrv"][i],
												 ReducedDicta["Rvl"][i]))  # Ri matrix
	Riid = Ri + I  # add back the identity
	# find Li+1 inverse through row reduction.
	Li1idinv = rv_cp_vl_m__inverse( deepcopy(ReducedDicta["Lrv"][i+1]),
									deepcopy(ReducedDicta["Lcp"][i+1]),
									deepcopy(ReducedDicta["Lvl"][i+1]), Ln)
	goodmatrix = dropzeros(Riid * Li1idinv)  # calculate Ri * (Li+1)^-1
	# find birth and death times for each pivot col. in reduced bdr i
	# find zero columns of bdr which are not nonzero rows of bdr1
	ntimesteps = maximum(grainvec[1])  # num of filtration steps is the maximum
	cols = ReducedDicta["Pcp"][i]
	cols1 = ReducedDicta["Pcp"][i+1]
	numCol = length(cols) - 1 # number of columns in the reduced matrix
	generators = Array{SparseMatrixCSC{Rational}}(undef, numCol)
	persistence = Array{Array{Rational{Int64},1}}(undef, numCol)
	simplexpersistence = Array{Array{Rational{Int64},1}}(undef, numCol)
	for j in 1 : numCol
		if cols[j] === cols[j+1] # if this is a zero column
			bi = ntimesteps + 1 - grainvec[i + 1][j] # birth time comes from grain vector for boundary i
			simplexbi = j
			# have to find the column with an entry in the jth row in order to find di
			if j in broadcast(abs, ReducedDicta["prowa"][i + 1])
				jval = findall(x -> abs(x) === j, ReducedDicta["prowa"][i + 1])
				jcol = maximum(broadcast(abs,ReducedDicta["pcola"][i + 1][jval]))
				di   = ntimesteps + 1 - grainvec[i + 2][jcol]
				simplexdi = jcol
			else
				di = Inf
			end
			if bi != di && di != Inf
				col = spzeros(Rational{Int64}, length(grainvec[i+1]))
				col[ReducedDicta["Prv"][i + 1][cran(ReducedDicta["Pcp"][i+1], jcol)]] = ReducedDicta["Pvl"][i + 1][cran(ReducedDicta["Pcp"][i + 1], jcol)]
				generators[j] = dropzeros(goodmatrix * sparse(col))
				persistence[j] = [bi, di]
				simplexpersistence[j] = [simplexbi, simplexdi]
			end
		end
	end
	defindices = filter(i -> isassigned(generators, i), 1:length(generators))
	generators = generators[defindices]
	persistence = persistence[defindices]
	simplexpersistence = simplexpersistence[defindices]
	return generators, persistence, simplexpersistence, goodmatrix
end

##############

function pointcloud_maxDim__bdrmatrices(pointcloud, distmat = false, maxdim = 1)
	"""
	takes in a pointcloud, returns bdr matrices dictionary, by default max dim = 1.
	"""
	if distmat
		C = eirene(pointcloud, maxdim = maxdim) # calls Eirene on the pointcloud and gets a dictionary C that contains useful information
	else
		C = eirene(pointcloud, maxdim = maxdim  , model = "pc") # calls Eirene on the pointcloud and gets a dictionary C that contains useful information
	end
	# retrieves information from C
	bdrMatrices = Dict{String, Vector}()  # storing bdr matrices info in a dictionary
	bdrMatrices["rv"] 	 	= 	Array{Array{Int64, 1}}(undef, maxdim + 1)
	bdrMatrices["cp"] 	 	= 	Array{Array{Int64, 1}}(undef, maxdim + 1)
	bdrMatrices["vl"] 	 	= 	Array{Array{Rational, 1}}(undef, maxdim + 1) # the value entries should be rational.
	bdrMatrices["m"] 	 	= 	Array{Int64}(undef, maxdim + 1) # m is the number of rows of the bdr matrix.
	bdrMatrices["grain"] 	= 	C["grain"][1:maxdim + 2] # grain vectors gives us information about the birth time of each simplex.
	bdrMatrices["lverts"] 	= 	Array{Array{Int64, 2}}(undef, maxdim + 1) # row labels
	bdrMatrices["hverts"] 	= 	Array{Array{Int64, 2}}(undef, maxdim + 1) # col labels
	bdrMatrices["distVec"] 	= 	C["ocg2rad"]
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

struct simplexhObject
		"""
		a custom object to store all the information we need for calculations below
		"""
		pointCloud :: Array{Float64, 2}
		bdrMatrices :: Dict{String, Vector}
		reducedbdr :: Dict{String, Vector}
		barCode :: Vector{Array{Array{Float64, 1}, 1}}
		simplexcode :: Vector{Array{Array{Float64, 1}, 1}}
		grainVec :: Vector{Array{Int64, 1}}
		distVec :: Vector{Float64}
		pcola :: Vector{Array{Int64, 1}}
		prowa :: Vector{Array{Int64, 1}}
		generators :: Vector{Array{SparseMatrixCSC, 1}}
		permutedlverts :: Vector{Array{Int64, 2}}
		permutedhverts :: Vector{Array{Int64, 2}}
end

function computeHomology(pointcloud, distmat = false, maxDim = 1)
	"""
	performs the full compute PH pipeline and creates a homologyObject (as above)
	which contains information about boundary matrices, generators, etc.
	"""
	bdrMatrices	 = pointcloud_maxDim__bdrmatrices(pointcloud, distmat, maxDim)
	filtrationorder  = bdrMatrices__filtrationorder(bdrMatrices, maxDim)
	reduced,grainvec,permutedlverts, permutedhverts = bdrMatrices__reducedbdr(filtrationorder)
	generators = Array{Array{SparseMatrixCSC}}(undef, maxDim)
	barCodes = Array{Array{Array{Rational{Int64},1}}}(undef, maxDim)
	simplexbarCodes = Array{Array{Array{Rational{Int64},1}}}(undef, maxDim)
	for i in 1: maxDim
		@time gens,bars, simplexbars= reduceddicta_i__gens_bars(reduced,grainvec,i)
		generators[i] = gens
		barCodes[i] = bars
		simplexbarCodes[i] = simplexbars
	end
	distVec = bdrMatrices["distVec"]
	# convert barcode to distances
	ldv1 = length(distVec) + 1
	distcode = Vector{Vector{Array{Float64,1}}}(undef, maxDim)
	simplexdistcode = Vector{Vector{Array{Float64,1}}}(undef, maxDim)
	for d in 1:maxDim
		distcode[d] = Vector{Array{Float64,1}}(undef,length(barCodes[d]))
		simplexdistcode[d] = Vector{Array{Float64,1}}(undef,length(simplexbarCodes[d]))
		for i in 1:length(barCodes[d])
			distcode[d][i] = [distVec[ldv1 - Int(barCodes[d][i][1])],distVec[ldv1 - Int(barCodes[d][i][2])]]
			simplexdistcode[d][i] = [simplexbarCodes[d][i][1], simplexbarCodes[d][i][2]]
		end
	end

	E = simplexhObject(
		pointcloud,
		filtrationorder,
		reduced,
		distcode,
		simplexdistcode,
		grainvec,
		distVec,
		reduced["pcola"],
		reduced["prowa"],
		generators,
		permutedlverts,
		permutedhverts)
	return E
end
