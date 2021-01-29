############ ouput functions

"""
* point cloud --> C, where C is some object with "all the information you need"
					about the computation (in Eirene this is implemented as a
					dictionary/struct, but there are other options worth thinking
					about, such as customized objects)
* (C,m) 	  --> the euclidean coordinates of the mth point in the point cloud
					(assuming the input is a point cloud)
* (C,d) 	  --> dth boundary matrix
* (C,d)	      --> dth barcode
* (C,d,n)	  --> birth time of nth cell in dimension d
* (C,d,n) 	  --> vertices of nth cell in dimension d
* (C,d,k) 	  --> row index of the kth pivot element for the boundary matrix in dim d
* (C,d,k) 	  --> col index of the kth pivot element for the boundary matrix in dim d
* (C,d,l) 	  --> generator for lth persistent homology class
					(expressed as a linear combination of simplices; concretely,
					this would probably take the form of a 1-column sparse matrix)
"""

function C_m__coord(C, m)
	"""
	* (C,m) 	  --> the euclidean coordinates of the mth point in the point cloud
						(assuming the input is a point cloud)
	"""
	pc = C.pointCloud
	if m > size(pc, 2) || m < 0
		println("Invalid point index. ")
		return
	end
	return C.pointCloud[:, m]
end

function C_d__bdrmatrix(C, d)
	"""
	* (C,d) 	  --> the boundary matrix for dimension d
	"""
    if d > length(C.bdrMatrices["m"]) || d <= 0
        println("Invalid dimension d.")
        return
    end
    bdrs = C.bdrMatrices
    return SparseMatrixCSC(bdrs["m"][d], length(bdrs["cp"][d]) - 1, bdrs["cp"][d], bdrs["rv"][d], bdrs["vl"][d])
end

function C_d_l__barCode(C, d, l, plotbars = true, plotcolor = "red")
	"""
	* (C,d,l)	      --> lth barcode in dimension d
	  plotbars 		  --> plot the barcode for dimension d,
	  					  can be turned off by passing in false.
	"""
	barCode = C.barCode
	if d > length(barCode) || d < 0
		println("Invalid barcode index. ")
		return
	end
	if plotbars
		p = Plots.plot(title = "Barcode", legend = false)
		for i in 1:length(barCode[d])
			Plots.plot!(p,barCode[d][i],[i,i],linecolor = plotcolor)
		end
		display(p)
	end
	return barCode[d][l]
end


function C_d__m_n_rv_cp_vl(C, d)
	"""
	* (C,d) 	  --> dth boundary matrix
	"""
	num = length(C.bdrMatrices["m"])
	if d > num || d < 0
		println("Invalid boundary matrix index. ")
		return
	end
	m = C.bdrMatrices["m"][d]
	cp = C.bdrMatrices["cp"][d]
	rv = C.bdrMatrices["rv"][d]
	vl = C.bdrMatrices["vl"][d]
	n = length(cp) - 1
	return SparseMatrixCSC(m, n, cp, rv, vl)
end

function C_d_n__birthTime(C, d, n)
	"""
	* (C,d,n)	  --> birth time of nth cell in dimension d
	"""
	grainVec = C.grainVec
	if d > length(grainVec)
		println("Invalid dimension input. ")
		return
	end
	d_grainVec = grainVec[d]
	if n > length(d_grainVec) || n < 0
		println("Invalid cell index. ")
		return
	end
	return grainVec[1][1] - grainVec[d][n]
end

function C_d_n__vertices(C, d, n)
	"""
	* (C,d,n) 	  --> vertices of nth cell in dimension d
	"""
	if d > length(C.bdrMatrices["m"]) || d < 0
		println("Invalid dimension input. ")
		return
	end
	verts = C.permutedlverts[d]
	if n > size(verts, 2) || n < 0
		println("Invalid cell index. ")
		return
	end
	if d < length(C.permutedlverts)
		return C.pointCloud[:,  C.permutedlverts[d][:,n]]
	end
	if d == length(C.permutedlverts)
		return C.pointCloud[:,  C.permutedhverts[d][:,n]]
	end
end


function C_d_k__prowind(C, d, k)
	"""
	* (C,d,k) 	  --> row index of the kth pivot element for the boundary matrix in dim d
	"""
	prowa = C.prowa
	if d > length(prowa) || d < 0
		println("Invalid input dimension number d. ")
		return
	end
	prowa_d = C.prowa[d]
	if k > length(prowa_d) || k < 0
		println("Invalid row index k.")
		return
	end
	return prowa_d[k]
end

function C_d_k__pcolind(C, d, k)
	"""
	* (C,d,k) 	  --> col index of the kth pivot element for the boundary matrix in dim d
	"""
	pcola = C.pcola
	if d > length(pcola) || d < 0
		println("Invalid input dimension number d. ")
		return
	end
	pcola_d = C.pcola[d]
	if k > length(pcola_d) || k < 0
		println("Invalid column index k.")
		return
	end
	return pcola_d[k]
end



function C_d_l__gen(C, d, l)
	"""
	* (C,d,l) 	  --> generator for lth persistent homology class in dim. d
					(expressed as a linear combination of simplices; concretely,
					this would probably take the form of a 1-column sparse matrix)
    """
	if d > length(C.generators)
		println("Invalid dimension d.")
		return
	end
	generators = C.generators[d]
	if l > length(generators) || l <= 0
		println("Invalid class number.")
		return
	end
	return C.generators[d][l]
end


function C_d_l__filt(C,d,l)
	"""
	* (C,d,l) 	  --> sparse matrix equal to the non-reduced boundary d+1 matrix,
					with entries before the birth time of the lth d-dimensional
					generator.
	"""
	# last nonzero entry of lth generator tells us its birth time
	ind = maximum(findall(x -> x != 0, Array(C_d_l__gen(C,d,l))))
	# we need columns of boundary d + 1 which are born before bgrain, the birth
	# time of the generator
	bgrain = C.grainVec[d+1][ind]
	ab = length(findall(x -> x >= bgrain, C.grainVec[d+2]))
	if ab == 0
		m = C.bdrMatrices["m"][d+1]
		return sparse(zeros(m,1))
	end
	lastcol = maximum(findall(x -> x >= bgrain, C.grainVec[d+2]))
	# remove columns born later than lastcol
	newcp = C.bdrMatrices["cp"][d+1][1:lastcol+1]
	newrv = C.bdrMatrices["rv"][d+1][1:newcp[lastcol + 1] - 1]
	newvl = C.bdrMatrices["vl"][d+1][1:newcp[lastcol + 1] - 1]
	m = C.bdrMatrices["m"][d+1]
	# return a sparse matrix with the columns of boundary i+1 born after the
	# birth time of the generator
	filt = SparseMatrixCSC(m,length(newcp) - 1,newcp,newrv,newvl)
	return filt
end


function C_d_l__filt_btw(C,d,l,q)
	"""
	* (C,d,l) 	  --> sparse matrix equal to the non-reduced boundary d+1 matrix,
					with entries after the birth time of the lth d-dimensional
					generator and before the death time of the lth d-dimensional
					generator.
	"""
	bdist, ddist = C.barCode[d][l]
	bgrain = findall(x -> x == bdist, C.distVec)[1]
	dgrain = findall(x -> x == ddist, C.distVec)[1]
	# we need columns of boundary d + 1 which are born after btime, and before dtime
	ab = length(findall(x -> dgrain < x < bgrain, C.grainVec[q+2]))
	if ab == 0
		m = C.bdrMatrices["m"][q+1]
		return sparse(zeros(m,1))
	end
	firstcol = minimum(findall(x -> x < bgrain, C.grainVec[q+2]))
	lastcol  = maximum(findall(x -> x > dgrain, C.grainVec[q+2]))
# keep columns born between first col and lastcol
# return a sparse matrix with the columns of boundary i+1 born after the
# birth time of the generator
	filt = filterSparseMatrix(C, q+1, firstcol, lastcol)
	return filt
end
