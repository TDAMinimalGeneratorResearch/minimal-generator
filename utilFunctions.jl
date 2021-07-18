

######### visualize in 3D
function plotGenerators(C, d, k, plotAllGens = false)
	"""
	Takes in a homologyObject C,
	d   --> dimension of the generator to visualize
	k   --> the class number of the generator you want to plot
	all --> true if you want to plot all generators in one plot.
	"""
	gens = C.generators[d]
	pc = C.pointCloud
	if size(pc, 1) == 2
		is2D = true
	else
		is2D = false
	end
	if plotAllGens
		total = Array{Int64}(undef, length(gens))
		for i in 1:length(gens)
			total[i]= length(C.generators[d][i].rowval)
		end
	else
		total = Array{Int64}(undef, 1)
		total[1]= length(C.generators[d][k].rowval)
	end
	T = Array{GenericTrace{Dict{Symbol, Any}}}(undef, sum(total) + 1)
	if is2D
		pttype = "scatter"
		zpt = ones(1, size(pc, 2))
	else
		pttype = "scatter3d"
		zpt = pc[3,:]
	end
	T[1] = Plotly.scatter(;x=pc[1,:], y=pc[2,:],
					 z = zpt,
					 mode = "markers",
					 type = pttype,
					 name = "point cloud",
					 hoverinfo="skip",
					 maker_size = 10)
	count = 2
	for j in 1: length(total)
		ab = transpose(C.permutedlverts[d + 1][:,findall(x -> x!=0, C.generators[d][k])])
		xpts = pc[1,ab]
		ypts = pc[2,ab]
		if !is2D
			zpts = pc[3,ab]
		end
		for i in 1: total[j]
			if is2D
				zp = ones(1, 3)
			else
				zp = zpts[i,:] #append!(zpts[i,:], sum(zpts[i,:])/2)
			end
			T[count] = Plotly.scatter(
					  x = xpts[i,:], #append!(xpts[i,:], sum(xpts[i,:])/2),
					  y = ypts[i,:], #append!(ypts[i,:], sum(ypts[i,:])/2),
					  z = zp,
		              mode="lines+text",
					  type=pttype,
		              name=round(C.distVec[C.grainVec[d+1][i]]; digits=4),
					  fill="toself",
					  maker_size = 10,
					  hovertext=round(C.distVec[C.grainVec[d+1][i]]; digits=4),
					  hoverinfo = ["skip", "skip", "text"],
					  )
			count += 1
		end
	end
	if is2D
		zrange = [-1,1]
	else
		zrange = [1.2 * (minimum(pc[3,:]) - 0.01), 1.2 * (maximum(pc[3,:]) + 0.01)]
	end
	layout = Layout(;title="Generators!",
	                 xaxis_range=[1.2 * (minimum(pc[1,:]) - 0.01), 1.2 * (maximum(pc[1,:]) + 0.01)],
	                 yaxis_range=[1.2 * (minimum(pc[2,:]) - 0.01), 1.2 * (maximum(pc[2,:]) + 0.01)],
					 zaxis_range= zrange)
	return PlotlyJS.plot(T, layout)
end

############ Visualize various kinds of minimal generators

### minimal length (L1 optimal)
function plotMinimalLengthGenerators(C, d, gens)
	return plotMinimalGenerators(C,d,gens, "Minimal Length Generators!")
end

### minimal egde number (L0 optimal)
function plotMinimalEdgeGenerators(C,d,gens)
	return plotMinimalGenerators(C, d, gens, "Minimal Edge Number Generators!")
end

### optimal volume
function plotMinimalVolGenerators(C,d,gens)
	return plotMinimalGenerators(C, d, gens, "Minimal Volume Generators!")
end

### minimal area
function plotMinimalAreaGenerators(C, d, gens)
	return plotMinimalGenerators(C,d,gens, "Minimal Area Generators!")
end


### helper function to plot generators
function plotMinimalGenerators(C, d, gens, title)
	"""
	Takes in a homologyObject C,
	d   --> dimension of the generator to visualize
	k   --> the class number of the generator you want to plot
	all --> true if you want to plot all generators in one plot.
	"""
	# gens = C.generators[d]
	pc = C.pointCloud
	if size(pc, 1) == 2
		is2D = true
	else
		is2D = false
	end
	total = Array{Int64}(undef, 1)
	total[1]= length(gens)
	T = Array{GenericTrace{Dict{Symbol, Any}}}(undef, length(findall(x -> x!=0, gens)) + 1)
	if is2D
		pttype = "scatter"
		zpt = ones(1, size(pc, 2))
	else
		pttype = "scatter3d"
		zpt = pc[3,:]
	end
	T[1] = Plotly.scatter(;x=pc[1,:], y=pc[2,:],
					 z = zpt,
					 mode = "markers",
					 type = pttype,
					 name = "point cloud",
					 hoverinfo=[1:length(C.permutedlverts[1])],
					 maker_size = 10)
	count = 2
	for j in 1: length(total)
		ab = transpose(C.permutedlverts[d + 1][:,findall(x -> x!=0, gens)])
		xpts = pc[1,ab]
		ypts = pc[2,ab]
		if !is2D
			zpts = pc[3,ab]
		end

		for i in 1:size(ab,1)
			if is2D
				zp = ones(1, 3)
			else
				zp = append!(zpts[i,:], sum(zpts[i,:])/2)
			end
			T[count] = Plotly.scatter(
					  x = append!(xpts[i,:], sum(xpts[i,:])/2),
					  y = append!(ypts[i,:], sum(ypts[i,:])/2),
					  z = zp,
		              mode="lines+text",
					  type=pttype,
		              name=euclidean([xpts[i,1],ypts[i,1]], [xpts[i,2],ypts[i,2]]),
					  fill="toself",
					  maker_size = 10,
					  hovertext=euclidean([xpts[i,1],ypts[i,1]], [xpts[i,2],ypts[i,2]]),
					  hoverinfo = ["skip", "skip", "text"],
					  )
			count += 1
		end
	end
	if is2D
		zrange = [-1,1]
	else
		zrange = [1.2 * (minimum(pc[3,:]) - 0.01), 1.2 * (maximum(pc[3,:]) + 0.01)]
	end
	layout = Layout(;title=title,
	                 xaxis_range=[1.2 * (minimum(pc[1,:]) - 0.5), 1.2 * (maximum(pc[1,:]) + 0.01)],
	                 yaxis_range=[1.2 * (minimum(pc[2,:]) - 0.01), 1.2 * (maximum(pc[2,:]) + 0.01)],
					 zaxis_range= zrange)
	return PlotlyJS.plot(T, layout)
end


############ plot barcode with hover effects
function plotBarCode(C, maxdim = 1)
	"""
	Plot the bar code of the homology object up to dimension maxdim.
	"""
	pc = C.pointCloud
	numDim1 = length(C.barCode[1])
	if maxdim == 2
		numDim2 = length(C.barCode[2])
		total = numDim1 + numDim2
	else
		total = numDim1
	end
	T = Array{GenericTrace{Dict{Symbol, Any}}}(undef, total)
	y = 1
	for j in 1: (length(C.barCode[1])-1)
		T[y] = Plotly.scatter(
				  x = C.barCode[1][j],
				  y = [y,y],
				  mode="lines+text",
				  type="scatter",
				  line=attr(color="blue", width=1.5),
				  fill="toself",
				  maker_size = 10,
				  showlegend=false,
				  hovertext=string("class", j),
				  hoverinfo = ["skip", "text"],
				  )
		y += 1
	end
	T[y] = Plotly.scatter(
			  x = C.barCode[1][length(C.barCode[1])],
			  y = [y,y],
			  mode="lines+text",
			  type="scatter",
			  line=attr(color="blue", width=1.5),
			  fill="toself",
			  maker_size = 10,
			  name = "dim 1",
			  hovertext=string("class", y),
			  hoverinfo = ["skip", "text"],
			  )
	t = 1
	if maxdim == 2
		for k in 1:(length(C.barCode[2]) - 1)
			T[t + y] = Plotly.scatter(
					  x = C.barCode[2][k],
					  y = [y + t, y + t],
					  mode="lines+text",
					  type="scatter",
					  line=attr(color="red", width=1.5),
					  fill="toself",
					  maker_size = 10,
					  showlegend=false,
					  hovertext=string("class", k),
					  hoverinfo = ["skip", "text"],
					  )
			t += 1
		end
		print("sth")
		printval(C.barCode[2][t], "start")
		T[t + y] = Plotly.scatter(
				  x = C.barCode[2][t],
				  y = [y + t,y + t],
				  mode="lines+text",
				  type="scatter",
				  line=attr(color="red", width=1.5),
				  fill="toself",
				  maker_size = 10,
				  name="dim 2",
				  hovertext=string("class", t),
				  hoverinfo = ["skip", "text"],
				  )
	end
	layout = Layout(;title="Bar Code!",
					 xaxis_range=[0, maximum(C.distVec)],
					 yaxis_range=[0, y + t+ 1],
					 )
	return PlotlyJS.plot(T, layout)
end

############ helper function to plot 2D plots

function plotGeneratorin2D(C, d, k, plotAllGens = false)
	"""
	Takes in a homologyObject C,
	d   --> dimension of the generator to visualize
	k   --> the class number of the generator you want to plot
	all --> true if you want to plot all generators in one plot.
	"""
	gens = C.generators[d]
	pc = C.pointCloud

	if plotAllGens
		total = Array{Int64}(undef, length(gens))
		for i in 1:length(gens)
			total[i]= length(C.generators[d][i].rowval)
		end
	else
		total = Array{Int64}(undef, 1)
		total[1]= length(C.generators[d][k].rowval)
	end
	T = Array{GenericTrace{Dict{Symbol, Any}}}(undef, sum(total) + 1)

	T[1] = Plotly.scatter(;x=pc[1,:], y=pc[2,:],
					 mode = "markers",
					 type = "scatter",
					 name = "point cloud",
					 hoverinfo="skip")
	count = 2
	for j in 1: length(total)
		ab = transpose(C.permutedlverts[d + 1][:,findall(x -> x!=0, gens[k])])
		xpts = pc[1,ab]
		ypts = pc[2,ab]
		for i in 1: total[j]
			T[count] = Plotly.scatter(
					  x = append!(xpts[i,:], sum(xpts[i,:])/2),
					  y = append!(ypts[i,:], sum(ypts[i,:])/2),
		              mode="lines+text",
		              name=round(C.distVec[C.grainVec[d+1]][i]; digits=4),
					  fill="toself",
					  hovertext=round(C.distVec[C.grainVec[2]][i]; digits=4),
					  hoverinfo = ["skip", "skip", "text"]
					  )
			count += 1
		end
	end
	layout = Layout(;title="Generators!",
	                 xaxis_range=[1.2 * (minimum(pc[1,:]) - 0.01), 1.2 * (maximum(pc[1,:]) + 0.01)],
	                 yaxis_range=[1.2 * (minimum(pc[2,:]) - 0.01), 1.2 * (maximum(pc[2,:]) + 0.01)],
					 # zaxis_range=[-10,10]
					 )
	PlotlyJS.plot(T, layout)
end
############ visualize in 2D

function plotGeneratorsin2D(C, d, k, labelEdge = true)
	ab = transpose(C.permutedhverts[1][:, C.generators[1][k].rowval])
	xpts = C.pointCloud[1,ab]
	ypts = C.pointCloud[2,ab]
	p = Plots.plot(xpts[1,:],ypts[1,:], legend = false, title = "Generator!")
	if labelEdge
		annotate!(p, [((xpts[1,:][1]+xpts[1,:][2])/2, (ypts[1,:][1] + ypts[1,:][2])/2 , 1)])
	end
	for i in 2:length(ab[:,1])
		plot!(p,xpts[i,:],ypts[i,:] )
		if labelEdge
			annotate!(p, [((xpts[i,:][1]+xpts[i,:][2])/2, (ypts[i,:][1] + ypts[i,:][2])/2 , i)])
		end
	end
	p
	scatter!(p,C.pointCloud[1,:],C.pointCloud[2,:], dpi=300)
end

############ PCA (not used yet)
# M = fit(PCA, ex; maxoutdim=100)
# Yte = MultivariateStats.transform(M, ex)
# ezplot_pjs(Yte)
# ezplot_pjs(ex)
# Xr = reconstruct(M, Yte)


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
	pc = C.pointcloud
	if m > size(pc, 2) || m < 0
		println("Invalid point index. ")
		return
	end
	return C.pointcloud[:, m]
end

function C_d__barCode(C, d, plotbars = true, plotcolor = "red")
	"""
	* (C,d)	      --> dth barcode
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
	return barCode[d]
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
	return grainVec[d][n]
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
	return C.pointCloud[:,  C.permutedlverts[d][:,n]]
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
	return prowa_D[k]
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
	generators = C.generators[d]
	if l > length(generators) || l <= 0
		println("Invalid class number.")
		return
	end
	return C.generators[d][l]
end


function C_d_l__filt(C,d,l, pivot = true)
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
	if pivot
		# get pivot column indices born before birth time of the generator
		pcola = C.pcola[d+1]
		pcola = broadcast(abs, pcola)
		pcolBeforeLastCol = pcola[findall(x-> x <= lastcol, pcola)]
		sort!(pcolBeforeLastCol)
		pivotcp = [1]
		pivotvl = Vector{Int64}(undef, 0)
		pivotrv = Vector{Int64}(undef, 0)
		for col in pcolBeforeLastCol
			colEnd = C.bdrMatrices["cp"][d+1][col+1]
			colStart =  C.bdrMatrices["cp"][d+1][col]
			append!(pivotcp, pivotcp[length(pivotcp)] + colEnd - colStart)
			append!(pivotvl, C.bdrMatrices["vl"][d+1][colStart: (colEnd - 1)])
			append!(pivotrv, C.bdrMatrices["rv"][d+1][colStart: (colEnd - 1)])
		end
		newcp = pivotcp
		newrv = pivotrv
		newvl = pivotvl
	else
	# remove columns born later than lastcol
		newcp = C.bdrMatrices["cp"][d+1][1:lastcol+1]
		newrv = C.bdrMatrices["rv"][d+1][1:newcp[lastcol + 1] - 1]
		newvl = C.bdrMatrices["vl"][d+1][1:newcp[lastcol + 1] - 1]
	end
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


function filterSparseMatrix(C, d, firstcol, lastcol)
	"""
	Given begin column index and end column index, filter out columns
	with a birth time between the first and last column, including the first
	and last column.
	C -> Homology object
	d -> dimension of bdr matrix to filter
	first col, last col -> indices of the bounding columns
	"""
# keep columns born between first col and lastcol, inclusive
	newcp = C.bdrMatrices["cp"][d][firstcol : lastcol+1]
# need to make sure colptr starts with 1.
	newcp = newcp .- (newcp[1] - 1)
# get new rv, vl, and m(row number) values.
	cpr = C.bdrMatrices["cp"][d]
	rowvalFirst = cpr[firstcol]
	rowvalLast  = cpr[lastcol + 1] - 1
	newrv  		= C.bdrMatrices["rv"][d][rowvalFirst : rowvalLast]
	newvl 		= C.bdrMatrices["vl"][d][rowvalFirst : rowvalLast]
	m = C.bdrMatrices["m"][d]
# return a sparse matrix with the columns of d-dimension boundary matrix born
# in [firstcol, lastcol]
	filt = SparseMatrixCSC(m,length(newcp) - 1,newcp,newrv,newvl)
end


function findArea(points)
	"""
	Takes in an array with three columns representing the coordinates
	of the triangle, and return the area of that triangle.
	"""
	if size(points, 2) != 3
		print("Invalid points input. ")
		return
	end
	a = euclidean(points[:,1], points[:,2])
	b = euclidean(points[:,1], points[:,3])
	c = euclidean(points[:,2], points[:,3])
	s = (a + b + c)/2
	area = sqrt(s * (s-a) * (s-b) * (s - c))
	return area
end


function findAreaOfPolygon(D, gens)
	"""
	Given a Homology object and one of its generators,
	return the area that the generator encloses.
	"""
	if length(findall(x->x!=0, gens)) == 0
		return 0.0
	end
	print(length(findall(x->x!=0, gens)))
	gen = findall(x->x!=0, gens)
	pts = D.permutedhverts[1][:, gen]
	orderedVerts = orderPolygonVertices(pts)
	x = D.pointCloud[:,unique(orderedVerts)][1,:]
	y = D.pointCloud[:,unique(orderedVerts)][2,:]
	areaA = 0.0
	j = length(gen)
	for i in 1: j
		areaA = areaA + (x[j] + x[i]) * (y[j] - y[i])
		j = i
	end
	return abs(areaA/2.0)
end

function swap(pts, i, j)
	"""
	Helper function to swap two columns of a matrix.
	"""
	temp = pts[:,i]
	pts[:,i] = pts[:,j]
	pts[:,j] = temp
	return pts
end

function orderPolygonVertices(pts)
	"""
	Order the vertices of a generator so that we can find the area of the polygon.
	"""
	target = pts[:,1][2]
	for i in 1: (size(pts, 2) - 1)
		print(target)
		idx = findall(x-> x == target, pts[:,i + 1:size(pts,2)])
		swap(pts, i + 1, i + idx[1][2])
		if idx[1][1] == 1
			target = pts[:, i + 1][2]
		else
			target = pts[:, i + 1][1]
		end
	end
	pts
	orderedVerts = unique(pts)
end

function checkCircuit(generator)
	"""
	Given a generator, return a boolean representing if the generator is a
	circuit or not.
	"""
	gen = findall(x->x!=0, generator)
	pts = C1.permutedhverts[1][:, gen]
	first = pts[:,1][1]
	target = pts[:,1][2]
	last = -1
	for i in 1: (size(pts, 2) - 1)
		idx = findall(x-> x == target, pts[:,i + 1:size(pts,2)])
		if length(idx) == 0
			return false
		end
		swap(pts, i + 1, i + idx[1][2])
		if idx[1][1] == 1
			target = pts[:, i + 1][2]
			last = pts[:, i + 1][1]
		else
			target = pts[:, i + 1][1]
			last = pts[:, i + 1][2]
		end
	end
	if last in pts[:, size(pts, 2)] && first in pts[:, size(pts, 2)]
		return true
	else
		return false
	end
	return orderedVerts = unique(pts)
end


function C_gen_checkLoops(C, g)
      D=SparseMatrixCSC(C.bdrMatrices["m"][1], length(C.bdrMatrices["cp"][1]) - 1,C.bdrMatrices["cp"][1],C.bdrMatrices["rv"][1],C.bdrMatrices["vl"][1])
      T = D[:,g.rowval]
      return length(g.rowval) - rank(Array(T))
end
