# Read in sparse formatted data
data = readdlm("Sociology_sparse_matrix.csv", ',')
rowIdx = Vector{Int64}(data[:, 1][2:size(data)[1]].+1)
colIdx = Vector{Int64}(data[:, 2][2:size(data)[1]].+1)
val = Vector{Float64}(data[:, 3][2:size(data)[1]])
distmat = sparse(rowIdx, colIdx, val, maximum(rowIdx), maximum(rowIdx))

# Format it into a sparse matrix and then convert to dense to change all 0.0 to 1.2
arr = Array(distmat)
arr[arr.==0] .= 1.2

# convert back the diagonal entries to be 0
for i in 1:size(arr)[1]
    arr[i,i] = 0
end

# compute persistence
C = eirene(arr, maxdim = 1) # calls Eirene on the pointcloud and gets a dictionary C that contains useful information

# plot bar code
plotbarcode_pjs(C, dim = 1)

# see the longer bars and then get the edges in the loop
i = ??  # set i to be the generator id you want to look at
C["cyclerep"][3][i]

# get the edges of the generator
edges = classrep(C, dim=1, class=i)

# read in indices file
concepts = readdlm("Sociology_concept_indices.csv", ',')

edgeList = Array{String,2}(concepts[:,2][edges.+1]) # adding 1 because Eirene index starts at 1 and concepts file index starts at 0.

writedlm("edges_sociology.csv", concepts[:,2][edges .+ 1])


writedlm("sociology_distance_mat.csv", arr)



################
# Find the (birth, death) edge pair, i.e., the edge born at the birth time of the cycle and the
# edge born at the death time of the cycle.
################
function get_birth_death_edge_pair(C, d, l)
    # get the index of the death triangle
    death_triangle_idx = Int(C.simplexcode[1][l][2])
    # get the edges of the death triangle
    edges_of_death_triangle = C.bdrMatrices["hverts"][2][:, death_triangle_idx]
    # get the edge born last in the death triangle, this is the edge that kills the cycle
    deathEdge = C.bdrMatrices["hverts"][1][:, maximum(edges_of_death_triangle)]
    # get the edge born at the birth time of the cycle
    birth_edge_idx = Int(C.simplexcode[1][l][1])
    birthEdge = C.bdrMatrices["hverts"][1][:, birth_edge_idx]
    return (birthEdge, deathEdge)
end

### This function helps to retrieve the birth and death concepts given a concepts file.
function get_birth_death_concept(C, d, l, concepts)
    birth_edge, death_edge = get_birth_death_edge_pair(C, d, l)
    birthConcept = Vector{String}(concepts[:,2][birth_edge.+1])
    deathConcept = Vector{String}(concepts[:,2][deathEdge.+1])
    return (birthConcept, deathConcept)
end


birth_edge, death_edge = get_birth_death_edge_pair(C, d, l)

concepts = readdlm("Sociology_concept_indices.csv", ',')
birthConcept, deathConcept = get_birth_death_concept(C, d, l, concepts)
