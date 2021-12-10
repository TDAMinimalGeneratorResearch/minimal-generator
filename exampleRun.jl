
"""
This file demonstrates a sample pipeline which takes in a pointcloud csv file as input,
we compute the homology as well as its dimension one minimal generators.
"""
##########################################
#              COMPUTE HOMOLOGY          #
##########################################
# 1. read in file:
pc = readdlm("data/Synthetic-data/Gamma/Gamma-pointcloud/2x100-Gamma-4.csv")
# 2. call computeHomology
C = computeHomology(pc, false, 1) # <- compute homology of pc in dimension 1
plotBarCode(C) # plots the bar code for a pointcloud
# 3.Some outputs we might be interested in:
m_th_point = C_m__coord(C, 1)
# gets the dth boundary matrix
d_bdrmatrix = C_d__bdrmatrix(C, 1)
# Plots the barcode
d_barcode = C_d_l__barCode(C, 1, 1)
# (C,d,n)	  --> birth time of nth cell in dimension d get the birth time of the nth cell in dimension d.
birth_time = C_d_n__birthTime(C, 2, 2)
# (C,d,n) 	  --> vertices of nth cell in dimension d
vertices = C_d_n__vertices(C, 2, 2)
# (C,d,k) 	  --> row index of the kth pivot element for the boundary matrix in dim d
rowind = C_d_k__prowind(C, 2, 14)
# (C,d,k) 	  --> col index of the kth pivot element for the boundary matrix in dim d
colind = C_d_k__pcolind(C, 2, 14)
C.pcola
# (C,d,l) 	  --> generator for lth persistent homology class
#           (expressed as a linear combination of simplices; concretely,
#           this would probably take the form of a 1-column sparse matrix)
gen = C_d_l__gen(C, 1, 10)

##########################################
#              MINIMAL GENERATOR         #
##########################################

Pkg.build("Gurobi") # run in repl: ENV["GUROBI_HOME"]="/Library/gurobi903/mac64/"
using Gurobi
d = 1  # dimension of the generator we hope to minimize
requireIntegralSol = false # whether we want to require the generator vector to have integral entries.



"""
output:
1. minimal generator vector
2. the length weighted cost of the optimal generator
3. the uniform weighted cost of the optimal generator
"""
##########################################
#              EDGE- LOSS                #
##########################################


## 1. Uniform-weighted edge loss methods, which minimizes the number of edges in the generators.

##################
# To optimize a basis of cycle representatives
##################
function prsb_unif_Edge(C)
    optimized = Array{SparseMatrixCSC{Float64, Int64}, 1}(undef, length(C.generators[1]))
    for l in 1: length(C.generators[1])
        uniform_weighted_minimal_gen, uniform_gen_len, uniform_zeroNorm = findUnifEdgeOptimalCycles_prs(C, d, l, optimized, requireIntegralSol, true, true)
        optimized[l] = vcat(uniform_weighted_minimal_gen, zeros(length(C.generators[1][1]) - length(uniform_weighted_minimal_gen.nzval)))
    end
    return optimized
end

## 2. Length-weighted edge loss methods, which minimizes the total length of edges in the generators.

function prsb_length_Edge(C)
    optimized = Array{SparseMatrixCSC{Float64, Int64}, 1}(undef, length(C.generators[1]))
    for l in 1: length(C.generators[1])
        uniform_weighted_minimal_gen, uniform_gen_len, uniform_zeroNorm = findLenEdgeOptimalCycles_prs(C, d, l, optimized, requireIntegralSol, true, true)
        optimized[l] = vcat(uniform_weighted_minimal_gen, zeros(length(C.generators[1][1]) - length(uniform_weighted_minimal_gen.nzval)))
    end
    return optimized
end

# Example calls:
len_weighted_edge = prsb_unif_Edge(C)
len_weighted_edge[l] # to access the lth optimized generator

unif_weighted_edge = prsb_length_Edge(C)
unif_weighted_edge[l] # to access the lth optimized generator


##################
# To optimize a single generator
##################
l = 17 # index of the generator we hope to minimize
len_weighted = false # can substitute with a vector of lengths to customize your own way of defining length of edges
uniform_weighted_minimal_gen, uniform_Len, uniform_zeroNorm = C_d_l__minimal(C,d,l, len_weighted, false, C.generators[1][1:l-1], false)

len_weighted = true
length_weighted_minimal_gen, len_Len, len_zeroNorm = C_d_l__minimal(C,d,l, true, len_weighted, C.generators[1][1:l-1], false)

# plot only works up to dimension 3.
plotMinimalEdgeGenerators(C,1,uniform_weighted_minimal_gen) # plots the optimal generator
plotGenerators(C,d,l) # plots the original generator


##########################################
#         TRIANGLE - LOSS                #
##########################################

## 3. Uniform-weighted triangle loss methods, which minimizes the number of triangles the generator bounds.

function unif_tri(C)
    optimized = Array{SparseMatrixCSC{Float64, Int64}, 1}(undef, length(C.generators[1]))
    for l in 1: length(C.generators[1])
        dd, triVerts = constructInput(C,1,l)
        optimal = findVolumeOptimalCycle(C, 1, l, dd, triVerts, true)
        optimized[l] = optimal[1]
    end
    return optimized
end

## 4. Area-weighted triangle loss methods, which minimizes the total area of triangles the generator bounds.

function area_tri(C)
    optimized = Array{SparseMatrixCSC{Float64, Int64}, 1}(undef, length(C.generators[1]))
    for l in 1: length(C.generators[1])
        dd, triVerts = constructInput(C,1,l)
        optimal = findAreaOptimalCycle(C, 1, l, dd, triVerts, true)
        optimized[l] = optimal[1]
    end
    return optimized
end

# Example calls
area_tri(C)

unif_tri(C)

# Optimize a single generator
dd, triVerts = constructInput(C,d,l)
optimal_cycle = findAreaOptimalCycle(C, d, l, dd, triVerts, true)[1]
optimal_volume = findAreaOptimalCycle(C, d, l, dd, triVerts, true)[2] # number of triangles this cycle bounds

# save this C object for future use.
save("sampleHObject.jld", "C", C)
# load saved C objects.
C = load("sampleHObject.jld")["C"]
