"""
This file demonstrates a sample pipeline which takes in a pointcloud csv file as input,
we compute the homology as well as its dimension one minimal generators.
"""
##########################################
#              COMPUTE HOMOLOGY          #
##########################################
# 1. read in file:
pc = readdlm("data/Synthetic-data/Gamma/Gamma-pointcloud/8x100-Gamma-8.csv")
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
l = 17 # index of the generator we hope to minimize
requireIntegralSol = true # whether we want to require the generator vector to have integral entries.
## 1. Uniform-weighted minimal generator: (minimizing the number of edges in the generators)
"""
output:
1. minimal generator vector
2. the length weighted cost of the optimal generator
3. the uniform weighted cost of the optimal generator
"""

function prsb_Edge(C)
    optimied = Array{SparseMatrixCSC{Float64, Int64}, 1}(undef, length(C.generators[1]))
    for l in 1: length(C.generators[1])
        uniform_weighted_minimal_gen, uniform_gen_len, uniform_zeroNorm = findUnifEdgeOptimalCycles_prs(C, d, l, optimied, requireIntegralSol, true, true)
        optimied[l] = vcat(uniform_weighted_minimal_gen, zeros(length(C.generators[1][1]) - length(uniform_weighted_minimal_gen.nzval)))
    end
    return optimied
end

findEdge(C)


plotMinimalEdgeGenerators(C,1,uniform_weighted_minimal_gen) # plots the optimal generator
plotGenerators(C,d,l) # plots the original generator
## requiring integral solutions
requireIntegralSol = false
uniform_weighted_minimal_gen, uniform_Len, uniform_zeroNorm = findUniformWeightedOptimalCycles(C, d, l, requireIntegralSol)

function tri_loss(C)
    optimied = Array{SparseMatrixCSC{Float64, Int64}, 1}(undef, length(C.generators[1]))
    for l in 1: length(C.generators[1])
        lower, upper = C.barCode[1][l]
        dd, hverts=  findbasis(C, C.pointCloud, lower, upper, false)
        optimal = findVolumeOptimalCycle(C,d,l,dd, hverts)
        optimied[l] = optimal[1]
    end
    return optimied
end

tri_loss(C)

plotMinimalEdgeGenerators(C,1,uniform_weighted_minimal_gen)

# 2. length-weighted minimal generator: (minimizing the length of the generators)
length_weighted_minimal_gen, length_Len, length_zeroNorm = findLengthWeightedOptimalCycles(C, d, l, requireIntegralSol)
plotMinimalLengthGenerators(C,1,length_weighted_minimal_gen)
# save this C object for future use.
save("sampleHObject.jld", "C", C)
# load saved C objects.
C = load("sampleHObject.jld")["C"]
