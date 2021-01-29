Pkg.add(["MultivariateStats", "Plots", "Plotly", "DelimitedFiles", "StatsBase", "Distances", "JuMP", "Clp", "GLPK", "PlotlyJS", "JLD", "ORCA", "LinearAlgebra", "Eirene", "SparseArrays", "Random", "Distributions", "Gurobi"])
using MultivariateStats, Plots, Plotly, DelimitedFiles, StatsBase, Distances, JuMP, PlotlyJS, JLD, ORCA, LinearAlgebra, Eirene, SparseArrays, Random, Distributions
include("implementComputeHomology.jl")
include("utilFunctions.jl")
include("UniformAndLengthOptimalCycles.jl")
include("dataGenerator.jl")
include("outputFunctions.jl")
