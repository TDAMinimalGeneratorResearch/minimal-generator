Pkg.add(["MultivariateStats", "Plots", "Eirene", "Plotly", "DelimitedFiles", "StatsBase", "Distances", "JuMP", "Clp", "GLPK", "PlotlyJS", "JLD", "ORCA", "LinearAlgebra", "Eirene", "SparseArrays", "Random", "Distributions", "Gurobi"])
using MultivariateStats,  Distances, Plotly, DelimitedFiles, StatsBase, Distances, JuMP, PlotlyJS, JLD, ORCA, LinearAlgebra, Eirene, SparseArrays, Random, Distributions, SparseArrays, Eirene, GLPK
using Plots
include("computePH.jl")
include("utilFunctions.jl")
include("edge-loss.jl")
include("triangle-loss.jl")
include("outputFunctions.jl")
