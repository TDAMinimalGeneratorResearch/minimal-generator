# minimal-generator


# Compute Homology

### Requirements:
1. Download Julia from https://julialang.org/downloads/ and double click to install.
2. Download Atom from https://atom.io/.
3. Open Atom, go to preference -> packages, install package: `uber-juno`.
4. Download Gurobi and obtain a free academic license at https://www.gurobi.com/academia/academic-program-and-licenses/. Open a Julia cones ENV["GUROBI_HOME"] = "/Library/gurobi912/mac64"
5. Open the folder containing all the scripts in Atom.
6. Open the file named `installRequirements.jl`, run this file to install required packages and load all the other scripts into the environment. I had to run Pkg.build("HDF5") to compile Eirene.


### To compute homology and optimize generators for a given pointcloud, look at an example pipeline in the file `exampleRun.jl`.


### Some output functions:
* point cloud --> C, where C is some object with "all the information you need"
					about the computation, we used a customized object called "homologyObject" to represent C.
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
