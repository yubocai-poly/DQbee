import Pkg

# add all the packages needed
Pkg.add("ReachabilityAnalysis")
Pkg.add("CarlemanLinearization")
Pkg.add("Plots")
Pkg.add("LaTeXStrings")
Pkg.add("LinearAlgebra")
Pkg.add("SparseArrays")
Pkg.add("LazySets")
Pkg.add("CarlemanLinearization")

# make sure all packages are up-to-date
Pkg.update()