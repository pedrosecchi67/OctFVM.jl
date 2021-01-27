module OctFVM

    using LinearAlgebra
    using SparseArrays
    using IterTools

    include("OctMesh.jl")
    include("AdjacentMats.jl")

end # module
