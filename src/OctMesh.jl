abstract type AbstractOctNode end

"""
Struct to hold an OctTree node

* `subnodes`: matrix of subnodes
* `father`: `AbstractOctNode` or `Nothing`, for the node from which the object branches
* `lims`: geometric limits of the node in each axis (shape `(N, 2)`)
* `adjacent`: matrix of `AbstractOctNode` objects referring to `OctTree` nodes and leaves,
    with the same shape as `lims`
"""
mutable struct OctNode{N} <: AbstractOctNode
    subnodes::Array{AbstractOctNode, N}
    father::Union{AbstractOctNode, Nothing}
    lims::Matrix{Float64}
    adjacent::Matrix{Union{AbstractOctNode, Nothing}}
end

"""
Struct to hold an OctTree leaf

* `ind`: index for a leaf's corresponding degree of freedom. Defaults to zero in constructors, 
    and is set by `set_indexing!`.
    Can also refer to an array of leaf indices resulting from the merging of a branch
* `father`: `AbstractOctNode` or `Nothing`, for the node from which the leaf branches
* `lims`: geometric limits of the node in each axis (shape `(N, 2)`):
* `adjacent`: matrix of `AbstractOctNode` objects referring to `OctTree` nodes and leaves,
    with the same shape as `lims`
"""
mutable struct OctLeaf <: AbstractOctNode
    ind::Union{Int64, AbstractArray}
    father::Union{AbstractOctNode, Nothing}
    lims::Matrix{Float64}
    adjacent::Matrix{Union{AbstractOctNode, Nothing}}
end

"""
Constructor for OctLeaf

* `lims`: geometric limits for leaf in each axis(shape `(ndims, 2)`)
"""
OctLeaf(lims::Matrix{Float64}) = OctLeaf(
    0,
    nothing,
    lims,
    Union{AbstractOctNode, Nothing}[
        nothing for i=1:size(lims, 1), j=1:2
    ]
)

"""
Get area for leaf

* `lf`: leaf at hand
* `dim`: dimension to get area in
"""
function leaf_area(lf::AbstractOctNode, dim::Int64)

    a=1.0

    for i=1:size(lf.lims, 1)
        if i!=dim
            a*=(lf.lims[i, 2]-lf.lims[i, 1])
        end
    end

    a

end

"""
Get volume for leaf

* `lf`: leaf at hand
"""
function leaf_volume(lf::AbstractOctNode)

    v=1.0

    for i=1:size(lf.lims, 1)
        v*=(lf.lims[i, 2]-lf.lims[i, 1])
    end

    v

end

"""
Check if a given position is part of a node

* `pt`: vector of floats with position
* `node`: `AbstractOctNode`
* `level`: integer, `>=0`. Indicates how many more levels should be fetched

* return: subdomain if succeeded, nothing if failed
"""
function get_subdomain(pt::Vector{Float64}, node::AbstractOctNode, level::Int64)

    if isa(node, OctLeaf)
        for ndim=1:size(node.lims, 1)
            if pt[ndim]<node.lims[ndim, 1] || pt[ndim]>node.lims[ndim, 2]
                return nothing
            end
        end

        return node
    end

    if isa(node, OctNode)
        if level==0
            for ndim=1:size(node.lims, 1)
                if pt[ndim]<node.lims[ndim, 1] || pt[ndim]>node.lims[ndim, 2]
                    return nothing
                end
            end

            return node
        else
            for sd in node.subnodes
                ret=get_subdomain(pt, sd, level-1)

                if !isnothing(ret)
                    return ret
                end
            end

            fth=node.father

            if isnothing(fth)
                return nothing
            else
                return get_subdomain(pt, fth, level+1)
            end
        end
    end

end

"""
Obtain subdomain within which a point is inserted, using a recursive search reaching 
any level of discretization

* `pt`: vector of floats with position
* `node`: `AbstractOctNode` from which to search

* return: `nothing` if not found, `AbstractOctNode` for the subdomain if found
"""
get_subdomain(pt::Vector{Float64}, node::AbstractOctNode) = get_subdomain(pt, node, typemax(Int64))

"""
Set adjacent attribute in `AbstractOctNode`
"""
function set_adjacent!(node::AbstractOctNode; recursive::Bool=false)

    fth=node.father

    if isnothing(fth)
        node.adjacent.=nothing

        return
    end

    ndim=size(node.lims, 1)
    coords=zeros(Float64, ndim)
    
    for cinds in CartesianIndices(node.adjacent)
        nd=cinds[1]
        naux=cinds[2]

        micro_offset=(naux==1 ? -1e-10 : 1e-10)

        for c=1:ndim
            if c!=nd
                coords[c]=(node.lims[c, 1]+node.lims[c, 2])/2
            else
                coords[c]=node.lims[c, naux]+micro_offset
            end
        end

        node.adjacent[cinds]=get_subdomain(coords, fth, 1)
    end

    if recursive && isa(node, OctNode)
        for sd in node.subnodes
            set_adjacent!(sd; recursive=true)
        end
    end

end

"""
Fragment a leaf into an `OctNode` and return the result

* `lf`: `OctLeaf` struct
"""
function fracture(lf::OctLeaf)

    nd=size(lf.lims, 1)

    original_halfdims=[(lf.lims[i, 2]-lf.lims[i, 1])/2 for i=1:nd]

    subdoms=Array{AbstractOctNode, nd}(undef, [2 for i=1:nd]...)

    newdom=OctNode{nd}(
        subdoms,
        lf.father,
        lf.lims,
        lf.adjacent
    )

    for cind in CartesianIndices(subdoms)
        lims=copy(lf.lims)

        for i=1:nd
            lims[i, 1]=(cind[i]==1 ? lf.lims[i, 1] : lf.lims[i, 1]+original_halfdims[i])
            lims[i, 2]=(cind[i]==2 ? lf.lims[i, 2] : lf.lims[i, 1]+original_halfdims[i])
        end

        newdom.subnodes[cind]=OctLeaf(
            lf.ind,
            newdom,
            lims,
            Matrix{Union{AbstractOctNode, Nothing}}(undef, nd, 2),
        )
    end

    for sd in newdom.subnodes
        set_adjacent!(sd)
    end

    newdom

end

"""
Function to fracture a cell and re-adjust the connectivity map

* `oct`: `AbstractOctNode` for an `OctNode`, father to the desired `OctLeaf`,
    indexed by `inds`. If `inds` is absent, `oct` is expected to be a leaf itself,
    and no update of the nearby cells' adjacency maps is done
* `inds`: indices to a leaf in `oct`. If absent, `oct` is assumed to be a leaf itself, but the
    adjacency maps aren't updated for nearby cells

* return: a new `OctNode` to replace `oct`
"""
function fracture!(oct::AbstractOctNode; inds::Union{Nothing, Vector{Int64}}=nothing)

    if isnothing(inds) && !isa(oct, OctLeaf)
        throw(
            error(
                "fracture!:index argument not provided, but octnode requested for fracture is not a leaf"
            )
        )
    end

    if !isnothing(inds) && isa(oct, OctLeaf)
        throw(
            error(
                "fracture!:index argument provided, but oct is a leaf"
            )
        )
    end

    if isnothing(inds)
        return fracture(oct)
    end

    oct.subnodes[inds...]=fracture(oct.subnodes[inds...])

    for adj in oct.subnodes[inds...].adjacent
        if !isnothing(adj)
            set_adjacent!(adj; recursive=true)
        end
    end

    oct

end

"""
Merge the leaves in a branch

* `nd`: `OctNode` struct for the branch

* return: an `OctLeaf` struct for the resulting branch
"""
function merge_branch(branch::OctNode{N}) where N

    inds=[
        begin
            @assert isa(sd, OctLeaf) "merge_branch:ERROR:attempting to merge a branch with more than one downstream ramification"

            sd.ind
        end for sd in branch.subnodes
    ]

    lf=OctLeaf(
        inds,
        branch.father,
        branch.lims,
        branch.adjacent
    )

    lf

end

"""
Function to merge a branch and re-adjust the connectivity map

* `oct`: `AbstractOctNode` for an `OctNode`, father to the desired branch,
    indexed by `inds`. If `inds` is absent, `oct` is expected to be the branch itself,
    and no update of the nearby cells' adjacency maps is done
* `inds`: indices to a branch in `oct`. If absent, `oct` is assumed to be the branch itself, but the
    adjacency maps aren't updated for nearby cells

* return: a new `OctNode` to replace `oct`
"""
function merge_branch!(oct::OctNode{N}; inds::Union{Nothing, AbstractVector}=nothing) where N

    if isa(inds, AbstractVector)
        @assert N==length(inds) "merge_branch!:DEBUG:wrong dimensionality passed for indexing vector"
    end

    if isnothing(inds)
        return merge_branch(oct)
    end

    oct.subnodes[inds...]=merge_branch(oct.subnodes[inds...])
    
    for adj in oct.subnodes[inds...].adjacent
        if !isnothing(adj)
            set_adjacent!(adj; recursive=true)
        end
    end

    oct

end

function _rec_nleaves(t::AbstractOctNode, alrdy::Int64)

    if isa(t, OctLeaf)
        return alrdy+1
    end

    n=alrdy

    for s in t.subnodes
        n=_rec_nleaves(s, n)
    end

    n

end

"""
Get number of leaves in octree (equivalent to the number of degrees of freedom for each 
variable in the solved differential equation)
"""
nleaves(t::AbstractOctNode) = _rec_nleaves(t, 0)

function _rec_setinds!(n, oldvs, newvs, upto)

    if isa(n, OctLeaf)
        upto+=1

        if isa(n.ind, AbstractArray)
            for (o, nw) in zip(oldvs, newvs)
                nw[upto]=sum(o[n.ind])/length(n.ind)
            end
        else
            for (o, nw) in zip(oldvs, newvs)
                nw[upto]=o[n.ind]
            end
        end

        n.ind=upto
    else
        for sd in n.subnodes
            upto=_rec_setinds!(sd, oldvs, newvs, upto)
        end
    end

    upto

end

"""
Function to set indices at leaves of the octree

* `msh`: `AbstractOctNode` corresponding to the first layer of resolution
* `oldvs`: arrays with values for variables as previously stocked 
"""
function set_indices!(msh::AbstractOctNode, oldvs::AbstractVector...)

    nlvs=nleaves(msh)

    newvs=[zeros(Float64, nlvs) for i=1:length(oldvs)]

    ret=_rec_setinds!(msh, collect(oldvs), newvs, 0)

    @assert ret==nlvs

    newvs

end

export OctLeaf, OctNode, fracture!, merge_branch!, nleaves, set_indices!, get_subdomain
