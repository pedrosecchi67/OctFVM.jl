function _rec_get_fr(oct, lfs, areas, ndim, n_aux)

    if isa(oct, OctLeaf)
        push!(lfs, oct)
        push!(areas, leaf_area(oct, ndim))
    else
        for cinds in CartesianIndices(oct.subnodes)
            if cinds[ndim]==n_aux
                _rec_get_fr(oct.subnodes[cinds], lfs, areas,  ndim, n_aux)
            end
        end
    end

end

"""
Function to get vectors for all leaves in the frontier of an `AbstractOctNode`

* `oct`: `AbstractOctNode` at hand
* `ndim`: number of the searched dimension
* `n_aux`: 1 if looking for rear interfaces, 2 if otherwise

* return: a vector to all interfaced leaves and a vector with their respective interface area
"""
function get_frontier_interfaces(oct::AbstractOctNode, ndim::Int64, n_aux::Int64)

    lfs=Vector{AbstractOctNode}()
    areas=Vector{Float64}()

    _rec_get_fr(oct, lfs, areas, ndim, n_aux)

    return lfs, areas

end

"""
Function to get interfaces from neighboring leaves

* `oct`: `AbstractOctNode` at hand
* `ndim`: number of the searched dimension
* `n_aux`: 1 if looking for rear interfaces, 2 if otherwise

* return: a vector to all interfaced leaves and a vector with their respective interface area
"""
function get_neighbor_interfaces(oct::AbstractOctNode, ndim::Int64, n_aux::Int64)

    if !isnothing(oct.adjacent[ndim, n_aux])
        lfs, areas=get_frontier_interfaces(
            oct.adjacent[ndim, n_aux],
            ndim,
            (n_aux==1 ? 2 : 1)
        )

        if length(areas)==1
            areas[1]=leaf_area(oct, ndim)
        end

        return lfs, areas
    end
    
    lfs=Vector{AbstractOctNode}()
    areas=Vector{Float64}()

    return lfs, areas

end

function _rec_maxind(oct, prev)

    if isa(oct, OctLeaf)
        return (prev > oct.ind ? prev : oct.ind)
    end

    mx=prev

    for sd in oct.subnodes
        mx=_rec_maxind(sd, mx)
    end

    mx

end

function _rec_setlvint(oct, ndim, n_aux, m)

    if isa(oct, OctLeaf)
        lfs, areas=get_neighbor_interfaces(oct, ndim, n_aux)

        for (lf, a) in zip(lfs, areas)
            m[oct.ind, lf.ind]+=a
        end
    else
        for sd in oct.subnodes
            _rec_setlvint(sd, ndim, n_aux, m)
        end
    end

end

"""
Function to get interfaces to get rear or aft interfaces to leaves in the form of 
a linear transformation.

**Pressuposes that `set_indexing!` has already been ran**

* `oct`: `AbstractOctNode` at hand
* `ndim`: number for the dimension for which the linear transformation is desired
* `n_aux`: 1 if expecting rear interfaces, 2 if otherwise
* `nlvs`: number of leaves, for matrix order. Surmised from leaves in `oct` if not provided

* return: a sparse matrix mapping values to `∫u_adj dS` surface integrals
"""
function neighbor_lintransf(oct::AbstractOctNode, ndim::Int64, n_aux::Int64;
    nlvs::Union{Int64, Nothing}=nothing)

    _nlvs=(isnothing(nlvs) ? _rec_maxind(oct, 0) : nlvs)

    m=spzeros(_nlvs, _nlvs)

    _rec_setlvint(oct, ndim, n_aux, m)

    m

end

function _rec_volset(oct, v)

    if isa(oct, OctLeaf)
        v[oct.ind]=leaf_volume(oct)
    else
        for sd in oct.subnodes
            _rec_volset(sd, v)
        end
    end

end

"""
Function to generate a vector of cell volumes through recursion through the tree

* `oct`: `AbstractOctNode` at hand
* `nlvs`: number of leaves, for vector order. Surmised from leaves in `oct` if not provided
"""
function volume_vector(oct::AbstractOctNode; nlvs::Union{Int64, Nothing}=nothing)

    _nlvs=(isnothing(nlvs) ? _rec_maxind(oct, 0) : nlvs)

    v=zeros(Float64, _nlvs)

    _rec_volset(oct, v)

    v

end

"""
Get a vector that implies on the additional fluxes to all cells in a given direction
when a unit Neumann condition is applied at a frontier

* `oct`: `AbstractOctNode` at hand
* `ndim`: number for the dimension for which the flux is desired
* `n_aux`: 1 if expecting rear interfaces, 2 if otherwise
* `nlvs`: number of leaves, for vector order. Surmised from leaves in `oct` if not provided
"""
function get_BC_flux(oct::AbstractOctNode, ndim::Int64, n_aux::Int64;
    nlvs::Union{Int64, Nothing}=nothing)

    _nlvs=(isnothing(nlvs) ? _rec_maxind(oct, 0) : nlvs)

    v=spzeros(_nlvs)

    lfs, areas=get_frontier_interfaces(oct, ndim, n_aux)

    for (lf, a) in zip(lfs, areas)
        v[lf.ind]=a
    end

    v

end

export get_frontier_interfaces, get_neighbor_interfaces, neighbor_lintransf, volume_vector, get_BC_flux
