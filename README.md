# OctFVM.jl

OctFVM.jl is a package built to assist in the assembly of finite volume solvers using octree meshes (and their equivalents with an arbitrary number of dimensions).

## Building Meshes

To start meshing with a 1 x 1 rectangle in a two-dimensional space:

```
msh = OctLeaf(
    [
        0.0 1.0;
        0.0 1.0
    ] # limits in each dimension
)
```

Alternatively, in three dimensions:

```
msh = OctLeaf(
    [
        0.0 1.0;
        0.0 1.0;
        1.0 2.0 # or any other definition for the z axis
    ] # limits in each dimension
)
```

To split the current cell into four (or any N-dimensional equivalent) leaves:

```
msh = fracture!(msh)

isa(msh, OctNode{2})
# true

nleaves(msh)
# 4
```

Or, to do it again with a leaf from the recently split branch:

```
msh = fracture!( msh; inds = [1, 2] )

nleaves(msh)
# 7

# even further!
msh.subnodes[1, 2] = fracture!( msh.subnodes[1, 2]; inds=[2, 1] )
```

Notice that only the second execution (with a specification for a father node and an index for the desired leaf) reconfigures the adjacency maps of nearby cells, and must be used when refining a specific leaf within a complex tree.

## Attributing Data to Leaves

To create vectors relating data to specific leaves:

```
data = zeros(Float64, nleaves(msh))
set_indices!(msh)
```

`set_indices!` ensures each leaf has a specific index to relate it to vectors storing data. It can also be used to reshape vectors accordingly when refining or coarsening a mesh:

```
data2 = similar(data)

fracture!(msh; inds = [2, 2])
data, data2 = set_indices!(msh, data, data2) # inherits data from fractured cells

merge_branch!(msh; inds = [1, 2])
data, data2 = set_indices!(msh, data, data2) # averages data between merged cells
```

## Calculating Fluxes

Suppose `f` is a vector with y-axis directed advection fluxes. x-axis directed fluxes of value `f_bc` are also applied in the aft, x-axis face of the mesh.
To calculate the expected variations in a volume-averaged variable of interest `u` for a forward Euler scheme, one could use:

```
# obtaining sparse matrices for linear transformations:
y_flux_lintransf_rear = neighbor_lintransf(msh, 2, 1) # respectively, a dimension and an index for the face.

# (2, 1), for example, stands for the first (rear) face in the second dimension

y_flux_lintransf_aft = neighbor_lintransf(msh, 2, 2)

# sparse vector for BC flux:
x_flux_lintransf = get_BC_flux(msh, 1, 2)

flux_integrals = (y_flux_rear * f .- y_flux_aft * f) .- (x_flux_lintransf .* f_bc)

volumes = volume_vector(msh)

u .+= dt .* (flux_integrals ./ volumes)
```

One could also "manually" obtain vectors of interface areas with leaves of neighboring regions:

```
lfs, areas = get_neighbor_interfaces(msh.subnodes[1, 2], 1, 1) # neighbors to the rear face of the first spatial dimension

neighbor_data = [data[lf.ind] for lf in lfs]
```

Or obtain cells within meshes or nodes that have contact with neighboring regions, and their respective interface areas:

```
lfs, areas = get_frontier_interfaces(msh, 1, 1) # "first" leaves when marching from the first (rear) face in the first dimension
```

## Further Possibilities

Get father for a node or leaf:

```
nd.father
```

Get adjacent domain of equal or higher refinement level:

```
nd.adjacent[1, 2] # neighbor to aft face in first dimension
# may point to a leaf or a node, depending on its refinement level
```

For further information, you can access the docstrings for each object and data type in Julia's help mode.

## Installation

You can install the package from its GitHub repository:

```
]add https://github.com/pedrosecchi67/OctFVM.jl
```
