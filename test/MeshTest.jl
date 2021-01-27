@test begin
    msh=OctLeaf(
        [
            0.0 1.0;
            0.0 1.0
        ]
    )

    frc=fracture!(msh)

    @assert frc.subnodes[1, 1].adjacent[1, 2]==frc.subnodes[2, 1]
    @assert frc.subnodes[1, 1].adjacent[1, 2].adjacent[1, 1]==frc.subnodes[1, 1]

    set_indices!(frc)

    xs=[1.0, 2.0, 3.0, 4.0]
    ys=[4.0, 3.0, 2.0, 1.0]

    xs, ys=set_indices!(frc, xs, ys)

    @assert nleaves(frc)==4

    frc=fracture!(frc; inds=[1, 1])
    
    @assert nleaves(frc)==7

    xs, ys=set_indices!(frc, xs, ys)
    
    @assert length(xs)==7 && length(ys)==7

    @assert frc.subnodes[1, 1].subnodes[2, 1].adjacent[1, 2]==frc.subnodes[2, 1]
    @assert frc.subnodes[2, 1].adjacent[1, 1]==frc.subnodes[1, 1]

    lfs, areas=get_frontier_interfaces(frc, 1, 1)
    @assert isapprox(areas, [0.25, 0.25, 0.5])

    lfs, areas=get_neighbor_interfaces(frc.subnodes[1, 1].subnodes[2, 1], 1, 2)
    @assert lfs[1]==frc.subnodes[2, 1]
    @assert isapprox(areas, [0.25])

    y_neighbors=neighbor_lintransf(frc, 2, 2)

    volset=volume_vector(frc)
    @assert isapprox(volset, [0.0625, 0.0625, 0.0625, 0.0625, 0.25, 0.25, 0.25])

    bc_flux=get_BC_flux(frc, 1, 1)

    frc=merge_branch!(frc; inds=[1, 1])

    xs, ys=set_indices!(frc, xs, ys)

    @assert isapprox(xs, [1.0, 2.0, 3.0, 4.0])
    @assert isapprox(ys, [4.0, 3.0, 2.0, 1.0])
    
    true
end
