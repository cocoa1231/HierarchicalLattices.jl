mutable struct DiamondLattice
    generation
    initial_state
    final_state
end

DiamondLattice(G::MetaGraph, order::Integer) = DiamondLattice(order, deepcopy(G), deepcopy(G))

"""
    Takes a lattice graph and an edge and applies b = 2 diamond transform
"""
function diamond_order_zero_transform!(lattice, edge; scale = 1)
    
    e = edge
    src = lattice.vprops[e.src][:loc]
    dst = lattice.vprops[e.dst][:loc]

    Δz = dst - src
    new_left  = (Δz * 1im) *  scale
    new_right = (Δz * -1im) * scale
    disp = src + Δz/2
    new_pts = [new_left, new_right]
    add_vertices!(lattice, 2)
    for (idx, new_addition) in enumerate(vertices(lattice)[end-1:end])
        set_props!(lattice, new_addition, Dict(
                :loc => new_pts[idx] + disp,
                :val => false
            ))
        add_edge!(lattice, e.src, new_addition)
        add_edge!(lattice, new_addition, e.dst)
    end

    rem_edge!(lattice, e) 

end

"""
    Creates an order 0, b = 2 diamond lattice
"""
function make_diamond_lattice0()
    order_zero = SimpleDiGraph(2, 1) |> MetaGraph
    locations = [ 0. + 1im, 0. - 1im ]
    values = [false, false]
    for (idx, vertex) in enumerate(vertices(order_zero))

        set_props!(order_zero, vertex, Dict(
                :loc => locations[idx],
                :val => values[idx]
            ))
    end
    return order_zero
end

"""
    Generates an arbitrary order diamond lattice
"""
function diamond_lattice(order::Integer)
    oz = make_diamond_lattice0()
    
    for i in 1:order
        for e in collect(edges(oz))
            diamond_order_zero_transform!(oz, e; scale = 0.5^i)
        end
    end
    
    return oz
end

"""
    Generates an arbitrary order diamond Ising lattice given an initial state
    where the initial state can be :zero or :infty
"""
function diamond_ising_lattice(order::Integer, state::Symbol)
    l = diamond_lattice(order)
    if state == :zero
        for v in vertices(l)
            set_prop!(l, v, :val, 1)
        end
    elseif state == :infty
        for v in vertices(l)
            set_prop!(l, v, :val, rand([+1, -1]))
        end
    end
    return l
end