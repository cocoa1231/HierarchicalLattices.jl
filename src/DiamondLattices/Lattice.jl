"""
Main type for a diamond lattice. Takes four fields
    - `generation`: The order or generation of the lattice.
    - `initial_state`: The initial state of the lattice at time of creation.
    - `final_state`: Final state of the lattice after MCMC simulation.
    - `interaction_weight`: Strength of the bonds of the lattice.

No type restrictions are imposed on the fields (expect to be changed). There are
two external constructors.
    - `DiamondLattice(G::MetaGraph, o::Integer)`: Construct a diamond lattice
    with initiial and final states `G` and `o`. The strength of the lattice will
    be `G.defaultweight`. Note: this will create a deepcopy of `G`.
    - `DiamondLattice(G::MetaGraph, o::Integer, latticeweight::Number)`: 
    Construct a diamond lattice with initiial and final states `G` and `o`. The
    strength of the lattice will be `latticeweight`. Note: this will create a
    deepcopy of `G`.
"""
mutable struct DiamondLattice
    generation
    initial_state
    final_state
    interaction_weight
end

DiamondLattice(G::MetaGraph, order::Integer) = DiamondLattice(order, deepcopy(G), deepcopy(G), G.defaultweight)
DiamondLattice(G::MetaGraph, order::Integer, latticeweight::Number) = DiamondLattice(order, deepcopy(G), deepcopy(G), latticeweight)

"""
Takes a lattice graph and an edge and applies b = 2 diamond transform. Thus
function will also calculate the physical locations of the new spins with
respect to locations of the old one, and thus should only be used for visual
demonstration purposes. It is, otherwise, quite memory inefficient. Inputs for
the function are:
    - lattice: The `MetaGraph` to apply the transformation to.
    - edge: The edge to apply the transformation to.
    - scale (default = 1): How much to scale the locations of the new bonds by.
    A scale of 1 indicates that the distance between the new bond and `edge`
    should be equal to the length of `edge`.
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
Creates an order 0, b = 2 diamond lattice with spin values set to `false`
and locations of the spins at `1im` and `-1im`. Returns a MetaGraph.
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
Generates an arbitrary order diamond lattice with b = 2. Thus function will
calculate the location of every spin. The location of every new spin will be
scaled by $0.5^o$ where $o$ is the order of the spin. Inputs:
    - order::Integer: Order of the lattice to generate
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
where the initial state can be :zero or :infty. This method will calculate the
location of every spin as well for plotting. Inputs:
    - `order::Integer`: Order of lattice to generate.
    - `state::Symbol`: Can either be `:zero` for all spins pointing up (+1) or
    `:infty` for all spins pointing randomly.
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

"""
Takes a lattice graph and an edge and applies diamond transform. This will apply
the diamond transform with arbitrary `b`. Inputs:
    - `lattice`: `MetaGraph` lattice to apply the transformation to.
    - `edge`: Edge to apply the transformation to.
    - `b`: Number of new branches created.
"""
function transform_edge_diamond!(lattice, edge, b)
    
    e = edge
    add_vertices!(lattice, b)
    for new_addition in vertices(lattice)[end-b+1:end]
        add_edge!(lattice, e.src, new_addition)
        add_edge!(lattice, new_addition, e.dst)
    end
    rem_edge!(lattice, e) 

end

"""
Creates an order 0 diamond lattice. Does not calculate the location of the
spins.
"""
function order_zero_diamond_lattice()
    order_zero = SimpleGraph(2, 1) |> MetaGraph
    return order_zero
end

"""
Raise order of diamond lattice. This will apply `transform_edge_diamond!` to all
edges of the lattice, given a `b`. Inputs:
    - `lattice`: `MetaGraph` to apply the transformation to.
    - `b`: Number of new branches created for every edge.
    - `showprogress` (default `false`): Add a progress bar.
    - `gc_freq` (default `0.01`): Since the function iterates over every edge
    and adds new nodes with temporary variables at each iteration, for larger
    lattices this parameter can be helpful to keep memory usage down. Setting
    this to anything above `0` disables multithreading. At every iteration, run
    the garbage collector with the probability `gc_freq`. Any value above `1`
    will indicate running this at every iteration.
"""
function raise_order_diamond!(lattice, b; showprogress = false, gc_freq = 0.01)
    elist = collect(edges(lattice))
    P = Progress(length(elist), enabled = showprogress)
    for e in elist
        transform_edge_diamond!(lattice, e, b)
        next!(P)
        rand() < gc_freq && GC.gc()
    end
end

"""
Generates an arbitrary order diamond lattice. This function does not calculate
physical locations of every new point. Inputs:
    - `order::Int64`: Order of the lattice to create.
    - `b::Int64`: Branching number of the lattice.
    - `showprogress` (default `false`): Add a progress bar. 
    - `gc_freq` (default `0.`): Since the function iterates over every edge
    and adds new nodes with temporary variables at each iteration, for larger
    lattices this parameter can be helpful to keep memory usage down. Setting
    this to anything above `0` disables multithreading. At every iteration, run
    the garbage collector with the probability `gc_freq`. Any value above `1`
    will indicate running this at every iteration.
"""
function diamond_lattice(order::Int64, b::Int64; showprogress = false, gc_freq = 0.)
    oz = order_zero_diamond_lattice()
    
    for _ in 1:order
        raise_order_diamond!(oz, b; showprogress = showprogress, gc_freq = gc_freq)
    end 
    return oz
end

"""
Generates an arbitrary order diamond Ising lattice given an initial state
where the initial state can be :zero or :infty. This function does not calculate
the location of every spin. Inputs:
    - `order::Int64`: Order of the lattice to generate.
    - `b::Int64`: Branching number of the lattice.
    - `showprogress` (default `false`): Add a progress bar. 
    - `gc_freq` (default `0.`): Since the function iterates over every edge
    and adds new nodes with temporary variables at each iteration, for larger
    lattices this parameter can be helpful to keep memory usage down. Setting
    this to anything above `0` disables multithreading. At every iteration, run
    the garbage collector with the probability `gc_freq`. Any value above `1`
    will indicate running this at every iteration.
"""
function diamond_ising_lattice(order::Int64, b::Int64, state::Symbol; showprogress = false, gc_freq = 0.)
    l = diamond_lattice(order, b, showprogress = showprogress, gc_freq = gc_freq)
    states = [+1, -1]
    if state == :zero
        for v in vertices(l)
            set_prop!(l, v, :val, 1)
        end
    elseif state == :infty
        for v in vertices(l)
            set_prop!(l, v, :val, rand(states))
        end
    end
    return l
end

function magnetization(lattice::DiamondLattice; state = :final)
    L = getproperty(lattice, Symbol(string(state)*"_state"))
    M = sum([L.vprops[v][:val] for v in vertices(L)])
    return M
end

function energy(lattice::DiamondLattice; state = :final)
    f = getproperty(lattice, Symbol(string(state)*"_state"))
    
    return sum( (e -> -lattice.interaction_weight*get_prop(f, e.src, :val)*get_prop(f, e.dst, :val)).(edges(f)) )
end

function ΔE(lattice::DiamondLattice, s::Integer, n::Vector{<:Integer}; state = :final)
    L = getproperty(lattice, Symbol(string(state)*"_state"))
    J = lattice.interaction_weight
    si = [ L.vprops[i][:val] for i in n ]
    sk = L.vprops[s][:val]
    return 2J*sk*sum(si)
end

function _fill_U_history!(::DiamondLattice, data; showprogress = false, progressoutput = stdout)
    lattice = DiamondLattice(data.lattice.initial_state, data.lattice.generation)
    P = Progress(length(data.spinflip_history), desc = "Filling Internal Energy History...", output = progressoutput)

    data.internalenergy_history = Float64[energy(lattice)]
    for s_k in data.spinflip_history
        if s_k == -1
            push!(data.internalenergy_history, data.internalenergy_history[end])
        else
            # Calculate new E
            push!(data.internalenergy_history, data.internalenergy_history[end] + ΔE(lattice, s_k, neighbors(lattice.final_state, s_k)))
            
            # Update spin for next calculation
            lattice.final_state.vprops[s_k][:val] *= -1
        end

        if showprogress
            next!(P)
        end
    end
end

function _fill_M_history!(::DiamondLattice, data; showprogress = false, progressoutput = stdout)
    lattice = DiamondLattice(data.lattice.initial_state, data.lattice.generation)
    P = Progress(length(data.spinflip_history), desc = "Filling Magnetization History...", output = progressoutput)

    data.magnetization_history = Float64[magnetization(lattice)]
    for s_k in data.spinflip_history
        if s_k == -1
            push!(data.magnetization_history, data.magnetization_history[end])
        else
            # Update lattice for next calculation
            lattice.final_state.vprops[s_k][:val] *= -1
            
            # Calculate change in magnetization then
            push!(data.magnetization_history, data.magnetization_history[end] + 2*lattice.final_state.vprops[s_k][:val])
        end

        if showprogress
            next!(P)
        end
    end
end
