import Base:in

mutable struct StackedDiamondLattice
    generation
    stackingweight
    initial_state
    final_state
    edgeweightmap
end

function in(v::Integer, E::Edge)
    if v == E.src || v == E.dst
        return true
    else
        return false
    end
end

function weight(f::T, e)::Number where T <: AbstractMetaGraph
    w = getkey(f.eprops, e, f.defaultweight)
    if w isa Edge
        return f.eprops[w][f.weightfield]
    else
        return w
    end
end

function isstacked(f::T) where T <: AbstractMetaGraph
    vdicts = f.vprops |> values |> collect
    if isempty(filter(x -> :level in x, keys.(vdicts)))
        return false
    elseif allequal((x -> x[:level]).(vdicts))
        return false
    else
        return true
    end
end

function level_length(f::T) where T <: AbstractMetaGraph
    if !isstacked(f)
        return length(vertices(f))
    else
        vdicts = f.vprops |> values |> collect
        return filter(x -> x[:level] == 1, vdicts) |> length
    end
end

function getmaxlevel(f::T) where T <: AbstractMetaGraph
    if !isstacked(f)
        return 1
    else
        vlevels = (x -> x[:level]).(f.vprops |> values |> collect)
        return maximum(vlevels)
    end
end

function diamondstack(f::T, g::T; stackingweight) where T <: AbstractMetaGraph
    vf = collect(vertices(f))
    vg = collect(vertices(g))

    # Assert that lattices must be unstacked
    @assert !isstacked(f) "Please use diamondconcat for concatinating two stacked lattices"
    @assert !isstacked(g) "Please use diamondconcat for concatinating two stacked lattices"

    offset = length(vf)
    h = T(length(vf) + length(vg))

    # Assign levels
    for v in vf
        set_prop!(h, v, :level, 1)
    end
    for v in vg
        set_prop!(h, offset + v, :level, 2)
    end
    
    # Copy over edges
    for e in edges(f)
        add_edge!(h, e.src, e.dst)
        set_prop!(h, e, h.weightfield, 1.)
    end
    for e in edges(g)
        add_edge!(h, e.src+offset, e.dst+offset)
        set_prop!(h, Edge(e.src+offset, e.dst+offset), h.weightfield, 1.)
    end
    
    # Connect two lattices
    for v in vf[end-length(vg)+1:end]
        add_edge!(h, v, offset+v)
        set_prop!(h, Edge(v, offset+v), h.weightfield, stackingweight)
    end

    # Copy over metadata
    for v in vf
        set_props!(h, v, f.vprops[v])
    end
    for v in vg
        set_props!(h, offset + v, g.vprops[v])
    end

    return h
end

function diamondconcat(f::T, g::T; stackingweight) where T <: AbstractMetaGraph
    vf = collect(vertices(f))
    vg = collect(vertices(g))
    
    @assert isstacked(f)
    @assert level_length(f) == level_length(g)

    offset = length(vf)
    fdepth = getmaxlevel(f)
    level_len = level_length(f)
    h = T(length(vf) + length(vg))

    # Copy metadata
    for v in vf
        set_props!(h, v, f.vprops[v])
    end
    for v in vg
        set_props!(h, offset + v, g.vprops[v])
    end
    
    # Set levels
    for v in vf
        set_prop!(h, v, :level, f.vprops[v][:level])
    end
    if isstacked(g)
        for v in vg
            set_prop!(h, offset + v, :level, fdepth + g.vprops[v][:level])
        end
    else
        for v in vg
            set_prop!(h, offset + v, :level, fdepth + 1)
        end
    end

    # Copy over edges
    for e in edges(f)
        add_edge!(h, e.src, e.dst)
        set_prop!(h, e, h.weightfield, weight(f, e))
    end
    for e in edges(g)
        add_edge!(h, e.src+offset, e.dst+offset)
        set_prop!(h, Edge(e.src+offset, e.dst+offset), h.weightfield, weight(g, e))
    end

    # Connect lattices
    for v in vf[end-level_len:end]
        add_edge!(h, v, v + level_len)
        set_prop!(h, Edge(v, v + level_len), h.weightfield, stackingweight)
    end

    return h
end

function enumerate_edges(f::T, v) where T <: AbstractGraph
    edgelist = Edge[]
    for e in edges(f)
        if v in e
            push!(edgelist, e)
        end
    end
    
    return edgelist
end

function enumerate_edgeweights(f::T, v) where T <: AbstractMetaGraph
    connected_edges = enumerate_edges(f, v)
    weightmap = Dict{Float64, Vector{Int64}}()
    
    for e in connected_edges
        if v == e.src
            si = e.dst
        else
            si = e.src
        end
        
        w = weight(f, e)
        if !(w in keys(weightmap))
            push!(weightmap, w => [si])
        else
            push!(weightmap[w], si)
        end
    end
    
    return weightmap
end

function stack(lattices::Vector{T}, stackingweight) where T <: AbstractMetaGraph
    @assert length(lattices) > 1 "Need more than 1 lattice to stack!"
    
    l = diamondstack(lattices[1], lattices[2]; stackingweight = stackingweight)
    for nextlattice in lattices[3:end]
        l = diamondconcat(l, nextlattice; stackingweight = stackingweight)
    end
    
    return l
end

function StackedDiamondLattice(order, depth, stackingweight, initstate)
    lattice = stack(
        [diamond_ising_lattice(order, initstate) for _ in Base.OneTo(depth)],
        stackingweight
        )
    edgeweightmap = Dict(vertices(lattice) .=> (x -> enumerate_edgeweights(lattice, x)).(vertices(lattice)))
    
    return StackedDiamondLattice(order, stackingweight, deepcopy(lattice), deepcopy(lattice), edgeweightmap)
end

function energy(l::StackedDiamondLattice; state = :final_state)
    f = getproperty(l, state)

    total = 0
    for e in edges(f)
        if e in keys(f.eprops)
            weight = f.eprops[e][:weight]
        else
            weight = f.defaultweight
        end
        total += -weight * f.vprops[e.src][:val]*f.vprops[e.dst][:val]
    end
    return total
end

function nnsum(f::T, vlist) where T <: AbstractMetaGraph
    total = 0
    for v in vlist
        total += f.vprops[v][:val]
    end
    return total
end

function ΔE(l::StackedDiamondLattice, sk; state = :final_state)
    f = getproperty(l, state)
    total = 0
    spink = f.vprops[sk][:val]
    for (w, nnspins) in l.edgeweightmap[sk]
        total += 2*w*spink*nnsum(f, nnspins)
    end
    return total
end

function magnetization(l::StackedDiamondLattice; state = :final_state)
    f = getproperty(l, state)
    sum([f.vprops[v][:val] for v in vertices(f)])
end

function _fill_M_history!(::StackedDiamondLattice, data; showprogress = false)
    lattice = deepcopy(data.lattice)
    lattice.final_state = lattice.initial_state

    P = Progress(length(data.spinflip_history), desc = "Filling Magnetization History...")

    data.magnetization_history = Float64[magnetization(lattice)]
    for s_k in data.spinflip_history
        if s_k == -1
            push!(data.magnetization_history, data.magnetization_history[end])
        else
            # Update lattice for next calculation
            lattice.final_state.vprops[s_k][:val] *= -1

            # Calculate change in magnetization
            push!(data.magnetization_history, data.magnetization_history[end] + 2*lattice.final_state.vprops[s_k][:val])
        end

        if showprogress
            next!(P)
        end
    end
end

function _fill_U_history!(::StackedDiamondLattice, data; showprogress = false)
    lattice = deepcopy(data.lattice)
    lattice.final_state = lattice.initial_state

    P = Progress(length(data.spinflip_history), desc = "Filling Internal Energy History...")

    data.internalenergy_history = Float64[energy(lattice)]
    for s_k in data.spinflip_history
        if s_k == -1
            push!(data.internalenergy_history, data.internalenergy_history[end])
        else
            push!(data.internalenergy_history, data.internalenergy_history[end] + ΔE(lattice, s_k))
            lattice.final_state.vprops[s_k][:val] *= -1
        end
        if showprogress
        next!(P)
    end
    end
end