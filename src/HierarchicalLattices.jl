module HierarchicalLattices

using Graphs
using MetaGraphs
using RecipesBase
using ProgressMeter
include("Metropolis.jl")


# Diamond Lattice data
export  diamond_order_zero_transform!, make_diamond_lattice0, diamond_lattice, diamond_ising_lattice

# Metropolis Stuff
export  metropolis!, energy, ΔE, magnetization, fill_data!

# Main data structure
export IsingData

@recipe function f(L::MetaGraph)
    aspect_ratio := :equal
    
    for e in edges(L)
        @series begin
            seriestype := :path
            label := false
            linecolor --> :green
            src_x, src_y = real(L.vprops[e.src][:loc]), imag(L.vprops[e.src][:loc])
            dst_x, dst_y = real(L.vprops[e.dst][:loc]), imag(L.vprops[e.dst][:loc])
            [src_x, dst_x], [src_y, dst_y]
        end
    end

    @series begin
        vx, vy = Float64[], Float64[]
        seriestype := :scatter
        label := false
        color --> :blue
        for v in vertices(L)
            if L.vprops[v][:val] == +1
                push!(vx, real(L.vprops[v][:loc]))
                push!(vy, imag(L.vprops[v][:loc]))
            end
        end
        vx, vy
    end

    @series begin
        vx, vy = Float64[], Float64[]
        seriestype := :scatter
        label := false
        color --> :black
        for v in vertices(L)
            if L.vprops[v][:val] == -1
                push!(vx, real(L.vprops[v][:loc]))
                push!(vy, imag(L.vprops[v][:loc]))
            end
        end
        vx, vy
    end
    
end

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

function diamond_lattice(order::Integer)
    oz = make_diamond_lattice0()
    
    for i in 1:order
        for e in collect(edges(oz))
            diamond_order_zero_transform!(oz, e; scale = 0.5^i)
        end
    end
    
    return oz
end

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

function sum_autocorr(magnetization_history, t)
    # t_max is the number of monte carlo steps we've taken in total
    t_max = length(magnetization_history)
    
    # Magnetization history m, assuming I've run fill_magnetization_history! already on the lattice
    m = magnetization_history
    
    # Normalization constant
    norm = 1/(t_max - t)
    
    mean_m² = norm * sum(m[1:t_max-t]) * norm * sum(m[t+1:t_max])
    
    return norm * sum(m[1:t_max - t] .* m[t+1:t_max]) - mean_m²
end

function generate_autocorr_data(array, N, nsweeps; showprogress = false)
    m = array
    t_max = floor(Int, nsweeps*N)

    χ = zeros( ceil(Integer, t_max / N) )

    if showprogress
        P = Progress(length(χ))
    end
    Threads.@threads for t in eachindex(χ)
        χ[t] = sum_autocorr(m, Int(t*N))
        if showprogress
            next!(P)
        end
    end

    return χ
end

function _fill_U_history!(lattice::IsingData; J = 1)
    l = deepcopy(lattice.initial_state)
    
    lattice.internalenergy_history = Float64[energy(l)]
    for s_k in lattice.spinflip_history
        if s_k == -1
            push!(lattice.internalenergy_history, lattice.internalenergy_history[end])
        else
            # Calculate new E
            push!(lattice.internalenergy_history, lattice.internalenergy_history[end] + ΔE(l, s_k, neighbors(l, s_k), J = J))
            
            # Update spin for next calculation
            l.vprops[s_k][:val] *= -1
        end
    end
end

function _fill_M_history!(lattice::IsingData)
    l = deepcopy(lattice.initial_state)
    
    lattice.magnetization_history = Float64[magnetization(l)]
    for s_k in lattice.spinflip_history
        if s_k == -1
            push!(lattice.magnetization_history, lattice.magnetization_history[end])
        else
            # Update lattice for next calculation
            l.vprops[s_k][:val] *= -1
            
            # Calculate change in magnetization then
            push!(lattice.magnetization_history, lattice.magnetization_history[end] + 2*l.vprops[s_k][:val])
        end
    end
end

"""
    Fill the internal energy or magnetization history of a lattice evolved using the
    `HierarchicalLattices.metropolis!` function.
"""
function fill_data!(lattice, data::Symbol)
    if data == :M
        return _fill_M_history!(lattice)
    elseif data == :U
        return _fill_U_history!(lattice)
    else
        throw(ArgumentError("Data parameter $(string(data)) not implimented!"))
    end
end

end # module HierarchicalLattices