module HierarchicalLattices

using Graphs
using MetaGraphs
using RecipesBase
using ProgressMeter
include("DiamondLattices/Lattice.jl")
include("DiamondLattices/Metropolis.jl")


# Diamond Lattice data
export  diamond_order_zero_transform!, make_diamond_lattice0, diamond_lattice, diamond_ising_lattice

# Metropolis Stuff
export  metropolis!, energy, ΔE, magnetization, fill_data!

# Main data structure
export IsingData, DiamondLattice

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

include("DiamondLattices/Lattice.jl")

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

function _fill_U_history!(data::IsingData; J = 1, showprogress = false)
    lattice = data.lattice
    l = deepcopy(lattice.initial_state)
    P = Progress(length(data.spinflip_history), desc = "Filling Internal Energy History...")
    data.internalenergy_history = Float64[energy(l)]
    for s_k in data.spinflip_history
        if s_k == -1
            push!(data.internalenergy_history, data.internalenergy_history[end])
        else
            # Calculate new E
            push!(data.internalenergy_history, data.internalenergy_history[end] + ΔE(l, s_k, neighbors(l, s_k), J = J))
            
            # Update spin for next calculation
            l.vprops[s_k][:val] *= -1
        end

        if showprogress
            next!(P)
        end
    end
end

function _fill_M_history!(data::IsingData; showprogress = false)
    lattice = data.lattice
    l = deepcopy(lattice.initial_state)
    P = Progress(length(data.spinflip_history), desc = "Filling Magnetization History...")
    data.magnetization_history = Float64[magnetization(l)]
    for s_k in data.spinflip_history
        if s_k == -1
            push!(data.magnetization_history, data.magnetization_history[end])
        else
            # Update lattice for next calculation
            l.vprops[s_k][:val] *= -1
            
            # Calculate change in magnetization then
            push!(data.magnetization_history, data.magnetization_history[end] + 2*l.vprops[s_k][:val])
        end

        if showprogress
            next!(P)
        end
    end
end

"""
    Fill the internal energy or magnetization history of a lattice evolved using the
    `HierarchicalLattices.metropolis!` function.
"""
function fill_data!(lattice, data::Symbol; showprogress = false)
    if data == :M
        return _fill_M_history!(lattice, showprogress = showprogress)
    elseif data == :U
        return _fill_U_history!(lattice, showprogress = showprogress)
    else
        throw(ArgumentError("Data parameter $(string(data)) not implimented!"))
    end
end

end # module HierarchicalLattices