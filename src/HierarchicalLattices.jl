module HierarchicalLattices

using Graphs
using MetaGraphs
using RecipesBase
using ProgressMeter

# Diamond Lattices
include("DiamondLattices/Lattice.jl")
include("DiamondLattices/StackedLattice.jl")
include("DiamondLattices/Metropolis.jl")
include("DiamondLattices/Cluster.jl")

# Diamond Lattice data
export  diamond_order_zero_transform!,
    make_diamond_lattice0,
    diamond_lattice,
    diamond_ising_lattice,
    numberofspins

# MCMC Stuff
export  metropolis!, energy, ΔE, magnetization, fill_data!
export wolffstep!, wolff!

# Main data structures
export IsingData, WolffData, DiamondLattice, StackedDiamondLattice

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

end # module HierarchicalLattices