export IsingData
export magnetization
export energy

mutable struct IsingData
    generation::Integer
    initial_state::MetaGraph
    final_state::MetaGraph
    magnetization_history::Vector{Float64}
    internalenergy_history::Vector{Float64}
    spinflip_history::Vector{Int64}
end

IsingData(L::MetaGraph, g::Integer) = IsingData(g, deepcopy(L), deepcopy(L), Float64[], Float64[], Int64[])

function IsingData(L::MetaGraph)
    nbonds = edges(L) |> collect |> length
    g = Int(log(4, nbonds))
    return IsingData(L, g)
end

function magnetization(L::MetaGraph)
    N = vertices(L) |> length
    M = sum([L.vprops[v][:val] for v in vertices(L)])
    
    return M
end

function energy(L::MetaGraph; J = 1)
    E = 0
    for e in edges(L)
        E += L.vprops[e.src][:val]*L.vprops[e.dst][:val]
    end
    return -J*E
end
energy(L::IsingData; J = 1) = energy(L.final_state, J = J)

function ΔE(L::MetaGraph, s::Integer, n::Vector{<:Integer}; J = 1)
    si = [ L.vprops[i][:val] for i in n ]
    sk = L.vprops[s][:val]
    return 2J*sk*sum(si)
end

function metropolis!(lattice::IsingData, steps::Integer, T::Float64; showprogress = false)
    p     = Progress(steps)
    vlist = vertices(lattice.final_state)
    β     = 1/T
    L     = lattice
    
    # Store possible exponential values
    z_max = 4^L.generation
    
    # Store a map of neighbors of each vertex
    neighbours_dict = Dict{Int64, Vector{Int64}}()
    for v in vlist
        N = neighbors(L.initial_state, v)
        push!(neighbours_dict, v => N)
    end

    # dE => exp(-β*dE)
    exponential = Dict( [-2z_max:2:2z_max;] .=> exp.(-β .* [-2z_max:2:2z_max;] ) )
    
    for _ in Base.OneTo(steps)
        # Pick a random vertex
        v = rand(vlist)
        # Calculate energy for flipping the spin
        dE = ΔE(L.final_state, v, neighbours_dict[v])
        
        # If new energy is not lower, probablistically flip it
        u = rand()
        if dE < 0
            L.final_state.vprops[v][:val] *= -1
            push!(lattice.spinflip_history, v)
        elseif (u < exponential[Integer(dE)])
            L.final_state.vprops[v][:val] *= -1
            push!(lattice.spinflip_history, v)
        else
            push!(lattice.spinflip_history, -1)
        end
        
        # push!(L.magnetization_history, magnetization(L))
        # push!(L.internalenergy_history, energy(L))
        
        if showprogress
            next!(p)
        end
    end
end
