mutable struct IsingData
    lattice
    magnetization_history::Vector{Float64}
    internalenergy_history::Vector{Float64}
    spinflip_history::Vector
end

IsingData(D::DiamondLattice) = IsingData(D, Float64[], Float64[], Int64[])
function IsingData(D::DiamondLattice)
    nbonds = edges(D.final_state) |> collect |> length
    g = Int(log(4, nbonds))
    return IsingData(D)
end


function metropolis!(data::IsingData, steps::Integer, T::Float64; showprogress = false)
    metropolis!(data.lattice, data, steps, T; showprogress = showprogress)
end

function metropolis!(::DiamondLattice, data::IsingData, steps::Integer, T::Float64; showprogress = false)
    p     = Progress(steps)
    lattice = data.lattice
    vlist = vertices(lattice.final_state)
    β     = 1/T
    D     = data
    
    # Store possible exponential values
    z_max = 4^lattice.generation
    
    # Store a map of neighbors of each vertex
    neighbours_dict = Dict{Int64, Vector{Int64}}()
    for v in vlist
        N = neighbors(lattice.final_state, v)
        push!(neighbours_dict, v => N)
    end

    # dE => exp(-β*dE)
    exponential = Dict( [-2z_max:2:2z_max;] .=> exp.(-β .* [-2z_max:2:2z_max;] ) )
    
    for _ in Base.OneTo(steps)
        # Pick a random vertex
        v = rand(vlist)
        # Calculate energy for flipping the spin
        dE = ΔE(lattice, v, neighbours_dict[v])
        
        # If new energy is not lower, probablistically flip it
        u = rand()
        if dE < 0
            lattice.final_state.vprops[v][:val] *= -1
            push!(data.spinflip_history, v)
        elseif (u < exponential[Integer(dE)])
            lattice.final_state.vprops[v][:val] *= -1
            push!(data.spinflip_history, v)
        else
            push!(data.spinflip_history, -1)
        end

        
        if showprogress
            next!(p)
        end
    end
end

function metropolis!(::StackedDiamondLattice, data::IsingData, steps::Integer, T::Float64; showprogress = false)
    lattice = data.lattice
    P = Progress(steps)
    β = 1/T
    K = 2
    z_max = 4^lattice.generation
    viter = vertices(lattice.final_state)
    exponential = Dict( [-2*(z_max-K):2:2*(z_max+K);] .=> exp.(-β .* [-2*(z_max-K):2:2*(z_max+K);]) )
    
    for _ in Base.OneTo(steps)
        # Pick a random vertex
        v = rand(viter)
        
        # Calculate dE
        dE = ΔE(lattice, v)
        
        # If new energy is not lower, flip it with probability exp(-βΔE)
        if dE < 0
            lattice.final_state.vprops[v][:val] *= -1
            push!(data.spinflip_history, v)
        elseif rand() < exponential[Integer(dE)]
            lattice.final_state.vprops[v][:val] *= -1
            push!(data.spinflip_history, v)
        else
            push!(data.spinflip_history, -1)
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
function fill_data!(data::IsingData, variable::Symbol; showprogress = false)
    if variable == :M
        return _fill_M_history!(data.lattice, data, showprogress = showprogress)
    elseif variable == :U
        return _fill_U_history!(data.lattice, data, showprogress = showprogress)
    else
        throw(ArgumentError("Data parameter $(string(data)) not implimented!"))
    end
end