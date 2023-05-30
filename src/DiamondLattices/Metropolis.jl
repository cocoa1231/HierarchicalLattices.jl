"""
Main type used by the `metropolis!` function. The fields are:

- `lattice`: Lattice to evolve.
- `magnetization_history::Vector{Float64}`: Every index `i` corresponds to the
  magnetization at the state `i`.
- `internalenergy_history::Vector{Float64}`: Every index `i` corresponds to the
  internal energy at the state `i`.
- `spinflip_history::Vector`:  Every index `i` corresponds to the spin flipped
  in the previous Metropolis step `i-1`.
"""
mutable struct IsingData
    lattice
    magnetization_history::Vector{Float64}
    internalenergy_history::Vector{Float64}
    spinflip_history::Vector
end

IsingData(D::DiamondLattice) = IsingData(D, Float64[], Float64[], Int64[])

"""
Evolves the system the given number of steps.

- `data::IsingData`: Lattice to evolve.
- `steps::Integer`: Number of steps to evolve the system for.
- `T::Float64`: Temperature at which to evolve the system.
- `showprogress` (default `false`): Add a progress bar.
- `progressoutput` (default `stdout`): Write the progress bar output to this IO.
"""
function metropolis!(data::IsingData, steps::Integer, T::Float64; showprogress = false, progressoutput = stdout)
    metropolis!(data.lattice, data, steps, T; showprogress = showprogress, progressoutput = progressoutput)
end

function metropolis!(::DiamondLattice, data::IsingData, steps::Integer, T::Float64; showprogress = false, progressoutput = stdout)
    p     = Progress(steps, enabled = showprogress, showspeed = true, output = progressoutput)
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

        
        next!(p)
    end
end

function metropolis!(::StackedDiamondLattice, data::IsingData, steps::Integer, T::Float64; showprogress = false, progressoutput = stdout)
    lattice = data.lattice
    P = Progress(steps, enabled = showprogress, showspeed = true, output = progressoutput)
    β = 1/T
    K = 2
    z_max = 4^lattice.generation
    viter = vertices(lattice.final_state)
    exponential = Dict{Any, Any}( [-2*(z_max-K):2:2*(z_max+K);] .=> exp.(-β .* [-2*(z_max-K):2:2*(z_max+K);]) )
    
    spinsflipped = zeros(Integer, steps)

    for idx in Base.OneTo(steps)
        # Pick a random vertex
        v = rand(viter)
        
        # Calculate dE
        dE = ΔE(lattice, v)
        if !(dE in keys(exponential))
            push!(exponential, dE => exp(-β*dE))
        end

        # If new energy is not lower, flip it with probability exp(-βΔE)
        if dE < 0
            lattice.final_state.vprops[v][:val] *= -1
            spinsflipped[idx] = v
        elseif rand() < exponential[dE]
            lattice.final_state.vprops[v][:val] *= -1
            spinsflipped[idx] = v
        else
            spinsflipped[idx] = -1
        end
        
        next!(P)
    end

    append!(data.spinflip_history, spinsflipped)
end

"""
Fill the internal energy or magnetization history of a lattice evolved using the
`metropolis!` function. Inputs:

- `data::IsingData`: Lattice to fill data for.
- `variable::Symbol`: Which variable to fill. Possible values are `:U` for
  internal energy and `:M` for magnetization.
- `showprogress` (default `false`): Add a progress bar.
- `progressoutput` (default `stdout`): Write the progress bar output to this IO.
"""
function fill_data!(data::IsingData, variable::Symbol; showprogress = false, progressoutput = stdout)
    if variable == :M
        return _fill_M_history!(data.lattice, data, showprogress = showprogress, progressoutput = progressoutput)
    elseif variable == :U
        return _fill_U_history!(data.lattice, data, showprogress = showprogress, progressoutput = progressoutput)
    else
        throw(ArgumentError("Data parameter $(string(data)) not implimented!"))
    end
end
