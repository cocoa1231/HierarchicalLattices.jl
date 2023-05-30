using Dictionaries
import Base:@kwdef

"""
Main type used by the `wolff!` function. Fields:

- `lattice::T`: Lattice to evolve using the Wolff algorithm.
- `magnetization_history::Vector{Float64}`: Every index `i` corresponds to the
  magnetization of the system at the wolff step `i*saveinterval`.
- `internalenergy_history::Vector{Float64}`: Every index `i` corresponds to the
  internal energy of the system at the wolff step `i*saveinterval`.
- `saveinterval::Int64`: Save the state of the system after this many steps.
- `thermalization_steps::Int64`: Steps the system for this many steps before 
  starting data collection.
"""
@kwdef mutable struct WolffData{T}
    lattice::T
    magnetization_history::Vector{Float64}
    internalenergy_history::Vector{Float64}
    saveinterval::Int64
    thermalization_steps::Int64
end

"""
Step the system for one step of the Wolff algorithm. Inputs:

- `lattice::DiamondLattice`: Lattice to evolve. This function is specialized for
  the diamond lattice.
- `P_add::Float64`: Probability of adding the a spin to the cluster.
- `adjacentmap::Dictionary{Int64, Vector{Int64}}`: Adjacency map of the system.
  Will be changed to use the adjacency matrix.
"""
function wolffstep!(lattice::DiamondLattice, P_add::Float64, adjacentmap::Dictionary{Int64, Vector{Int64}})
    state = lattice.final_state
    viter = vertices(state)
    buff = zeros(Int64, 0)
    seed = rand(viter)
    cluster = Int64[seed]
    clustersign = get_prop(state, seed, :val)


    set_prop!(state, seed, :val, -clustersign)

    for n in adjacentmap[seed]
        if get_prop(state, n, :val) == clustersign &&
            rand() < P_add
            
            push!(buff, n)
            push!(cluster, n)
        end
    end

    while !isempty(buff)
        target = pop!(buff)
        set_prop!(state, target, :val, -clustersign)
        for n in filter(!in(cluster), adjacentmap[target])
            if get_prop(state, n, :val) == clustersign &&
                rand() < P_add
                push!(buff, n)
                push!(cluster, n)
            end
        end
    end

    return cluster
end

"""
Step the system for one step of the Wolff algorithm. Inputs:

- `lattice::StackedDiamondLattice`: Lattice to evolve. This function is specialized for
  the stacked diamond lattice.
- `probdict::Dictionary{Float64, Float64}`: Probability of adding the a spin to
  the cluster. Maps the bond (or edge) strength to the probability.
- `adjacentmap::Dictionary{Int64, Vector{Int64}}`: Adjacency map of the system.
  Will be changed to use the adjacency matrix.
"""
function wolffstep!(lattice::StackedDiamondLattice, probdict::Dictionary{Float64, Float64}, 
        adjacentmap::Dictionary{Int64, Vector{Int64}})

    state = lattice.final_state
    viter = vertices(state)
    buff = zeros(Int64, 0)
    seed = rand(viter)
    cluster = Int64[seed]
    clustersign = get_prop(state, seed, :val)
    
    set_prop!(state, seed, :val, -clustersign)
    
    for n in adjacentmap[seed]
        if get_prop(state, n, :val) == clustersign &&
            rand() < probdict[get_prop(state, seed, n, :weight)]
            
            push!(buff, n)
            push!(cluster, n)
        end
    end
    
    spinsflipped = 1
    while !isempty(buff)
        target = pop!(buff)
        set_prop!(state, target, :val, -clustersign)
        spinsflipped += 1
        for n in filter(!in(cluster), adjacentmap[target])
            if get_prop(state, n, :val) == clustersign &&
                rand() < probdict[get_prop(state, target, n, :weight)]
                push!(buff, n)
                push!(cluster, n)
            end
        end
    end
    
    return cluster, spinsflipped
end

"""
Run the Wolff algorithm for a given number of steps at a temperature. Inputs:

- `WD::WolffData{L}` `WolffData` with lattice type `L`.
- `nsteps`: Number of steps to evolve the system for.
- `T`: Temperature at which to evolve the system.
- `showprogress` (default `false`): Add a progress bar.
- `verbose` (default `false`): Look at the argument name.
- `progressoutput` (default `stdout`): Write the progress output to the given
  IO.
"""
function wolff!(WD::WolffData{L}, nsteps, T; showprogress = false, verbose = false, progressoutput = stdout) where L
    β = 1/T

    if L == DiamondLattice
        P_add = 1 - exp(-2β)
    elseif L == StackedDiamondLattice
        P_add = Dict(1.0 => 1 - exp(-2β), WD.lattice.stackingweight => 1 - exp(-2β*WD.lattice.stackingweight)) |> Dictionary
    end

    if verbose
        @info "Creating adjacent spin map"
    end

    adjacentmap = Dictionary(Dict(
        vertices(WD.lattice.final_state) .=>  (x -> neighbors(WD.lattice.final_state, x)).(vertices(WD.lattice.final_state))
    ))

    Pthermalize = Progress(WD.thermalization_steps, enabled = showprogress, showspeed = true, output = progressoutput)
    if verbose
        @info "Thermalizing system for $(WD.thermalization_steps) steps."
    end
    for _ in 1:WD.thermalization_steps
        wolffstep!(WD.lattice, P_add, adjacentmap)
        next!(Pthermalize)
    end

    if verbose
        @info "Starting simulation for $nsteps steps."
        @info "T = $T."
        @info "Save interval = $(WD.saveinterval) steps."
    end

    Ptrial = Progress(nsteps, enabled = showprogress, showspeed = true, output = progressoutput)
    energies = zeros(floor(Int64, nsteps / WD.saveinterval))
    magnetizations = zeros(floor(Int64, nsteps / WD.saveinterval))
    for i in 1:nsteps
        wolffstep!(WD.lattice, P_add, adjacentmap)
        
        if i % WD.saveinterval == 0
            idx = Int(i / WD.saveinterval)
            energies[idx] = energy(WD.lattice)
            magnetizations[idx] = magnetization(WD.lattice)
        end
        next!(Ptrial)
    end

    append!(WD.internalenergy_history, energies)
    append!(WD.magnetization_history, magnetizations)
    
    return energies, magnetizations
end
