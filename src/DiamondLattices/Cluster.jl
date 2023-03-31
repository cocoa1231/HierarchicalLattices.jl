using Dictionaries
import Base:@kwdef

@kwdef mutable struct WolffData{T}
    lattice::T
    magnetization_history::Vector{Float64}
    internalenergy_history::Vector{Float64}
    saveinterval::Int64
    thermalization_steps::Int64
end

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

function wolff!(WD::WolffData{L}, nsteps, T; showprogress = false, verbose = false) where L
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

    Pthermalize = Progress(WD.thermalization_steps, enabled = showprogress, showspeed = true)
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

    Ptrial = Progress(nsteps, enabled = showprogress, showspeed = true)
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