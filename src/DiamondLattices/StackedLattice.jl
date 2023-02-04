mutable struct StackedDiamondLattice
    coupling
    initial_state
    final_state
    adjacentspins
end
function StackedDiamondLattice(coupling::Number, depth::Number, order::Integer, initstate::Symbol)
    lattices = [
        DiamondLattice(diamond_ising_lattice(order, initstate), order) for _ in 1:depth]
    adjacentspins = Dict(vertices(lattices[1].final_state) .=> (x -> neighbors(lattices[1].final_state, x)).(vertices(lattices[1].final_state)))
    D = StackedDiamondLattice(coupling, deepcopy(lattices), deepcopy(lattices), adjacentspins)
end

function energy(S::StackedDiamondLattice; state = :final)
    if state == :final
        D = S.final_state
    elseif state == :inital
        D = S.initial_state
    end

    lattice = (x -> x.final_state).(D)
    total_K = 0.
    total_J = 0.
    depth = length(lattice)

    for level in 1:depth
        L = lattice[level]
        spins = vertices(L)
        bonds = edges(L)
        
        for spin in spins
            total_K += L.vprops[spin][:val]*lattice[mod1(level-1, depth)].vprops[spin][:val]
            total_K += L.vprops[spin][:val]*lattice[mod1(level+1, depth)].vprops[spin][:val]
        end
        
        for bond in bonds
            total_J += L.vprops[bond.src][:val]*L.vprops[bond.dst][:val]
        end
    end
    
    return -(S.coupling*total_K/2) - total_J
end

function ΔE(S::StackedDiamondLattice, site, n; state = :final)
    level, spin = site
    if state == :final
        lattice = S.final_state
    elseif state == :initial
        lattice = S.initial_state
    end
    depth = length(lattice)
    si  = [ lattice[level].vprops[i][:val] for i in n ]
    sKi = [ lattice[mod1(level+i, depth)].vprops[spin][:val] for i in [+1, -1] ]
    sk = lattice[level].vprops[spin][:val]
    ΔJ = 2J*sk*sum(si)
    ΔK = 2K*sk*sum(sKi)
    
    return ΔJ + ΔK * S.coupling
end
function magnetization(S::StackedDiamondLattice; state = :final)
    if state == :final
        lattice = S.final_state
    elseif state == :UndefInitializer
        lattice = S.initial_state
    end

    M = 0.
    for L in lattice
        M += sum([L.final_state.vprops[v][:val] for v in vertices(L.final_state)])
    end
    
    return M
end
numberofspins(S::StackedDiamondLattice) = sum([length(D.final_state.vprops) for D in S.final_state])