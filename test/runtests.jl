using HierarchicalLattices
using Test
using StatsBase

@testset "Critical point test with plottable lattices" begin
    g = DiamondLattice(diamond_ising_lattice(5, :zero), 5)
    ID = IsingData(g, Float64[], Float64[], Int64[])
    nsteps = 10000 * 687
    N = 687

    T = [1., 1.6, 3.]
    U = zeros(length(T))
    U2 = zeros(length(T))

    for (idx, t) in enumerate(T)
        metropolis!(ID, nsteps, t)
        fill_data!(ID, :U)

        U[idx] = mean(ID.internalenergy_history[end-1000N:2N:end])
        U2[idx] = mean(ID.internalenergy_history[end-1000N:2N:end] .^2)

        ID.internalenergy_history = Float64[]
        ID.spinflip_history = Int64[]
    end

    specheat = U2 .- (U .^2)
    maxidx = findfirst(==(maximum(specheat)), specheat)

    @test maxidx == 2
end

@testset "Critical point test with non-plottable b = 2 lattices" begin
    g = DiamondLattice(diamond_ising_lattice(5, 2, :zero), 5)
    ID = IsingData(g, Float64[], Float64[], Int64[])
    nsteps = 10000 * 687
    N = 687

    T = [1., 1.6, 3.]
    U = zeros(length(T))
    U2 = zeros(length(T))

    for (idx, t) in enumerate(T)
        metropolis!(ID, nsteps, t)
        fill_data!(ID, :U)

        U[idx] = mean(ID.internalenergy_history[end-1000N:2N:end])
        U2[idx] = mean(ID.internalenergy_history[end-1000N:2N:end] .^2)

        ID.internalenergy_history = Float64[]
        ID.spinflip_history = Int64[]
    end

    specheat = U2 .- (U .^2)
    maxidx = findfirst(==(maximum(specheat)), specheat)

    @test maxidx == 2
end

