# HierarchicalLattices.jl

This library provides several functions to generate hierarchical lattices
[(Griffiths,
Kaufman)])(<https://journals.aps.org/prb/abstract/10.1103/PhysRevB.26.5022>).
The primary idea to do this is to define a transformation on a part of a graph
and then iterate said transformation across all the subgraphs. This library
utilizes Graphs.jl to construct graphs and MetaGraphs.jl to store spin metadata
on the vertices and bond strength on the edges. At the moment, this library also
includes implementations of the Metropolis and Wolff algorithms for simulating
spin-1/2 Ising systems on these lattices, however this will be moved to a
separate library in the future. Implementations and interfaces are in a nascent
stage, thus expect breaking changes in both. At the moment, only the diamond
hierarchical lattice has been implemented.

## Code structure

For every lattice type, it's implementation will be stored in `src/NameLattice`.
Any MCMC algorithm can also be stored in this directory. At the moment the
`DiamondLattice` implements lattice generation methods in `Lattice.jl`

## Usage

> **Warning**<br>
This documentaiton may be out of date if breaking commits are made after
29 May, 2023. Raise an issue in such a case.

> **Warning**<br>
This library currently has a plotting recipe for a `MetaGraph`, which means it
can clash with another graph plotting library implementing the same. Use with
caution.

There are two main data structures exported for the diamond lattice, the
`DiamondLattice` and `StackedDiamondLattice` types. Refer to the doc strings for
how to use each function. A typical use case of simulating a lattice can be as
follows.

```julia
julia> using HierarchicalLattices

julia> order = 5; b = 2; initstate = :zero
:zero

julia> lattice = DiamondLattice(diamond_ising_lattice(order, b, initstate), order)
DiamondLattice(5, {684, 1024} undirected Int64 metagraph with Float64 weights defined by :weight (default weight 1.0), {684, 1024} undirected Int64 metagraph with Float64 weights defined by :weight (default weight 1.0), 1.0)

julia> ID = IsingData(lattice, Float64[], Float64[], Int64[])
IsingData(DiamondLattice(5, {684, 1024} undirected Int64 metagraph with Float64 weights defined by :weight (default weight 1.0), {684, 1024} undirected Int64 metagraph with Float64 weights defined by :weight (default weight 1.0), 1.0), Float64[], Float64[], Int64[])

julia> metropolis!(ID, 1000, 1.64)

julia> ID.lattice.initial_state.vprops
Dict{Int64, Dict{Symbol, Any}} with 684 entries:
  319 => Dict(:val=>1)
  185 => Dict(:val=>1)
  420 => Dict(:val=>1)
  525 => Dict(:val=>1)
  365 => Dict(:val=>1)
  638 => Dict(:val=>1)
  263 => Dict(:val=>1)
  422 => Dict(:val=>1)
  242 => Dict(:val=>1)

julia> ID.lattice.final_state.vprops
Dict{Int64, Dict{Symbol, Any}} with 684 entries:
  319 => Dict(:val=>1)
  185 => Dict(:val=>1)
  420 => Dict(:val=>-1)
  525 => Dict(:val=>1)
  365 => Dict(:val=>1)
  638 => Dict(:val=>1)
  263 => Dict(:val=>1)
  422 => Dict(:val=>1)
  242 => Dict(:val=>-1)
```

Similar usage for a `StackedDiamondLattice` can be seen (this time with the
Wolff algorithm).

```julia
julia> using HierarchicalLattices

julia> order = 5; depth = 15; stackingweight = 1.5; initstate = :zero
:zero

julia> stackedlattice = StackedDiamondLattice(order, depth, stackingweight, initstate)
StackedDiamondLattice(5, 1.5, {10260, 24936} undirected Int64 metagraph with Float64 weights defined by :weight (default weight 1.0), {10260, 24936} undirected Int64 metagraph with Float64 weights defined by :weight (default weight 1.0), Dict(4986 => Dict(1.5 => [4302, 5670], 1.0 => [4789, 4845]), 7329 => Dict(1.5 => [6645, 8013], 1.0 => [6860, 6975]), 4700 => Dict(1.5 => [4016, 5384], 1.0 => [4137, 4250]), 4576 => Dict(1.5 => [3892, 5260], 1.0 => [4122, 4160]), 7144 => Dict(1.5 => [6460, 7828], 1.0 => [6845, 6950]), 6073 => Dict(1.5 => [5389, 6757], 1.0 => [5506, 5619]), 2288 => Dict(1.5 => [1604, 2972], 1.0 => [2054, 2128]), 1703 => Dict(1.5 => [1019, 2387], 1.0 => [1375, 1494]), 1956 => Dict(1.5 => [1272, 2640], 1.0 => [1400, 1492]), 8437 => Dict(1.5 => [7753, 9121], 1.0 => [8210, 8281])…))

julia> WD = WolffData(lattice = stackedlattice, magnetization_history=Float64[], internalenergy_history=Float64[], thermalization_steps=1000, saveinterval=20)
WolffData{StackedDiamondLattice}(StackedDiamondLattice(5, 1.5, {10260, 24936} undirected Int64 metagraph with Float64 weights defined by :weight (default weight 1.0), {10260, 24936} undirected Int64 metagraph with Float64 weights defined by :weight (default weight 1.0), Dict(4986 => Dict(1.5 => [4302, 5670], 1.0 => [4789, 4845]), 7329 => Dict(1.5 => [6645, 8013], 1.0 => [6860, 6975]), 4700 => Dict(1.5 => [4016, 5384], 1.0 => [4137, 4250]), 4576 => Dict(1.5 => [3892, 5260], 1.0 => [4122, 4160]), 7144 => Dict(1.5 => [6460, 7828], 1.0 => [6845, 6950]), 6073 => Dict(1.5 => [5389, 6757], 1.0 => [5506, 5619]), 2288 => Dict(1.5 => [1604, 2972], 1.0 => [2054, 2128]), 1703 => Dict(1.5 => [1019, 2387], 1.0 => [1375, 1494]), 1956 => Dict(1.5 => [1272, 2640], 1.0 => [1400, 1492]), 8437 => Dict(1.5 => [7753, 9121], 1.0 => [8210, 8281])…)), Float64[], Float64[], 20, 1000)

julia> wolff!(WD, 100, 1.64; showprogress = true)
Progress: 100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| Time: 0:02:31 ( 0.15  s/it)
Progress: 100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| Time: 0:00:14 ( 0.15  s/it)
([-29469.0, -29387.0, -29486.0, -29497.0, -29386.0], [10194.0, 10176.0, 10204.0, 10202.0, 10174.0])
```

The `StackedDiamondLattice` constructor can construct a `MetaGraph` that is
stacking diamond lattices and provides you with a single resulting `MetaGraph`.
This has the advantage that, for the most part, the algorithms do not change,
and one can just implement the algorithms for a general `MetaGraph` with spin
stored in the `:val` property and interaction strength stored in the edge
weight (if nonexistent, then this is the default weight of the graph edges).
