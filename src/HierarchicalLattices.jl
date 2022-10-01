module HierarchicalLattices

using Graphs
using MetaGraphs
using RecipesBase
include("Metropolis.jl")


# Diamond Lattice data
export  diamond_order_zero_transform!,
        make_diamond_lattice0,
        diamond_lattice

# Metropolis Stuff
export  metropolis!,
        energy,
        ΔE,
        magnetization

# Main data structure
export IsingData

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

function diamond_order_zero_transform!(lattice, edge; scale = √3/2)
    
    e = edge
    src = lattice.vprops[e.src][:loc]
    dst = lattice.vprops[e.dst][:loc]

    Δz = dst - src
    new_left  = (Δz * 1im) *  scale
    new_right = (Δz * -1im) * scale
    disp = src + Δz/2
    new_pts = [new_left, new_right]
    add_vertices!(lattice, 2)
    for (idx, new_addition) in enumerate(vertices(lattice)[end-1:end])
        set_props!(lattice, new_addition, Dict(
                :loc => new_pts[idx] + disp,
                :val => false
            ))
        add_edge!(lattice, e.src, new_addition)
        add_edge!(lattice, new_addition, e.dst)
    end

    rem_edge!(lattice, e) 

end

function make_diamond_lattice0()
    order_zero = SimpleDiGraph(2, 1) |> MetaGraph
    locations = [ 0. + 1im, 0. - 1im ]
    values = [false, false]
    for (idx, vertex) in enumerate(vertices(order_zero))

        set_props!(order_zero, vertex, Dict(
                :loc => locations[idx],
                :val => values[idx]
            ))
    end
    return order_zero
end

function diamond_lattice(order::Integer)
    oz = make_diamond_lattice0()
    
    for i in 1:order
        for e in collect(edges(oz))
            diamond_order_zero_transform!(oz, e; scale = 0.5^i)
        end
    end
    
    return oz
end

end # module HierarchicalLattices