mutable struct Vertex
    id_vertex::Int
    pos_x::Float64
    pos_y::Float64
    profit::Int
end

# Directed graph
mutable struct InputGraph
    V::Array{Vertex} # set of vertices
    A::Array{Tuple{Int64,Int64}} # set of arcs
    time::Dict{Tuple{Int64,Int64},Float64} # time for each arc
end

mutable struct DataTOP
    n::Int
    G::InputGraph
    T::Float64
    F::Int
end

# Euclidian distance
function distance(v1::Vertex, v2::Vertex)
    x_sq = (v1.pos_x - v2.pos_x)^2
    y_sq = (v1.pos_y - v2.pos_y)^2
    return sqrt(x_sq + y_sq)
end

function readTOPData(path_file::String)
    # STEP 1 : pushing data in a vector.
    data = Array{Any,1}()
    open(path_file) do file
        for line in eachline(file)
            for peaceofdata in split(line)
                push!(data, String(peaceofdata))
            end
        end
    end

    n = parse(Int, data[1])
    F = parse(Int, data[2])
    T = parse(Float64, data[3])

    vertices = Vertex[]
    for i in 0:(n + 1)
        offset = 3 + i * 3
        x = parse(Float64, data[offset + 1])
        y = parse(Float64, data[offset + 2])
        p = parse(Int, data[offset + 3])
        push!(vertices, Vertex(i, x, y, p))
    end

    A = Tuple{Int64,Int64}[]
    time = Dict{Tuple{Int64,Int64},Float64}()

    function add_arc!(i, j)
        push!(A, (i, j))
        return time[(i, j)] = distance(vertices[i + 1], vertices[j + 1])
    end
    for i in 1:n
        #arc from depot
        add_arc!(0, i)
        #arc to depot
        add_arc!(i, n + 1)
        for j in 1:n
            if (i != j)
                add_arc!(i, j)
            end
        end
    end

    return DataTOP(n, InputGraph(vertices, A, time), T, F)
end

arcs(data::DataTOP) = data.G.A # return set of arcs
function t(data, a)
    if !(haskey(data.G.time, a))
        return Inf
    end
    return data.G.time[a]
end
num_customers(data::DataTOP) = data.n # return number of requests
p(data::DataTOP, i) = data.G.V[i + 1].profit # return profit of i
fleet_size(data::DataTOP) = data.F
max_duration(data::DataTOP) = data.T
