mutable struct Solution
    cost::Float64
    routes::Array{Array{Int}}
end

# build Solution from the variables x
function getsolution(data::DataTOP, x, objval, optimizer)
    A, dim = arcs(data), num_customers(data) + 2
    adj_list = [[] for i in 1:dim]
    for a in A
        val = get_value(optimizer, x[a])
        if val > 0.5
            push!(adj_list[a[1] + 1], a[2])
        end
    end
    visited, routes = [false for i in 2:dim], []
    for i in adj_list[1]
        if !visited[i]
            r, prev = [], 0
            next, prev = i, 0
            maxit, it = dim - 2, 0
            while next != dim - 1 && it < maxit
                length(adj_list[next + 1]) != 1 && error(
                    "Problem trying to recover the route from the x values. " *
                    "Customer $next has $(length(adj_list[next+1])) outcoming arcs.",
                )
                visited[next] && error("Customer $next is visited more than once")
                visited[next] = true
                push!(r, next)
                aux = next
                next, prev = adj_list[next + 1][1], aux
                it += 1
            end
            (next != dim - 1) && error(
                "Problem trying to recover the route from the x values. " *
                "Route cannot be recovered because the return to depot is never reached",
            )
            push!(routes, r)
        elseif i != dim - 1
            error("Customer $i is visited more than once")
        end
    end

    return Solution(objval, routes)
end

function print_routes(solution)
    for (i, r) in enumerate(solution.routes)
        print("Route #$i: ")
        for j in r
            print("$j ")
        end
        println()
    end
end
