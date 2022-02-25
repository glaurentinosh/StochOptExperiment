#include("bellman-example-max-julia.jl")
include("eval-bellmans.jl")
using Profile
using Distributed

#gr()
#if (nprocs() == 1) addprocs(4) end
#printstyled("NUM_PROCS : $(nprocs())\n"; color = :blue)
#println("NUM_THREADS : $(Threads.nthreads())")

#= params = [
[100,40,40,5,10,4],
[100,40,50,5,10,4],
[100,40,60,5,10,4],
[100,40,70,5,10,4],
[100,40,80,5,10,4],
[100,40,90,5,10,4],
[100,50,40,5,10,4],
[100,50,50,5,10,4],
[100,50,60,5,10,4],
[100,50,70,5,10,4],
[100,50,80,5,10,4],
[100,50,90,5,10,4],
]

for param in params
    call(param...)
end =#

#= for s in 100:50:400
    for u in 20:10:0.9*s
        for w in 20:10:0.9*s
            println("s, u, w : $(s), $(u), $(w)")
            call(s, u, w, 5, 10, 8)
        end
    end
end
=#

# 100,80,80
function calling(args...)
    diff = 0
    smax, umax, wmax, pmax, Pmax = 0,0,0,0,0
    x = 0
    ustep = 20
    sstart = length(args) == 0 ? 50 : parse(Int64, args[1])

    bellman_function = Array{Real}[]
    hd_bellman_function = Array{Real}[]
    
    for s in sstart:50:300
        for u in 0.2*s:ustep:0.9*s
            for w in 0.2*s:ustep:3*s
                for p in 2:2:10
                    for P in 5:5:25
                        println("s, u, w, p, P : $(s), $(u), $(w), $(p), $(P)")
                        x = @time call(s, u, w, p, P, 8, bellman_function, hd_bellman_function)
                        if (diff < x)
                            diff = x
                            smax, umax, wmax, pmax, Pmax = s, u, w, p, P
                            printstyled("Found diff $(x) for s, u, w, p, P = $(s), $(u), $(w), $p, $P\n"; color = :green )
                        else
                            printstyled("Found diff $(x) for s, u, w, p, P = $(s), $(u), $(w), $p, $P "; color = :red )
                            printstyled("still $(diff) for s, u, w, p, P = $(smax), $(umax), $(wmax), $pmax, $Pmax\n"; color = :green )
                        end
                    end
                end
            end
        end
    end
end

@time calling(ARGS...)