echo `date`
julia bellman-example-max-julia.jl 200 80 50 5 10 8 &
julia bellman-example-max-julia.jl 100 40 50 5 10 8 &
wait # do not return before background tasks are complete
echo `date`