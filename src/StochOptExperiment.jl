module StochOptExperiment

Smax = 200
#Smin = 0

# Possible decisions
Umax = 80
Umin = 0

# Possible scenarios
Wmax = 50
Wmin = 0

# Instantaneous cost function
p = 5
Pf = 10

function L(u, s)
    return p*u*u
end



end # module
