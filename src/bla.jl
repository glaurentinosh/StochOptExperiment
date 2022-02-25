function laplacian_bad(lap_x::Array{Float64,2}, x::Array{Float64,2})
    nr, nc = size(x)
    for ir = 2:nr-1, ic = 2:nc-1 # bad loop nesting order
        lap_x[ir, ic] =
            (x[ir+1, ic] + x[ir-1, ic] +
            x[ir, ic+1] + x[ir, ic-1]) - 4*x[ir, ic]
    end
end

function laplacian_good(lap_x::Array{Float64,2}, x::Array{Float64,2})
    nr,nc = size(x)
    for ic = 2:nc-1, ir = 2:nr-1 # good loop nesting order
        lap_x[ir,ic] =
            (x[ir+1,ic] + x[ir-1,ic] +
            x[ir,ic+1] + x[ir,ic-1]) - 4*x[ir,ic]
    end
end

function laplacian_good_nocheck(lap_x::Array{Float64,2}, x::Array{Float64,2})
    nr,nc = size(x)
    for ic = 2:nc-1, ir = 2:nr-1 # good loop nesting order
        @inbounds begin lap_x[ir,ic] = # no array bounds checking
            (x[ir+1,ic] +  x[ir-1,ic] +
            x[ir,ic+1] + x[ir,ic-1]) - 4*x[ir,ic]
        end
    end
end

function main_test(nr, nc)
    field = zeros(nr, nc)
    for ic = 1:nc, ir = 1:nr
        if ir == 1 || ic == 1 || ir == nr || ic == nc
            field[ir,ic] = 1.0
        end
    end
    lap_field = zeros(size(field))

    t = @elapsed laplacian_bad(lap_field, field)
    println(rpad("laplacian_bad", 30), t)
    
    t = @elapsed laplacian_good(lap_field, field)
    println(rpad("laplacian_good", 30), t)
    
    t = @elapsed laplacian_good_nocheck(lap_field, field)
    println(rpad("laplacian_good no check", 30), t)
end