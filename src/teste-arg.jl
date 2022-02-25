print("Oi"*ARGS[1]*"\n")
res = 0
for arg in ARGS
    global res += parse(Int64, arg)
end
@time print(res)