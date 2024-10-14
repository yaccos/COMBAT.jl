include("array_interface.jl")
using BenchmarkTools

function create_simulation_variables()
    A = 1.0u"m"   # meters
    T = 2.0u"s"   # seconds
    AT = 3.0u"m/s"  # meters per second
    B = [1.0, 2.0, 3.0]u"m"  # vector of meters

    return DiscreteSimulationVariables(A, T, AT, B)
end

val = create_simulation_variables()

val .+= val


function f(val,N)
    for i in 1:N
        val .+ val
    end
end

function g(val,N)
    for i in 1:N
        val .+= val
    end
end


N = 1e6

g(val,N)

f(val,N)

@btime g(val, N)

@btime f(val, N)

@btime g(val.B, N)

@btime g(val.B, N)

@profview g(val,N) C=true

using Cthulhu
using Revise

@descend broadcast!(+,val,val)

using JET

@report_opt g(val, N)



@btime f(val,N)

@profview f(val,N)

Profile.init()

@profview_allocs g(val,N)

using ProfileView

Profile.init()

Profile.init(n=1000,delay=1e-6)

ProfileView.@profview g(val,N)

@profview_allocs g(val,N)

@profview_allocs f(val,N)

@time f(val,N)

@time g(val,N)

@profview_allocs g(val,N)

using Profile

profile = Profile.Allocs.@profile sample_rate=1 g(val,N)

using PProf

PProf.Allocs.pprof(from_c = true)
