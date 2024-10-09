include("array_interface.jl")

function create_simulation_variables()
    A = 1.0u"m"   # meters
    T = 2.0u"s"   # seconds
    AT = 3.0u"m/s"  # meters per second
    B = @SVector [1.0u"m", 2.0u"m", 3.0u"m"]  # vector of meters

    return DiscreteSimulationVariables(A, T, AT, B)
end

val = create_simulation_variables()

val .+ val

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

N = 1e10

f(val, N)

g(val,N)

@time f(val, N)

@time g(val, N)

Profile.init()

@profview_allocs f(val,N)

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
