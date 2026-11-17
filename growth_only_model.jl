using PhysicalConstants.CODATA2018
using Unitful
using DifferentialEquations
using BenchmarkTools
using Revise
using HeterogeneousArrays

include("growth_only_params.jl")

# Initial conditions

function initialize_system(params)
    B_0 = params.B_0
    Ae_0 = params.Ae_0
    Ai_0 = params.Ai_0
    AT_0 = params.AT_0
    AU_0 = params.AU_0
    HeterogeneousVector(Ae=Ae_0, Ai=Ai_0, AT=AT_0, AU=AU_0, B=B_0)
end

u0 = initialize_system(model_params)

abstol_struct = abstol .* oneunit.(u0)

function ode_system!(du, u, params, t)    
    Ai = u.Ai
    Ae = u.Ae
    AT = u.AT
    AU = u.AU
    B = u.B
    
    T = (B * params.T_0 - AT) / (B * params.V_i * params.N_A)

    function r(f_free)
        if f_free > params.f_c
            f_free / (oneunit(params.f_c) - params.f_c)
        else
            zero(params.f_c)
        end
        
    end

    function r(AT, f_c)
        f_free = T / (T + AT)
        r(f_free)
    end

    du.B = r(params.T_0, AT, params.f_c) * B * (1 - B / params.K)
    du.A_e = -params.p * (A_e - Ai)
    du.A_i = params.p * B * (A_e - Ai)- params.k_f * Ai * T + params.k_r * AT - k_uf * Ai * B + k_ur * AU
    du.AT = k_f * Ai * T - params.k_r * AT
    du.AU = params.k_uf * k_f * Ai * T - params.k_r * AT
    du
end

problem = ODEProblem(ode_system!,u0,(zero(model_params.tspan),model_params.tspan),model_params)

# Does not work as autodiff is not yet supported
# sol = solve(problem,AutoTsit5(Rosenbrock23()))
# This does not work either because the solver wants to evaluate real(eltype(prob.b)) and
# eltype(prob.b) is not a concrete type
# sol = solve(problem,AutoTsit5(Rosenbrock23(autodiff=false, concrete_jac=false)))


# Does not work because the solver wants to evaluate similar(u, Complex{eltype(u)})
# sol = solve(problem, RadauIIA3(autodiff=false))

# Do not work for similar reasons
# sol = solve(problem, Trapezoid(autodiff=false))
# sol = solve(problem, ImplicitMidpoint(autodiff=false), dt = 1.0u"s")
# sol = solve(problem, Hairer42(autodiff=false))
# sol = solve(problem, RKC(), saveat=tsave)
# sol = solve(problem, QNDF(autodiff=false), saveat=tsave)
# sol = solve(problem, QBDF(autodiff=false), saveat=tsave)
# sol = solve(problem, EPIRK4s3A(autodiff=false), saveat=tsave, dt = 1.0u"s")
# sol = solve(problem, ESERK5(), saveat=tsave)

# Does work, but is does not provide substantial performance benefits over RK4()
# sol = solve(problem, Tsit5())
# @btime solve(problem, Tsit5())

@time sol = solve(problem,RK4(); abstol=abstol_struct, saveat=tsave)
@btime solve(problem,RK4(); abstol=abstol_struct, saveat=tsave)
@profview solve(problem,RK4();abstol=abstol_struct, saveat=tsave)
