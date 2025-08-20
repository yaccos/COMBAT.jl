using PhysicalConstants.CODATA2018
using Unitful
using DifferentialEquations
using BenchmarkTools
using Revise
include("heterogeneous_vector.jl")
include("params.jl")

# Initial conditions

function initialize_system(params)
    A_0 = params.A
    T_0 = zero(A_0)
    AT_0 = zero(A_0)
    B_0 = params.B_0
    B_start = zeros(typeof(B_0),params.n+1)
    B_start[1] = B_0
    HeterogeneousVector(A=A_0,T=T_0,AT=AT_0,B=B_start)
end

u0 = initialize_system(model_params)

abstol_struct = abstol .* oneunit.(u0)

function ode_system!(du, u, p, t)    
    B = u.B
    # We can keep dB as a separate variable since it points to an array,
    # but we must explitly reference dA, du.T, du.AT since they would be
    # copied instead of being pointed to
    # For the variables in u, it does not hurt to make a copy since we
    # are not going to assign back into this structure anyway
    dB = du.B
    A = u.A
    T = u.T
    AT = u.AT

    index_start = firstindex(B)
    index_end = lastindex(B)
    carrying_coefficient = (p.C - sum(B)) / p.C
    binding_coefficient = p.k_f/(p.V*p.N_A)
    # Translates between the B array indicies and the number of antibiotic molecules bound in each compartment
    # This is an abstraction intended to avoid off-by-one index errors and support both one and zero indexed versions
    # of arrays
    bound_target_number = @inline  x -> x - index_start
    free_target_number = @inline  x -> p.n - (x - index_start)

    # Compartment-wise rates 
    binding_rate = @inline  x -> binding_coefficient * free_target_number(x) * A * B[x]
    unbinding_rate =  @inline x -> p.k_r * bound_target_number(x) * B[x]
    division_rate = @inline x -> p.r_x[x]*B[x]*carrying_coefficient
    death_rate = @inline x -> p.d_x[x]*B[x]*carrying_coefficient
    
    dB_fun = @inline x -> binding_rate(x-1) - unbinding_rate(x) - binding_rate(x) + unbinding_rate(x+1) -
     division_rate(x)- death_rate(x)
    # For the beginning and end of the B array, we must have custom assignment in order to avoid
    # index-out-of-bounds errors
    dB[begin] =  -binding_rate(index_start) + unbinding_rate(index_start+1) -
     division_rate(index_start) - death_rate(index_start)
    dB[index_start+1:index_end-1] .= dB_fun.(index_start+1:index_end-1)
    dB[end] = binding_rate(index_end-1) -
     unbinding_rate(index_end) -
     division_rate(index_end) - death_rate(index_end)

    # The growth rates and the constant factor 2 are already baked into
    # the matrix
    # The B vector and the carrying coefficient cannot be precomputed, so
    # we have to multiply with these inside the ODE function
    # Computes the rho function for each of the components and adds it to the answer
    mul!(dB, p.f_scaled, B, carrying_coefficient, one(eltype(dB)))

    # sum() with a function inside is far more efficient than a broadcasted call over the array followed by a sum operation
    # This approach makes sure there are no or minimal allocations, considerably reducing overhead
    unbound_targets = sum(x -> free_target_number(x) / unit(eltype(B)) * B[x], eachindex(B))
    bound_targets = sum(x -> bound_target_number(x) / unit(eltype(B)) * B[x], eachindex(B))
    free_targets_released = sum(x -> p.d_x[x]*free_target_number(x) / unit(eltype(B)) * B[x], eachindex(B))
    bound_targets_released = sum(x -> p.d_x[x]*bound_target_number(x) / unit(eltype(B)) * B[x], eachindex(B))
    du.A = -binding_coefficient * (A*T + A*unbound_targets) +  
    p.k_r * (AT +  bound_targets)
    du.T = -binding_coefficient * A*T + p.k_r * AT + free_targets_released
    du.AT = binding_coefficient * A*T - p.k_r * AT + bound_targets_released
    du
end

problem = ODEProblem(ode_system!,u0,(zero(model_params.t_span),model_params.t_span),model_params)

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

sol = solve(problem,RK4(); abstol=abstol_struct, saveat=tsave)
@btime solve(problem,RK4(); abstol=abstol_struct, saveat=tsave)
@profview solve(problem,RK4();abstol=abstol_struct, saveat=tsave)

using Cthulhu
@descend solve(problem,RK4())
