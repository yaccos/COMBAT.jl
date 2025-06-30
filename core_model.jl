using PhysicalConstants.CODATA2018
using Unitful
using DifferentialEquations
using BenchmarkTools
include("array_interface.jl")
include("params.jl")

# Initial conditions

function initialize_system(params)
    A_0 = params.A
    T_0 = zero(A_0)
    AT_0 = zero(A_0)
    B_0 = params.B_0
    B_start = zeros(typeof(B_0),params.n+1)
    B_start[1] = B_0
    DiscreteSimulationVariables(A_0,T_0,AT_0,B_start)
end

u0 = initialize_system(model_params)

abstol_struct = similar(u0)

abstol_struct.A = abstol*oneunit(abstol_struct.A)
abstol_struct.T = abstol*oneunit(abstol_struct.T)
abstol_struct.AT = abstol*oneunit(abstol_struct.AT)
fill!(abstol_struct.B,abstol*oneunit(eltype(abstol_struct.B)))

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
    @inline bound_target_number = x -> x - index_start
    @inline free_target_number = x -> p.n - (x - index_start)

    

    # Compartment-wise rates 
    @inline binding_rate = x -> binding_coefficient * free_target_number(x) * A * B[x]
    @inline unbinding_rate = x -> p.k_r * bound_target_number(x) * B[x]
    @inline division_rate = x -> p.r_x[x]*B[x]*carrying_coefficient
    @inline death_rate = x -> p.d_x[x]*B[x]*carrying_coefficient
    
    @inline dB_fun = x -> binding_rate(x-1) - unbinding_rate(x) - binding_rate(x) + unbinding_rate(x+1) -
     division_rate(x)- death_rate(x)
    # For the beginning and end of the B array, we must have custom assignment in order to avoid
    # index-out-of-bounds errors
    dB[begin] =  -binding_rate(index_start) + unbinding_rate(index_start+1) -
     division_rate(index_start) - death_rate(index_start)
    dB[begin+1:end-1] .= dB_fun.(index_start+1:index_end-1)
    dB[end] = binding_rate(index_end-1) -
     unbinding_rate(index_end) -
     division_rate(index_end) - death_rate(index_end)

    # The growth rates and the constant factor 2 are already baked into
    # the matrix
    # The B vector and the carrying coefficient cannot be precomputed, so
    # we have to multiply with these inside the ODE function
    #Computes the rho function for each of the components and add it to the answer
    mul!(dB, p.f_scaled, B, carrying_coefficient, one(eltype(dB)))

    # mapreduce is far more efficient than a broadcasted call over the array followed by a sum operation
    # This approach makes sure there are no or minimal allocations, considerably reducing overhead
    unbound_targets = sum(x -> free_target_number(x) / oneunit(eltype(B)) * B[x], eachindex(B))
    bound_targets = sum(x -> bound_target_number(x) / oneunit(eltype(B)) * B[x], eachindex(B))
    free_targets_released = sum(x -> p.d_x[x]*free_target_number(x) / oneunit(eltype(B)) * B[x], eachindex(B))
    bound_targets_released = sum(x -> p.d_x[x]*bound_target_number(x) / oneunit(eltype(B)) * B[x], eachindex(B))
    du.A = -binding_coefficient * (A*T + A*unbound_targets) +  
    p.k_r * (AT +  bound_targets)
    du.T = -binding_coefficient * A*T + p.k_r * AT + free_targets_released
    du.AT = binding_coefficient * A*T - p.k_r * AT + bound_targets_released
    
    du
end

problem = ODEProblem(ode_system!,u0,(zero(model_params.t_span),model_params.t_span),model_params)

# sol = solve(problem,AutoTsit5(Rosenbrock23()))
sol = solve(problem,RK4(); abstol=abstol_struct, saveat=tsave)
@btime solve(problem,RK4(); abstol=abstol_struct, saveat=tsave)
@profview solve(problem,RK4();abstol=abstol_struct, saveat=tsave)

using Cthulhu
@descend solve(problem,RK4())
