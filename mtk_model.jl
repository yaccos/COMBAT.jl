using Unitful
using ModelingToolkit
using DifferentialEquations
using DynamicQuantities

include("params.jl")

strip_units(x) = x isa Unitful.AbstractQuantity ? Float64(Unitful.ustrip(x)) : Float64(x)
strip_units_array(x) = strip_units.(x)

function build_mtk_problem(params; check_units = true)
    n = Int(params.n)
    nB = n + 1

    A0 = strip_units(params.A)
    B0 = strip_units(params.B_0)
    C = strip_units(params.C)
    V = strip_units(params.V)
    N_A = strip_units(params.N_A)
    k_f = strip_units(params.k_f)
    k_r = strip_units(params.k_r)
    beta = k_f / (V * N_A)
    r_x = strip_units_array(params.r_x)
    d_x = strip_units_array(params.d_x)
    f_dimless = strip_units_array(params.f)
    t_end = strip_units(params.t_span)

    @independent_variables t [unit = DynamicQuantities.u"s"]
    @variables A(t) [unit = DynamicQuantities.u"1"]
    @variables T(t) [unit = DynamicQuantities.u"1"]
    @variables AT(t) [unit = DynamicQuantities.u"1"]
    @variables B(t)[1:nB] [unit = DynamicQuantities.u"1"]
    @parameters beta_p [unit = DynamicQuantities.u"s^-1"]
    @parameters k_r_p [unit = DynamicQuantities.u"s^-1"]
    @parameters C_p [unit = DynamicQuantities.u"1"]
    @parameters rate_p [unit = DynamicQuantities.u"s^-1"]
    D = Differential(t)

    carrying = (C_p - sum(B[i] for i in 1:nB)) / C_p
    rho = [rate_p * carrying * sum(2 * f_dimless[i, j] * r_x[j] * B[j] for j in 1:nB) for i in 1:nB]

    eqs = Equation[]

    for i in 1:nB
        free_targets = n - (i - 1)
        bound_targets = i - 1

        division_rate = rate_p * r_x[i] * B[i] * carrying
        death_rate = rate_p * d_x[i] * B[i] * carrying

        if i == 1
            binding_here = beta_p * free_targets * A * B[i]
            unbinding_next = k_r_p * (bound_targets + 1) * B[i + 1]
            dBi = -binding_here + unbinding_next - division_rate - death_rate + rho[i]
        elseif i == nB
            binding_prev = beta_p * (free_targets + 1) * A * B[i - 1]
            unbinding_here = k_r_p * bound_targets * B[i]
            dBi = binding_prev - unbinding_here - division_rate - death_rate + rho[i]
        else
            binding_prev = beta_p * (free_targets + 1) * A * B[i - 1]
            unbinding_here = k_r_p * bound_targets * B[i]
            binding_here = beta_p * free_targets * A * B[i]
            unbinding_next = k_r_p * (bound_targets + 1) * B[i + 1]
            dBi = binding_prev - unbinding_here - binding_here + unbinding_next - division_rate - death_rate + rho[i]
        end

        push!(eqs, D(B[i]) ~ dBi)
    end

    unbound_targets = sum((n - (i - 1)) * B[i] for i in 1:nB)
    bound_targets = sum((i - 1) * B[i] for i in 1:nB)
    free_targets_released = sum(rate_p * d_x[i] * (n - (i - 1)) * B[i] for i in 1:nB)
    bound_targets_released = sum(rate_p * d_x[i] * (i - 1) * B[i] for i in 1:nB)

    push!(eqs, D(A) ~ -beta_p * (A * T + A * unbound_targets) + k_r_p * (AT + bound_targets))
    push!(eqs, D(T) ~ -beta_p * A * T + k_r_p * AT + free_targets_released)
    push!(eqs, D(AT) ~ beta_p * A * T - k_r_p * AT + bound_targets_released)

    if check_units && !ModelingToolkit.validate(eqs)
        error("ModelingToolkit unit validation failed. Inspect equation warnings above for offending terms.")
    end

    states = [A, T, AT, B...]
    params_sym = [beta_p, k_r_p, C_p, rate_p]
    @named combat_sys = ODESystem(eqs, t, states, params_sym)
    combat_sys_compiled = mtkcompile(combat_sys)

    u0 = [A => A0, T => 0.0, AT => 0.0]
    append!(u0, [B[i] => (i == 1 ? B0 : 0.0) for i in 1:nB])
    pmap = [beta_p => beta, k_r_p => k_r, C_p => C, rate_p => 1.0]

    prob = ODEProblem(combat_sys_compiled, u0, (0.0, t_end), pmap)
    return prob, combat_sys_compiled
end

function solve_mtk_model(params; solver = RK4(), abstol = abstol, reltol = 1e-6, saveat = strip_units_array(tsave), check_units = true)
    prob, sys = build_mtk_problem(params; check_units = check_units)
    sol = solve(prob, solver; abstol = abstol, reltol = reltol, saveat = saveat)
    return sol, prob, sys
end

@time sol_mtk, prob_mtk, sys_mtk = solve_mtk_model(model_params)
