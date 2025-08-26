using Unitful
using DifferentialEquations
f(u, p, t) = u ./ u"s"
f!(du,u, p, t) = (du .= u ./ u"s")
u0 = [1.0]
t_span = (0.0u"s", 1.0u"s")
problem_const = ODEProblem(f, u0, t_span)
problem_mut = ODEProblem(f!, u0, t_span)
solve(problem_const, ESERK5())
solve(problem_mut, ESERK5())




sol = solve(problem_mut, RK4())
sol = solve(problem_mut, Tsit5())
sol = solve(problem_const, Vern6())
sol = solve(problem_const, Rosenbrock23(autodiff=AutoFiniteDiff()))
