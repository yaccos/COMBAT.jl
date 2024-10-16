include("core_model.jl")
d_u0 = similar(u0) ./ 1u"s"

ode_system!(d_u0, u0, model_params, 0u"s")
