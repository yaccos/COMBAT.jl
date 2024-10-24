using Distributions
using LinearAlgebra

@refunit cell "cell" Cells Unitful.𝐍 false

n_targets = 100
starting_population = 1e6cell
treatment_length = 86400.0u"s" # 
# initial_antibiotic_level = 1e6/cell
initial_antibiotic_level = 1u"mg/L"
maximum_kill_rate = 0.001u"1/s"
killing_threshold = 60
replication_threshold = 50
r_max = 0.00025u"1/s"
max_kill_rate = 0.001u"1/s"
total_volume = 1u"L"
intracellular_volume = 1e-15u"L"/cell
unbinding_rate = 0.01u"1/s"
carrying_capacity = 1e9cell
molecular_weight = 555.5u"g/mol"
binding_rate = 1u"L/mol/s"
N_A = AvogadroConstant
binding_coefficient = binding_rate / (total_volume * N_A)



r_fun = x -> ifelse(x <= replication_threshold, r_max - r_max*(x / replication_threshold), zero(r_max))
r_x = r_fun.(0:n_targets)

    d_fun = x -> ifelse(x >= killing_threshold, max_kill_rate, zero(max_kill_rate))
d_x = d_fun.(0:n_targets)

hypergeom_density_fun = (i,x) -> pdf(Hypergeometric(x,2*n_targets - x, n_targets),i)
hypergeom_density_mat = hypergeom_density_fun.(0:n_targets,(0:n_targets)')


scaled_hypergeom_density_mat = diagm(r_x) * hypergeom_density_mat
function hypergeom_density(i,x)
    binomial(x, i) * binomial(2*n_targets-x, n_targets-i)/ binomial(2*n_targets, n_targets)
end

# This is intended to serve as a preallocated pool of memory when computing the rho function
cache_rho_fun = similar(starting_population .* r_x)

# Convert from mass per volume to number of molecules in entire volume
# Also makes sure we cancel out the mass units
A_n_molecules = uconvert(NoUnits, initial_antibiotic_level / molecular_weight * total_volume * N_A)

model_params = (n=n_targets,B_0=starting_population,t_span=treatment_length,A=A_n_molecules,
D_0=maximum_kill_rate,r_T=replication_threshold,k_f=binding_rate,k_r=unbinding_rate,
W=molecular_weight,V=total_volume,C=carrying_capacity,k_T=killing_threshold,N_A=N_A,
f=hypergeom_density_mat, f_scaled=scaled_hypergeom_density_mat,
d_x = d_x, r_x=r_x, rho_cache = cache_rho_fun
)
