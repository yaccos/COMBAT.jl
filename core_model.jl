using PhysicalConstants.CODATA2018
using Unitful
using Distributions


@refunit cell "cell" Cells Unitful.𝐍 false

n_targets = 100
starting_population = 1e6cell
treatment_length = 86400u"s" # 
initial_antibiotic_level = 1e6/cell
maximum_kill_rate = 0.001u"1/s"
killing_threshold = 60
replication_threshold = 50
r_max = 0.00025u"1/s"
max_kill_rate = 0.001u"1/s"
intracellular_volume = 1e-15u"L"/cell
unbinding_rate = 0.01u"1/s"
carrying_capacity = 1e9cell
molecular_weight = 555.5u"g/mol"
binding_rate = 1u"1/mol/s"
n_A = AvogadroConstant
binding_coefficient = binding_rate / (intracellular_volume * n_A)

r_fun = x -> ifelse(x <= replication_threshold, r_max - r_max*(x / replication_threshold), zero(r_max))
r_x = r_fun.(0:n_targets)

d_fun = x -> ifelse(x >= killing_threshold, max_kill_rate, zero(max_kill_rate))
d_x = d_fun.(0:n_targets)

hypergeom_density_fun = (i,x) -> pdf(Hypergeometric(x,2*n_targets - x, n_targets),i)
hypergeom_density_mat = hypergeom_density_fun.(0:n_targets,(0:n_targets)')
function hypergeom_density(i,x)
    binomial(x, i) * binomial(2*n_targets-x, n_targets-i)/ binomial(2*n_targets, n_targets)
end

# Initial conditions
B_0 = zeros(typeof(1.0cell),n_targets+1)
B_0[1] =  1.0e+6cell

A_0 = initial_antibiotic_level

AT = zero(A_0)

T_0 = zero(A_0)

function ode_system!(du, u, p, t)
    B = @view u[1:(n_targets + 1)]
    dB = @view du[1:(n_targets + 1)]
    A = u[n_targets + 2]
    dA = @view du[n_targets + 2]
    T = u[n_targets + 3]
    dT = @view du[n_targets + 3]
    AT = u[n_targets + 4]
    dAT = @view du[n_targets + 4]
    carrying_coefficient = (carrying_capacity - sum(B)) / carrying_capacity
    rho_fun = x -> 2*sum(hypergeom_density_mat[x:n_targets,x] .* B[x:n_targets] .* r_x[x:n_targets]) *
     carrying_capacity
    
    # Julia uses one number indexing,
    # while we would be would be better off using zero number indexing
    # Hence, we must be careful to get things straight
    # Therefore, we always encolse (x-1) in parenthesis to make it clear
    # that we are dealing with shifted indices even if it looks silly
    # Fortunately, this is optimized out by the compiler anyway
    dB_fun = x -> binding_coefficient * (n - (x-1) + 1) * A * B[x-1] -
     unbinding_rate * (x-1) * B[x] -
     binding_coefficient * (n_targets - (x-1)) * A * B[x] +
    unbinding_rate * ((x-1)+1) * B[x + 1] + rho_fun(x) -
     r_x[x]*B[x]*carrying_coefficient-
     d_x[x]*B[x]
    bB[1] =  -binding_coefficient * n_targets * A * B[1] +
     unbinding_rate*((1-1)+1)*B[2] +
      rho_fun(1) -
      r_x[1]*B[1]*carrying_coefficient - d_x[1]*B[1]
    dB[2:n_targets-1] = dB_fun.(2:n_targets-1)
    bB[n_targets] = binding_coefficient*A*B[n_targets-1] -
     unbinding_rate * n_targets * B[n_targets] +
     rho_fun(n_targets) -
     r_x[n_targets]*B[n_targets]*carrying_coefficient - d_x[n_targets]*B[n_targets]
    unbound_molecules = (x -> (n_targets - (x-1))*A*B[x]).(1:n_targets+1) |> sum
    bound_molecules = (x -> (x-1)*B[x]).(1:n_targets+1) |> sum
    targets_released = (x -> d_x[x]*(n_targets - (x-1))*B[x]).(1:n_targets+1) |> sum
    targets_bound = (x -> d_x[x]*(x-1)*B[x]).(1:n_targets+1) |> sum
    dA[] = -binding_coefficient * (AT + unbound_molecules) +  
    unbinding_rate * (AT +  bound_molecules)
    dT[] = -binding_coefficient * AT + unbinding_rate * AT + targets_released
    dAT[] = binding_coefficient * AT - unbinding_rate * AT + targets_bound
end
