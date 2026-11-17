using PhysicalConstants.CODATA2018
using Distributions
using LinearAlgebra
using Unitful


@refunit cell "cell" Cells Unitful.𝐍 false

# Tetracycline model as described by Wiesch et al. 2015
# From describtion of models fitted from experimental data
critical_ribosome_fraction = 2.6u"percent" |> NoUnits
cell_membrane_permability = 2.0e-11u"m/s"
unspecific_association_rate = 50.0u"1/s"
unspecific_dissociation_rate = 0.99u"1/s"

# From Table S3
k_f = 3.0e5u"mol/s"
k_r = 0.5u"1/s" # The range 0.33-0.96 1/s is given in the table
n_targets = 55000/cell

# Parameter to be tuned when making the plot
initial_antibiotic_concentration = 4.0u"mg/L" # Vary from 1-24 mg/L
initial_bacterial_density = 1e6cell/u"mL" #Vary from 0-1.0e+8 cell/mL

# Times, not directly mentioned in the paper
treatment_length = 24.0u"hr" |> u"s" # we convert it right now in order to avoid type conversion when simulating
tsave = 0u"s":1u"minute":treatment_length .|> u"s" # Saves the evolution of the system for each minute

# Conditions at start
initial_internal_antibiotic_cencentraion = zero(initial_antibiotic_concentration)
AT_0 = zero(n_targets)
AU_0 = zero(initial_antibiotic_concentration)

# Constants
N_A = AvogadroConstant
molecular_weight = 444.43u"g/mol" # Specific for tetracycline

# From Table S4
intracellular_volume = 1e-15u"L"/cell # Same as in the paper by Tran et al. 2022
carrying_capacity = 1e+9cell/u"mL"

model_params = (f_c=critical_ribosome_fraction, p=cell_membrane_permability, k_ur=unspecific_association_rate, k_uf=unspecific_dissociation_rate, k_f=k_f, k_r=k_r,
 T0=n_targets, B_0=initial_bacterial_density, tspan=treatment_length, tsave=tsave, K=carrying_capacity, T_0=n_targets,
Ae_0 = initial_antibiotic_concentration, Ai_0 = initial_internal_antibiotic_cencentraion, AT_0=AT_0, AU_0=AT_0, Vi=intracellular_volume)


abstol = 1.0e-3
