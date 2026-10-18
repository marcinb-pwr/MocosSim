
const StrainInfectivityTable = SVector{NUM_STRAINS, Float64}

"""Literature-calibrated protection against reinfection at immune peak.

Rows are the infecting strain and columns the challenge strain in the order
ancestral, Alpha, Delta, Omicron BA.1/BA.2. These are population-level
infection-protection proxies, not protection against severe disease. Published
studies do not all use the same (any versus symptomatic infection) endpoint.
The pre-Omicron values are anchored to the >82% pooled protection in
Lancet 2023 (doi:10.1016/S0140-6736(22)02465-5); the Omicron column is anchored
to its ≈45% pooled estimate and Qatar/England variant-specific studies
(doi:10.1056/NEJMoa2203965, doi:10.1016/S1473-3099(22)00729-0).
"""
function calibrated_cross_immunity_table()
  [0.90 0.88 0.85 0.45;
   0.88 0.90 0.86 0.45;
   0.82 0.84 0.90 0.50;
   0.40 0.42 0.45 0.60]
end

"""Heuristic uncertainty scale for scenario-level sensitivity draws.

This is intentionally broader for sparsely observed reverse-chronology cells;
it is not presented as a pooled sampling standard error.
"""
function calibrated_cross_immunity_std()
  [0.04 0.04 0.05 0.08;
   0.04 0.04 0.05 0.08;
   0.05 0.05 0.04 0.08;
   0.10 0.10 0.10 0.08]
end

"""Recent mRNA-booster protection proxy by challenge strain.

Delta and Omicron values follow the England test-negative study
(doi:10.1056/NEJMoa2119451). A vaccination programme with a different product,
dose count, or time-since-dose should override these values.
"""
calibrated_vaccination_immunity() = [0.90, 0.90, 0.95, 0.67]

function _uncertainty_shape(std, dims, label)
  values = std isa Real ? fill(Float64(std), dims) : Float64.(std)
  size(values) == dims || throw(ArgumentError("$label must be a scalar or have size $dims"))
  all(values .>= 0) || throw(ArgumentError("$label must be non-negative"))
  values
end

"""Draw bounded scenario probabilities around a mean table.

The draw is made once while constructing `SimParams`; it represents epistemic
calibration uncertainty and therefore must not be redrawn per infection.
"""
function sample_immunity_calibration(rng::AbstractRNG, means, std; label="immunity_std")
  means_array = Float64.(means)
  deviations = _uncertainty_shape(std, size(means_array), label)
  all(0 .<= means_array .<= 1) || throw(ArgumentError("immunity means must be in [0, 1]"))
  map(means_array, deviations) do mean, deviation
    deviation == 0 ? mean : rand(rng, truncated(Normal(mean, deviation), 0.0, 1.0))
  end
end

"""Protection probabilities against each challenge strain.

Rows of `infection` are the strain of the most recent prior infection and
columns are the challenge strain. `vaccination` applies while the existing
infection-immunity vaccination flag is active. Probabilities may optionally
wane exponentially; `Inf` keeps them constant. Keyword defaults use the
literature calibration above; explicit positional values remain available for
scenario-specific calibration.
"""
struct CrossImmunityParams
  infection::Matrix{Float64}
  vaccination::Vector{Float64}
  infection_half_life::Float64
  vaccination_half_life::Float64
  function CrossImmunityParams(infection, vaccination, infection_half_life::Real, vaccination_half_life::Real)
    size(infection) == (NUM_STRAINS, NUM_STRAINS) || throw(ArgumentError("infection cross-immunity must be a $NUM_STRAINS×$NUM_STRAINS matrix"))
    length(vaccination) == NUM_STRAINS || throw(ArgumentError("vaccination protection must contain $NUM_STRAINS values"))
    all(0 .<= infection .<= 1) || throw(ArgumentError("cross-immunity probabilities must be in [0, 1]"))
    all(0 .<= vaccination .<= 1) || throw(ArgumentError("vaccination probabilities must be in [0, 1]"))
    (infection_half_life > 0 || infection_half_life == Inf) || throw(ArgumentError("infection-immunity half-life must be positive"))
    (vaccination_half_life > 0 || vaccination_half_life == Inf) || throw(ArgumentError("vaccination-immunity half-life must be positive"))
    new(Matrix{Float64}(infection), Vector{Float64}(vaccination),
      Float64(infection_half_life), Float64(vaccination_half_life))
  end
end

CrossImmunityParams(; infection=calibrated_cross_immunity_table(),
  vaccination=calibrated_vaccination_immunity(), infection_half_life=500.0,
  vaccination_half_life=100.0) = CrossImmunityParams(
    infection, vaccination, infection_half_life, vaccination_half_life)

# Compatibility with the initial cross-immunity API: the third positional
# argument is the post-infection half-life. Vaccine protection remains constant
# until its scheduled loss event unless a fourth argument is supplied.
CrossImmunityParams(infection, vaccination, infection_half_life::Real) =
  CrossImmunityParams(infection, vaccination, infection_half_life, Inf)

wanedprotection(base::Real, age::Real, half_life::Real) =
  half_life == Inf ? Float64(base) : Float64(base) * exp2(-max(0.0, Float64(age)) / half_life)

function make_infectivity_table(;base_multiplier::Real=1.0, british_multiplier::Real=1.70, delta_multiplier::Real=1.7*1.5, omicron_multiplier::Real=1.7*1.5*2.0)::StrainInfectivityTable
  # needs validation with real data
  immunity = StrainInfectivityTable(base_multiplier, british_multiplier, delta_multiplier, omicron_multiplier)
  @assert all( 0 .<= immunity )
  immunity
end

struct ImmunizationEvents
  subjects::Vector{PersonIdx}
  immunity_uptake_times::Vector{TimePoint}
  lost_immunity_times::Vector{TimePoint}
  immunity_types::Vector{ImmunityState}
end

function make_immunity_table(state::AbstractSimState, level::Real)
  N = numindividuals(state)
  subjects = PersonIdx[]
  immunity_uptake_times = TimePoint[]
  lost_immunity_times = TimePoint[]
  immunity_types = ImmunityState[]
  for id in 1:floor(N*level)
    time_begin = 0
    time_half = rand(state.rng)*90
    time_end = time_half + rand(state.rng)*90
    push!(subjects, id)
    push!(immunity_uptake_times, time_begin)
    push!(lost_immunity_times, time_half)
    push!(immunity_types, against_infection)
    push!(subjects, id)
    push!(immunity_uptake_times, time_half)
    push!(lost_immunity_times, time_end)
    push!(immunity_types, against_severe_progression)
  end
  ImmunizationEvents(subjects, immunity_uptake_times, lost_immunity_times, immunity_types)
end

straininfectivity(table::StrainInfectivityTable, strain::StrainKind) = table[UInt(strain)]


function immunize!(state::SimState, immunization::ImmunizationEvents)::Nothing
  @info "Immunizing"
  N = length(immunization.immunity_uptake_times)
  current_time = time(state)
  count = 0
  enqueued = 0
  for i in 1:N
    if immunization.immunity_uptake_times[i] <= current_time
      if immunization.lost_immunity_times[i] <= current_time
        continue
      end
      setimmunity!(state, immunization.subjects[i], immunization.immunity_types[i])
      push!(state.queue,
        Event( Val(ImmunizationEvent),
          immunization.lost_immunity_times[i],
          immunization.subjects[i],
          immunization.immunity_types[i]
        )
      )
      count += 1
    else
      push!(state.queue,
        Event( Val(ImmunizationEvent),
          immunization.immunity_uptake_times[i],
          immunization.subjects[i],
          immunization.immunity_types[i]
        )
      )
      push!(state.queue,
        Event( Val(ImmunizationEvent),
          immunization.lost_immunity_times[i],
          immunization.subjects[i],
          immunization.immunity_types[i]
        )
      )
      enqueued += 1
    end
  end

  @info "Executed $N entires: immunized $count individuals and enqueued $enqueued "


  nothing
end
