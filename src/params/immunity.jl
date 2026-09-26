
const StrainInfectivityTable = SVector{NUM_STRAINS, Float64}

"""Protection probabilities against each challenge strain.

Rows of `infection` are the strain of the most recent prior infection and
columns are the challenge strain. `vaccination` applies while the existing
infection-immunity vaccination flag is active. Probabilities may optionally
wane exponentially; `Inf` keeps them constant. Defaults preserve the old
behaviour (vaccination blocks infection, prior infection adds no protection),
so calibrated values must be supplied explicitly for an Omicron scenario.
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

CrossImmunityParams(; infection=zeros(NUM_STRAINS, NUM_STRAINS),
  vaccination=ones(NUM_STRAINS), infection_half_life=Inf,
  vaccination_half_life=Inf) = CrossImmunityParams(
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
