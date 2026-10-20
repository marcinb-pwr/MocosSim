using DataFrames

@testset "Infection while quarantined preserves observation" begin
  population = DataFrame(household_index=[1,1], attending_school=[1,0],
    school_index=[1,0], class_index=[1,0], age=[10,40], gender=[0,1],
    ishealthcare=[false,false])
  params = MocosSim.make_params(MersenneTwister(1); individuals_df=population,
    progression_params=MocosSim.make_progression_params(1.0,1.0,1.0),
    household_kernel_param=0.0, constant_kernel_param=0.0)
  state = MocosSim.SimState(2)
  push!(state.forest, MocosSim.Event(Val(MocosSim.OutsideInfectionEvent), 0, 1,
    MocosSim.ChineseStrain))
  MocosSim.sethealth!(state, 1, MocosSim.Infectious)
  # The second person is quarantined but can still become infected in the household.
  first_quarantine = MocosSim.Event(Val(MocosSim.QuarantinedEvent), 0, 2, false)
  @test MocosSim.execute!(state, params, first_quarantine)
  @test MocosSim.detected(state, 2) == MocosSim.UnderObservation
  infection = MocosSim.Event(Val(MocosSim.TransmissionEvent), 1, 2, 1,
    MocosSim.HouseholdContact, MocosSim.ChineseStrain)
  @test MocosSim.execute!(state, params, infection)
  @test MocosSim.detected(state, 2) == MocosSim.UnderObservation
  @test MocosSim.isquarantined(state, 2)
  @test MocosSim.freedom(state, 2) == MocosSim.HomeQuarantine
  @test !MocosSim.execute!(state, params,
    MocosSim.Event(Val(MocosSim.QuarantinedEvent), 2, 2, false))
  @test MocosSim.quarantine_level(state, 2) == 2
  @test !MocosSim.execute!(state, params,
    MocosSim.Event(Val(MocosSim.QuarantineEndEvent), 14, 2))
  @test MocosSim.execute!(state, params,
    MocosSim.Event(Val(MocosSim.QuarantineEndEvent), 16, 2))
  @test !MocosSim.isquarantined(state, 2)
  @test MocosSim.freedom(state, 2) == MocosSim.Free
end

@testset "Forward tracing accepts a recovered susceptible contact" begin
  population = DataFrame(household_index=[1,1], attending_school=[1,0],
    school_index=[1,0], class_index=[1,0], age=[10,40], gender=[0,1],
    ishealthcare=[false,false])
  params = MocosSim.make_params(MersenneTwister(1); individuals_df=population,
    progression_params=MocosSim.make_progression_params(1.0,1.0,1.0),
    household_kernel_param=0.0, constant_kernel_param=0.0, forward_tracing_prob=1.0)
  state = MocosSim.SimState(2)
  push!(state.forest, MocosSim.Event(Val(MocosSim.OutsideInfectionEvent), 0, 1,
    MocosSim.ChineseStrain))
  push!(state.forest, MocosSim.Event(Val(MocosSim.TransmissionEvent), 1, 2, 1,
    MocosSim.HouseholdContact, MocosSim.ChineseStrain))
  MocosSim.setprogression!(state, 2, MocosSim.Progression(MocosSim.Mild, 2, 4, missing, 8, missing))
  MocosSim.sethealth!(state, 2, MocosSim.Recovered)
  MocosSim.finish_recovery!(state, 2)
  @test MocosSim.health(state, 2) == MocosSim.Healthy
  @test isnothing(MocosSim.forwardtrace!(state, params, 1; trace_household_connections=true))
  @test !isempty(state.queue)
  traced = pop!(state.queue)
  @test MocosSim.kind(traced) == MocosSim.TracedEvent
  @test MocosSim.subject(traced) == 2
end

@testset "Detection reset is specific to the new infection episode" begin
  state = MocosSim.SimState(1)
  for previous in (MocosSim.UnderObservation, MocosSim.TestPending, MocosSim.Detected)
    MocosSim.reset!(state)
    MocosSim.setdetected!(state, 1, previous)
    MocosSim.resetdetected!(state, 1)
    @test MocosSim.detected(state, 1) == MocosSim.Undetected
    MocosSim.setdetected!(state, 1, previous)
    MocosSim.quarantine_advance!(state, 1, 1)
    MocosSim.setfreedom!(state, 1, MocosSim.HomeQuarantine)
    MocosSim.resetdetected!(state, 1)
    @test MocosSim.detected(state, 1) == MocosSim.UnderObservation
    @test MocosSim.quarantine_level(state, 1) == 1
  end
end
