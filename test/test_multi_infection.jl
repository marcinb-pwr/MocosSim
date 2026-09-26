using MocosSim: RobinForest, Event, OutsideInfectionEvent, TransmissionEvent,
  OutsideContact, HouseholdContact, ChineseStrain, DeltaStrain, OmicronStrain,
  recentbackwardinfectionof, recentforwardinfectionsof, recentstrainof,
  CrossImmunityParams, subject, saveparams
using MocosSim: SimState, Progression, Mild, Recovered, Healthy, UndefinedSeverity,
  sethealth!, setprogression!, setstrain!, finish_recovery!, health, progressionof

@testset "Multiple infection history" begin
  forest = RobinForest(2)
  first = Event(Val(OutsideInfectionEvent), 10, 1, ChineseStrain)
  secondary = Event(Val(TransmissionEvent), 12, 2, 1, HouseholdContact, ChineseStrain)
  reinfection = Event(Val(OutsideInfectionEvent), 400, 1, OmicronStrain)
  omicron_secondary = Event(Val(TransmissionEvent), 402, 2, 1, HouseholdContact, OmicronStrain)

  foreach(event -> push!(forest, event), (first, secondary, reinfection, omicron_secondary))

  @test length(forest.infections) == 4
  @test recentbackwardinfectionof(forest, 1) == reinfection
  @test recentstrainof(forest, 1) == OmicronStrain
  @test collect(recentforwardinfectionsof(forest, 1)) == [omicron_secondary]

  output = Dict{String,Any}()
  saveparams(output, forest)
  @test output["infection_subjects"] == UInt32[1, 2, 1, 2]
  @test output["strains"] == UInt8.(Int.([ChineseStrain, ChineseStrain, OmicronStrain, OmicronStrain]))
end

@testset "Recovery enables another episode" begin
  state = SimState(1)
  setprogression!(state, 1, Progression(Mild, 2, 4, missing, 8, missing))
  setstrain!(state, 1, DeltaStrain)
  sethealth!(state, 1, Recovered)
  finish_recovery!(state, 1)

  @test health(state, 1) == Healthy
  @test state.progressions[1].severity == UndefinedSeverity
  @test state.individuals[1].strain == MocosSim.NullStrain
end

@testset "Cross-immunity configuration" begin
  weights = zeros(4, 4)
  weights[Int(DeltaStrain), Int(OmicronStrain)] = 0.2
  params = CrossImmunityParams(weights, [0.9, 0.8, 0.7, 0.3], 180)
  @test params.infection[Int(DeltaStrain), Int(OmicronStrain)] == 0.2
  @test params.vaccination[Int(OmicronStrain)] == 0.3
  @test params.half_life == 180
  @test_throws ArgumentError CrossImmunityParams(zeros(3, 4), ones(4), Inf)
  @test_throws ArgumentError CrossImmunityParams(fill(1.1, 4, 4), ones(4), Inf)
end
