import Base: iterate, length, eltype, sizehint!, push!

const InfectionIdx = UInt32
const EdgeIdx = UInt32

struct EdgeDictKey
  infection_id::InfectionIdx
  edge_id::EdgeIdx
end

"""A transmission forest indexed by infection episode rather than by person.

One person may occur in `infections` many times. `recent` makes the common
current-episode lookup O(1), while the append-only arrays retain the complete
history needed to analyse reinfections.
"""
struct RobinForest
  infections::Vector{Event}
  outdegrees::Vector{EdgeIdx}
  outedgedict::RobinDict{EdgeDictKey,InfectionIdx}
  recent::Vector{InfectionIdx}
end

RobinForest(num_individuals::Integer) = RobinForest(
  Event[], EdgeIdx[], RobinDict{EdgeDictKey,InfectionIdx}(),
  fill(InfectionIdx(0), num_individuals)
)

function sizehint!(forest::RobinForest, num_infections::Integer)
  sizehint!(forest.infections, num_infections)
  sizehint!(forest.outdegrees, num_infections)
  sizehint!(forest.outedgedict, num_infections)
  forest
end

function reset!(forest::RobinForest)
  empty!(forest.infections)
  empty!(forest.outdegrees)
  empty!(forest.outedgedict)
  fill!(forest.recent, InfectionIdx(0))
  nothing
end

function push!(forest::RobinForest, infection::Event)
  push!(forest.infections, infection)
  push!(forest.outdegrees, EdgeIdx(0))
  infection_id = InfectionIdx(length(forest.infections))
  subject_id = subject(infection)
  forest.recent[subject_id] = infection_id

  source_id = source(infection)
  if source_id == 0
    @assert contactkind(infection) == OutsideContact
    return nothing
  end

  source_infection_id = forest.recent[source_id]
  @assert source_infection_id != 0 "source $source_id has no registered infection"
  source_outdegree = forest.outdegrees[source_infection_id] += EdgeIdx(1)
  key = EdgeDictKey(source_infection_id, source_outdegree)
  @assert !haskey(forest.outedgedict, key)
  forest.outedgedict[key] = infection_id
  nothing
end

recentinfectionof(forest::RobinForest, person_id::Integer)::InfectionIdx = forest.recent[person_id]
function recentbackwardinfectionof(forest::RobinForest, person_id::Integer)::Event
  infection_id = recentinfectionof(forest, person_id)
  infection_id == 0 ? Event() : forest.infections[infection_id]
end
recentforwardinfectionsof(forest::RobinForest, person_id::Integer) = OutEdgeList(forest, person_id)
recentstrainof(forest::RobinForest, person_id::Integer)::StrainKind = strainkind(recentbackwardinfectionof(forest, person_id))

# Compatibility names now deliberately mean the most recent infection episode.
backwardinfection(forest::RobinForest, person_id::Integer) = recentbackwardinfectionof(forest, person_id)
forwardinfections(forest::RobinForest, person_id::Integer) = recentforwardinfectionsof(forest, person_id)
strainof(forest::RobinForest, person_id::Integer) = recentstrainof(forest, person_id)

struct OutEdgeList
  forest::RobinForest
  infection_id::InfectionIdx
  count::EdgeIdx
end

function OutEdgeList(forest::RobinForest, person_id::Integer)
  infection_id = recentinfectionof(forest, person_id)
  count = infection_id == 0 ? EdgeIdx(0) : forest.outdegrees[infection_id]
  OutEdgeList(forest, infection_id, count)
end

iterate(oel::OutEdgeList, state=1) = state > oel.count ? nothing :
  (oel.forest.infections[oel.forest.outedgedict[EdgeDictKey(oel.infection_id, state)]], state + 1)
eltype(::Type{OutEdgeList}) = Event
length(oel::OutEdgeList) = Int(oel.count)

function saveparams(dict, forest::RobinForest, prefix::AbstractString="")
  N = length(forest.infections)
  infection_times = Vector{Float32}(undef, N)
  infection_subjects = Vector{UInt32}(undef, N)
  infection_sources = Vector{UInt32}(undef, N)
  contact_kinds = Vector{UInt8}(undef, N)
  strains = Vector{UInt8}(undef, N)
  @simd for i in 1:N
    event = forest.infections[i]
    infection_times[i] = Float32(time(event))
    infection_subjects[i] = subject(event)
    infection_sources[i] = source(event)
    contact_kinds[i] = UInt8(contactkind(event))
    strains[i] = UInt8(strainkind(event))
  end
  dict[prefix*"infection_times"] = infection_times
  dict[prefix*"infection_subjects"] = infection_subjects
  dict[prefix*"infection_sources"] = infection_sources
  dict[prefix*"contact_kinds"] = contact_kinds
  dict[prefix*"strains"] = strains
end
