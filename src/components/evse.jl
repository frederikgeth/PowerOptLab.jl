"""
    EV(; id, energy_max, onboard_charge_max, kwargs...)

Vehicle specification for fixed-assignment AC charging. Capacity/bounds are Wh;
`onboard_charge_max` and `onboard_discharge_max=0` are AC-referred W. One-way
`eff_charge=1`, `eff_discharge=1` include all losses between the AC port and the
stored-energy account; do not add converter losses a second time. `phase_count=1`
is the number of simultaneously used phases, not a switchable phase choice.
`q_min=q_max=0` is the onboard charger's aggregate reactive capability (var).
A specification is not a battery-state trajectory; the session supplies the
initial energy. Optional `charge_acceptance::PWLFunction` is a nonnegative, nonincreasing
AC power cap versus PE energy fraction (`:unitless` → `:W`), covering the usable
energy domain. `phase_policy=:equal_power` constrains phase sharing unless both
vehicle and EVSE explicitly select `:independent`. Heat and aging are not modeled.
"""
Base.@kwdef struct EV
    id::String
    energy_max::Float64
    onboard_charge_max::Float64
    onboard_discharge_max::Float64 = 0.0
    energy_min::Float64 = 0.0
    eff_charge::Float64 = 1.0
    eff_discharge::Float64 = 1.0
    phase_count::Int = 1
    phase_policy::Symbol = :equal_power
    charge_acceptance::Union{Nothing,PWLFunction} = nothing
    q_min::Float64 = 0.0
    q_max::Float64 = 0.0
end

"""
    EVSE(; id, bus, s_max, i_max, kwargs...)

One AC charging outlet with a fixed conductor connection. `s_max` is aggregate
VA; `i_max` is A per phase conductor. `p_charge_max=s_max` and
`p_discharge_max=0` are AC W. Default `phase_terminals=["1"]`, `neutral="n"`.
Optional `neutral_current_max` bounds return current (A). `phase_policy` is
`:equal_power` (equal phase P and Q) or explicitly `:independent`.
`q_min=q_max=0` are aggregate equipment var bounds. `available=nothing` means
always operational; otherwise supply one Bool per interval. Equipment outages
do not free an occupied outlet for another session.

This is a continuous-power scheduling envelope, not an AC/DC converter, a pilot
response model, or a DC charging cabinet. It has no independent grid injection
when empty. Shared cabinets and standby services require distinct models.
"""
Base.@kwdef struct EVSE
    id::String
    bus::String
    s_max::Float64
    i_max::Float64
    p_charge_max::Float64 = s_max
    p_discharge_max::Float64 = 0.0
    phase_terminals::Vector{String} = ["1"]
    neutral::Union{String,Nothing} = "n"
    neutral_current_max::Union{Float64,Nothing} = nothing
    phase_policy::Symbol = :equal_power
    q_min::Float64 = 0.0
    q_max::Float64 = 0.0
    available::Union{Nothing,Vector{Bool}} = nothing
end

"""
    ChargingSession(; id, ev, evse, available, energy_init, departure_energy, kwargs...)

A fixed assignment of an [`EV`](@ref) to one [`EVSE`](@ref). `available` is a
horizon-length occupancy mask. `energy_init` is stored energy (Wh) at the first
horizon boundary; it is held before arrival and after departure. A single
contiguous occupied interval is required. `departure_period=nothing` means the
last occupied interval; the target is checked at its END. A supplied departure
period must equal that interval. `departure_energy` is a stored-energy floor,
not AC energy delivered. `energy_final=nothing` optionally fixes horizon-end
energy. `allow_v2g=false` and `allow_reactive_power=false` are session permissions.

`operation=:relaxed` uses a normalized charge/discharge product bound with
`complementarity_tolerance=1e-8`; `:complementarity` requires an MPCC backend.
`:independent` is an explicit outer relaxation for research comparisons.
`acceptance_formulation=:auto` lowers concave PWL acceptance bounds exactly;
smooth formulations or `ComplementarityGraph` are explicit alternatives.
`acceptance_conservative=true` corrects smooth upper bounds using their error
contract. Endpoint bounds conservatively constrain each piecewise-constant
interval. Conservative smoothing may make zero-cap states infeasible.
Use the same EVSE object for successive sessions. One EV may have only one
session per solve; repeated journeys need a persistent mobility-state model.
"""
Base.@kwdef struct ChargingSession <: AbstractDevice
    id::String
    ev::EV
    evse::EVSE
    available::Vector{Bool}
    energy_init::Float64
    departure_energy::Float64
    departure_period::Union{Int,Nothing} = nothing
    energy_final::Union{Float64,Nothing} = nothing
    allow_v2g::Bool = false
    allow_reactive_power::Bool = false
    operation::Symbol = :relaxed
    complementarity_tolerance::Float64 = 1e-8
    acceptance_formulation::Union{Symbol,AbstractPWLFormulation} = :auto
    acceptance_conservative::Bool = true
end

function _session_device(s::ChargingSession)
    v, e = s.ev, s.evse
    operational = e.available === nothing ? s.available : s.available .& e.available
    last_occupied = findlast(s.available)
    dp = s.departure_period === nothing ? last_occupied : s.departure_period
    EVDevice(id=s.id, bus=e.bus, phase_terminals=e.phase_terminals, neutral=e.neutral,
        p_charge_max=min(v.onboard_charge_max, e.p_charge_max),
        p_discharge_max=s.allow_v2g ? min(v.onboard_discharge_max, e.p_discharge_max) : 0.0,
        energy_max=v.energy_max, energy_min=v.energy_min, energy_init=s.energy_init,
        eff_charge=v.eff_charge, eff_discharge=v.eff_discharge,
        s_max=e.s_max, i_max=e.i_max, neutral_current_max=e.neutral_current_max,
        phase_policy=v.phase_policy == :equal_power ? :equal_power : e.phase_policy,
        q_min=s.allow_reactive_power ? max(v.q_min,e.q_min) : 0.0,
        q_max=s.allow_reactive_power ? min(v.q_max,e.q_max) : 0.0,
        available=operational, departure_period=dp, departure_energy=s.departure_energy,
        energy_final=s.energy_final, operation=s.operation,
        complementarity_tolerance=s.complementarity_tolerance)
end

function validate_device(s::ChargingSession, nets; periods::Integer=length(nets))
    v, e = s.ev, s.evse
    isempty(v.id) && throw(ArgumentError("EV id must not be empty"))
    isempty(e.id) && throw(ArgumentError("EVSE id must not be empty"))
    length(s.available) == periods || throw(ArgumentError("session occupancy must match horizon"))
    e.available === nothing || length(e.available) == periods ||
        throw(ArgumentError("EVSE availability must match horizon"))
    occupied = findall(s.available)
    isempty(occupied) && throw(ArgumentError("a session needs at least one occupied interval"))
    length(occupied) == last(occupied)-first(occupied)+1 || throw(ArgumentError(
        "session occupancy must be contiguous; multiple visits require mobility-state linking"))
    s.departure_period === nothing || s.departure_period == last(occupied) ||
        throw(ArgumentError("session departure_period must equal the last occupied interval"))
    v.phase_count == length(e.phase_terminals) || throw(ArgumentError(
        "EV phase_count must match the fixed EVSE connection"))
    e.phase_policy in (:equal_power,:independent) || throw(ArgumentError("invalid EVSE phase_policy"))
    v.phase_policy in (:equal_power,:independent) || throw(ArgumentError("invalid EV phase_policy"))
    isfinite(v.energy_max) && v.energy_max > 0 || throw(ArgumentError("EV energy_max must be positive"))
    f = v.charge_acceptance
    if f !== nothing
        f.input_unit == :unitless && f.output_unit == :W || throw(ArgumentError(
            "charge_acceptance must map energy fraction (:unitless) to AC watts (:W)"))
        first(f.breakpoints) <= v.energy_min/v.energy_max && last(f.breakpoints) >= 1 ||
            throw(ArgumentError("charge_acceptance must cover the entire usable energy fraction domain"))
        all(>=(0),f.values) && all(<=(0),f.slopes) || throw(ArgumentError(
            "charge_acceptance must be nonnegative and nonincreasing"))
        s.acceptance_formulation === :auto || s.acceptance_formulation isa AbstractPWLSmoothing ||
            s.acceptance_formulation isa ComplementarityGraph || throw(ArgumentError(
                "acceptance supports :auto, smooth NLP or ComplementarityGraph only"))
        plan = plan_pwl_relation(f,(v.energy_min/v.energy_max,1.0);relation=:upper,
            formulation=s.acceptance_formulation,
            conservative=s.acceptance_conservative && s.acceptance_formulation isa AbstractPWLSmoothing)
        plan.strategy == :unresolved && throw(UnsupportedFormulation(plan.reason))
    end
    # Validate raw capabilities even when an intersection or permission would
    # otherwise hide an invalid input (e.g. a negative disabled V2G rating).
    for (label, value) in ((:onboard_charge_max,v.onboard_charge_max),
            (:onboard_discharge_max,v.onboard_discharge_max),
            (:p_charge_max,e.p_charge_max),(:p_discharge_max,e.p_discharge_max))
        isfinite(value) && value >= 0 || throw(ArgumentError("$label must be finite and nonnegative"))
    end
    for (label, value) in ((:s_max,e.s_max),(:i_max,e.i_max))
        isfinite(value) && value > 0 || throw(ArgumentError("$label must be finite and positive"))
    end
    for (lo,hi) in ((v.q_min,v.q_max),(e.q_min,e.q_max))
        isfinite(lo) && isfinite(hi) && lo <= 0 <= hi || throw(ArgumentError(
            "EV and EVSE reactive capability bounds must be finite and include zero"))
    end
    validate_device(_session_device(s), nets; periods)
end

stamp_device!(ctx, s::ChargingSession; period::Integer=1,
              active::Bool=_device_active(_session_device(s),period)) =
    stamp_device!(ctx, _session_device(s); period, active)

# Collection-level constraints cannot be validated by an isolated session.
function _validate_charging_assignments(devices)
    sessions = filter(d -> d isa ChargingSession, devices)
    vehicles = [s.ev.id for s in sessions]
    allunique(vehicles) || throw(ArgumentError(
        "one session per EV per solve is supported; do not reset a vehicle's state across sessions"))
    equipment = Dict{String,EVSE}()
    occupancy = Dict{String,Vector{Bool}}()
    for s in sessions
        e = s.evse
        if haskey(equipment,e.id)
            # Value equality on every field (including array contents).
            prior = equipment[e.id]
            all(getfield(prior,k) == getfield(e,k) for k in fieldnames(EVSE)) ||
                throw(ArgumentError("conflicting definitions for EVSE '$(e.id)'"))
            any(occupancy[e.id] .& s.available) && throw(ArgumentError(
                "overlapping sessions occupy EVSE '$(e.id)'"))
            occupancy[e.id] .|= s.available
        else
            equipment[e.id] = e
            occupancy[e.id] = copy(s.available)
        end
    end
    nothing
end
