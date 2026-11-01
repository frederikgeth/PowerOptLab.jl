# Fixed hardware, with the same explicitly declared operating ceilings for
# local control and coordinated dispatch. These are study assumptions, not
# manufacturer peak ratings or additional sizing optimisation variables.
include("thermal_energy.jl")
# Declared before the corrected campaign (design rule in rated_capability.jl):
# the fixture must deliver its 6 kVA nameplate at the AS/NZS 4777.2 Australia A
# volt-var corners (207-258 V, -60 %/+44 % of rated S) under the model's own
# conservative carrier bounds, with phase current as the limiting quantity.
# Model requirement at those corners: 6.22 A capacitor RMS and 20.1 A
# phase/neutral peak; +10 % margin, rounded up to whole amperes. The rated
# unity-PF capacitor current cross-checks against Kolar & Round (2006):
# about 0.59 x 9.2 A = 5.4 A without filter ripple, versus 5.7 A modelled.
const ENGINEERING_OPERATING_LIMITS=(dc_min_V=735.,dc_max_V=765.,phase_peak_A=23.,neutral_peak_A=23.,
    capacitor_rms_A=7.)
# The earlier fixture: 1 A capacitor budget (chosen when the model only had a
# 2-omega term) and 14 A peak (about sqrt(2) x 10 A, no ripple headroom). It
# limits every inverter to about 1.2 kVA; retained only as a labelled sensitivity.
const ENGINEERING_LEGACY_LIMITS=(dc_min_V=735.,dc_max_V=765.,phase_peak_A=14.,neutral_peak_A=14.,
    capacitor_rms_A=1.)
const ENGINEERING_OPTIONS=merge(Dict(VA_OPTIONS),Dict("mu_strategy"=>"monotone"))

function engineering_plant(base;limits=ENGINEERING_OPERATING_LIMITS,dc_min=limits.dc_min_V,
        dc_max=limits.dc_max_V,dc_guard=true)
    fields=NamedTuple{fieldnames(AdvancedInverter)}(Tuple(getfield(base,k) for k in fieldnames(AdvancedInverter)))
    inv=AdvancedInverter(;merge(fields,(pwm_strategy=:CENTERED,f_sw=10000.,dc_energy_model=:nonlinear,
        pwm_current_envelope=true,pwm_dc_voltage_envelope=dc_guard,v_dc_min=dc_min,v_dc_max=dc_max,
        i_cap_max=limits.capacitor_rms_A,
        i_peak_max=limits.phase_peak_A,i_grid_peak_max=limits.phase_peak_A,
        In_peak_max=base.topology===:FOUR_LEG ? limits.neutral_peak_A : nothing))...)
    allowed=Set((:pwm_strategy,:f_sw,:dc_energy_model,:pwm_current_envelope,:pwm_dc_voltage_envelope,
        :v_dc_min,:v_dc_max,:i_cap_max,:i_peak_max,:i_grid_peak_max,:In_peak_max))
    for key in fieldnames(AdvancedInverter)
        key in allowed || isequal(getfield(base,key),getfield(inv,key)) || error("hardware changed: $key")
    end
    validate_device(inv)
    inv
end

function engineering_fleet(net,case;limits=ENGINEERING_OPERATING_LIMITS)
    mixed=mixed_thermal_fleet(net,case)
    devices=Dict(id=>ControlledDevice(engineering_plant(c.device;limits),c.controller) for (id,c) in mixed.spec.devices)
    merge(mixed,(spec=ControlledInverterFleetSpec(devices,mixed.spec.requests),))
end
