# Fixed hardware, with the same explicitly declared operating ceilings for
# local control and coordinated dispatch. These are study assumptions, not
# manufacturer peak ratings or additional sizing optimisation variables.
include("thermal_energy.jl")
const ENGINEERING_OPERATING_LIMITS=(dc_min_V=735.,dc_max_V=765.,phase_peak_A=14.,neutral_peak_A=14.)
const ENGINEERING_OPTIONS=merge(Dict(VA_OPTIONS),Dict("mu_strategy"=>"monotone"))

function engineering_plant(base;dc_min=ENGINEERING_OPERATING_LIMITS.dc_min_V,
        dc_max=ENGINEERING_OPERATING_LIMITS.dc_max_V,dc_guard=true)
    fields=NamedTuple{fieldnames(AdvancedInverter)}(Tuple(getfield(base,k) for k in fieldnames(AdvancedInverter)))
    inv=AdvancedInverter(;merge(fields,(pwm_strategy=:CENTERED,f_sw=10000.,dc_energy_model=:nonlinear,
        pwm_current_envelope=true,pwm_dc_voltage_envelope=dc_guard,v_dc_min=dc_min,v_dc_max=dc_max,
        i_peak_max=ENGINEERING_OPERATING_LIMITS.phase_peak_A,i_grid_peak_max=ENGINEERING_OPERATING_LIMITS.phase_peak_A,
        In_peak_max=base.topology===:FOUR_LEG ? ENGINEERING_OPERATING_LIMITS.neutral_peak_A : nothing))...)
    allowed=Set((:pwm_strategy,:f_sw,:dc_energy_model,:pwm_current_envelope,:pwm_dc_voltage_envelope,
        :v_dc_min,:v_dc_max,:i_peak_max,:i_grid_peak_max,:In_peak_max))
    for key in fieldnames(AdvancedInverter)
        key in allowed || isequal(getfield(base,key),getfield(inv,key)) || error("hardware changed: $key")
    end
    validate_device(inv)
    inv
end

function engineering_fleet(net,case)
    mixed=mixed_thermal_fleet(net,case)
    devices=Dict(id=>ControlledDevice(engineering_plant(c.device),c.controller) for (id,c) in mixed.spec.devices)
    merge(mixed,(spec=ControlledInverterFleetSpec(devices,mixed.spec.requests),))
end
