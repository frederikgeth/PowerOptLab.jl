# Export averaged reference points for the separately coded EMT circuit.
include("four_leg.jl")

function emt_reference(inv,voltage,current,converter_voltage,pdc;id,scope)
    report=inverter_physical_feasibility(inv;grid_voltage=voltage,grid_current=current,
        converter_voltage,converter_current=current,dc_min=0.,dc_max=5000.,
        pwm_thermal_headroom=inv.pwm_strategy===:NONE ? nothing : PowerOptLab._PWM_ALLOCATION_HEADROOM)
    report.checked && report.passed || error("unaccepted EMT reference: $id")
    PowerOptLab._lcl_active(inv) && error("this representative EMT matrix covers L filters only")
    pairs(v)=[[real(z),imag(z)] for z in v]
    losses=inv.p_loss_fixed+inv.a_loss*(sum(abs,current)+(inv.topology===:FOUR_LEG ? abs(sum(current)) : 0.))+
        inv.c_loss*(sum(abs2,current)+(inv.topology===:FOUR_LEG ? abs2(sum(current)) : 0.))
    Dict("case_id"=>id,"scope"=>scope,"topology"=>String(inv.topology),
        "dc_energy_model"=>String(inv.dc_energy_model),"reference_pwm_strategy"=>String(inv.pwm_strategy),
        "pwm_current_envelope"=>inv.pwm_current_envelope,"pwm_dc_voltage_envelope"=>inv.pwm_dc_voltage_envelope,
        "grid_voltage"=>pairs(voltage),"converter_voltage"=>pairs(converter_voltage),"current"=>pairs(current),
        "p_dc_W"=>pdc,"p_available_W"=>5000.,"mean_converter_loss_W"=>losses,
        "phase_R_ohm"=>inv.r_filter,"phase_L_H"=>inv.x_filter/(2pi*inv.f),
        "neutral_R_ohm"=>inv.r_filter_neutral,"neutral_L_H"=>inv.x_filter_neutral/(2pi*inv.f),
        "frequency_Hz"=>inv.f,"dc_voltage_V"=>inv.v_dc,"dc_capacitance_F"=>inv.c_dc,
        "phase_current_rating_A"=>inv.i_max,"neutral_current_rating_A"=>inv.In_max,
        "phase_peak_rating_A"=>inv.i_peak_max,"grid_peak_rating_A"=>inv.i_grid_peak_max,
        "neutral_peak_rating_A"=>inv.In_peak_max,"grid_current_rating_A"=>inv.i_grid_max,
        "dc_voltage_min_rating_V"=>inv.v_dc_min,"dc_voltage_max_rating_V"=>inv.v_dc_max,
        "dc_switching_voltage_bound_V"=>inv.pwm_dc_voltage_envelope ?
            PowerOptLab._dc_switching_voltage_bound(inv,real.(current),imag.(current)).total : nothing,
        "capacitor_current_rating_A"=>inv.i_cap_max,"dc_2omega_peak_rating_V"=>inv.dv2_max,
        "modulation_utilization"=>inv.m_max,
        "fundamental_reference_accepted"=>true,"reference_physical_check_count"=>length(report.checks))
end
