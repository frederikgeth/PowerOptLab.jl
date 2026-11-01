# Four-leg plant and per-phase VV/VW controller used by the mixed fleet.
include("voltage_aware.jl")

function four_leg_fleet(net,availability)
    devices=Dict{String,ControlledDevice}();requests=Dict{String,InverterControlRequest}()
    vv=PiecewiseLinearLaw(AUS_A.vv_v,AUS_A.vv_q;smoothing_epsilon=.05)
    vw=PiecewiseLinearLaw(AUS_A.vw_v,AUS_A.vw_p;smoothing_epsilon=.05)
    for (id,record) in net["ibr"]
        # A declared topology change, not a retrofit or a sizing optimization.
        record["topology"]="FOUR_LEG";record["terminal_map"]=["a","b","c","n"]
        hw=HARDWARE
        inv=AdvancedInverter(id=id,bus=record["bus"],topology=:FOUR_LEG,
            phase_terminals=["a","b","c"],neutral="n",s_max=hw.s_max,
            power_limit_location=:poc,i_max=hw.i_max,In_max=hw.i_max,
            v_dc=hw.v_dc,c_dc=hw.c_dc,dv2_max=hw.dv2_max,i_cap_max=hw.i_cap_max,
            m_max=hw.m_max,r_filter=filter_r(hw),x_filter=filter_x(hw),
            r_filter_neutral=filter_r(hw),x_filter_neutral=filter_x(hw),
            p_loss_fixed=hw.p_loss_fixed,a_loss=hw.a_loss,c_loss=hw.c_loss)
        devices[id]=ControlledDevice(inv,PerPhaseVoltVarWatt(volt_var=vv,volt_watt=vw,
            volt_watt_basis=AUS_A.vw_basis,current_target=:converter))
        requests[id]=InverterControlRequest(p_available=5000availability,p_rated=hw.s_max,q_scale=hw.s_max)
    end
    ControlledInverterFleetSpec(devices,requests)
end
