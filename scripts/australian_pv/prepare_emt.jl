# Export independently checkable averaged reference points for the EMT circuit.
include("four_leg.jl")

function emt_reference(inv,voltage,current,converter_voltage,pdc;id,scope)
    report=inverter_physical_feasibility(inv;grid_voltage=voltage,grid_current=current,
        converter_voltage,converter_current=current,dc_min=0.,dc_max=5000.)
    report.checked && report.passed || error("unaccepted EMT reference: $id")
    PowerOptLab._lcl_active(inv) && error("this representative EMT matrix covers L filters only")
    pairs(v)=[[real(z),imag(z)] for z in v]
    losses=inv.p_loss_fixed+inv.a_loss*(sum(abs,current)+(inv.topology===:FOUR_LEG ? abs(sum(current)) : 0.))+
        inv.c_loss*(sum(abs2,current)+(inv.topology===:FOUR_LEG ? abs2(sum(current)) : 0.))
    Dict("case_id"=>id,"scope"=>scope,"topology"=>String(inv.topology),
        "grid_voltage"=>pairs(voltage),"converter_voltage"=>pairs(converter_voltage),"current"=>pairs(current),
        "p_dc_W"=>pdc,"p_available_W"=>5000.,"mean_converter_loss_W"=>losses,
        "phase_R_ohm"=>inv.r_filter,"phase_L_H"=>inv.x_filter/(2pi*inv.f),
        "neutral_R_ohm"=>inv.r_filter_neutral,"neutral_L_H"=>inv.x_filter_neutral/(2pi*inv.f),
        "frequency_Hz"=>inv.f,"dc_voltage_V"=>inv.v_dc,"dc_capacitance_F"=>inv.c_dc,
        "phase_current_rating_A"=>inv.i_max,"neutral_current_rating_A"=>inv.In_max,
        "capacitor_current_rating_A"=>inv.i_cap_max,"dc_2omega_peak_rating_V"=>inv.dv2_max,
        "modulation_utilization"=>inv.m_max,
        "fundamental_reference_accepted"=>true,"reference_physical_check_count"=>length(report.checks))
end

function prepare_emt(output=joinpath(ROOT,"results","emt_inputs.json"))
    request=InverterControlRequest(p_available=5000.,p_rated=HARDWARE.s_max,q_scale=HARDWARE.s_max)
    inv=plant("emt","poc";settings=TERMINAL_CONTROL_SETTINGS)
    a=cis(2pi/3);v=230 .* ComplexF64[1,a^2,a]
    cases=Dict{String,Any}[]
    for (id,voltage,q,fraction,negative) in (("balanced_L",v,0.,1.,0im),
            ("binding_dc_capacitor",v.+4.6.*[1,a,a^2],-1800.,.7,2.0+0im))
        r=allocate_inverter_services(inv,request;grid_voltage=voltage,reactive_request=q,
            active_fraction=fraction,negative_sequence_request=negative,solver_options=VA_OPTIONS)
        solve_status(r).publishable || error("EMT reference allocation failed: $id")
        p=r.allocation.plant
        push!(cases,emt_reference(inv,voltage,collect(r.allocation.grid_phase_current),
            collect(r.allocation.converter_terminal.phase_voltage),p.p_dc;id,
            scope=id=="balanced_L" ? "balanced fixed 230 V boundary, unchanged Australian hardware" :
                "synthetic fixed 2% U2 boundary, Q=-1800 var, 70% rated VW ceiling and 2 A I2 request; unchanged hardware"))
    end
    folder=joinpath(ROOT,"results","four_leg")
    env=JSON3.read(read(joinpath(folder,"environment.json"),String),Dict{String,Any})
    repo=normpath(joinpath(ROOT,"..",".."))
    for (path,digest) in env["source_sha256"]
        bytes2hex(sha256(read(joinpath(repo,path))))==digest || error("stale four-leg reference: $path")
    end
    audit=JSON3.read(read(joinpath(folder,"audit.json"),String),Dict{String,Any})
    audit["all_planned_accepted"] || error("four-leg saved-point audit required")
    case=only(c for c in four_leg_cases() if c.case_id=="vuf2")
    base=parse_bmopf(joinpath(ROOT,"networks","LV3_55bus.bmopf.json"))
    net=voltage_aware_network(base,case);spec=four_leg_fleet(net,1.)
    point=JSON3.read(read(joinpath(folder,"points","vuf2.json"),String),Dict{String,Any})
    phasors(v)=ComplexF64[complex(p[1],p[2]) for p in v]
    id=argmax(Dict(id=>abs(sum(phasors(r["converter_current"]))) for (id,r) in point["devices"]))
    c=spec.devices[id];record=point["devices"][id];ts=point["bus"][c.device.bus]
    voltage=[phasor(ts[p])-phasor(ts["n"]) for p in c.device.phase_terminals]
    pdc=audit["cases"]["vuf2"]["devices"][id]["p_dc_W"]
    push!(cases,emt_reference(c.device,voltage,phasors(record["grid_current"]),
        phasors(record["converter_voltage"]),pdc;id="four_leg_unbalanced",
        scope="saved Australian vuf2 point, greatest neutral-current customer $id, physical POC PN voltage held fixed"))
    source=study_source_sha256()
    # Reference provenance is separate from the independently executed circuit,
    # whose own source and exact input digest are recorded by emt_circuit.run.
    for file in ("four_leg.jl","voltage_aware.jl","prepare_emt.jl")
        path=joinpath(@__DIR__,file);source[relpath(path,repo)]=bytes2hex(sha256(read(path)))
    end
    mkpath(dirname(output))
    writejson(output,Dict("cases"=>cases,"source_sha256"=>source,"hardware"=>HARDWARE,
        "git_revision"=>strip(read(Cmd(["git","-C",repo,"rev-parse","HEAD"]),String)),
        "scope"=>"independent switched L-filter circuit at fixed POC voltage; not feeder EMT or service-controller dynamics",
        "dc_source"=>"constant mean input power at regulated DC-link boundary; no 2omega input power, nonzero source current ripple allowed; not a validation of the optional R-L small-signal source model",
        "converter_loss"=>"constant lumped mean loss from accepted RMS fundamental operating point; no device-level switching-loss or thermal model",
        "timing"=>Dict("switching_Hz"=>10000.,"cycles"=>10,"discard_cycles"=>6,
            "steps_per_carrier"=>[32,64,128],"pwm"=>"centered carrier PWM, sampled sinusoidal feedforward with measured DC compensation; no current/PLL/service feedback"),
        "validation_bands"=>Dict("fundamental_relative"=>.01,"fundamental_absolute_A"=>.01,
            "dc_2omega_relative"=>.02,"dc_2omega_absolute_V"=>.05,
            "mesh_current_A"=>.01,"mesh_dc_2omega_V"=>.01,"energy_relative"=>1e-6)))
    cases
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    prepare_emt(isempty(ARGS) ? joinpath(ROOT,"results","emt_inputs.json") : ARGS[1])
end
