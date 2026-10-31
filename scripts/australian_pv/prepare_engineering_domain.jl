# A predeclared carrier-domain experiment. Engineering-limit failures are kept;
# passing the carrier approximation checks does not make those points feasible.
include("switching_engineering.jl")

const ENGINEERING_DOMAIN=(p_W=(0.,600.,5000.),q_var=(-2000.,0.,2000.),
    negative_sequence_A=(0.,1.),sequence_angle_rad=(0.,2pi/3,4pi/3),
    modulation_fraction=(.3,.7,.98),carrier_line_angle_rad=(0.,pi/3,2pi/3),
    voltage_negative_sequence_fraction=.02,four_leg_zero_sequence_A=.1,
    sampled_steps=(8192,16384),event_rms_tolerance_A=.005,event_dc_tolerance_V=.001)

function prepare_engineering_domain(output)
    ispath(output) && error("domain input already exists; preserve previous evidence")
    refs=Any[];a=cis(2pi/3);positive=ComplexF64[1,a^2,a];negative=conj.(positive)
    for four in (false,true)
        fixture=Dict{String,Any}("ibr"=>Dict{String,Any}("pv"=>Dict{String,Any}("bus"=>"poc")))
        base=four ? first(values(four_leg_fleet(fixture,1.).devices)).device :
            plant("pv","poc";settings=TERMINAL_CONTROL_SETTINGS)
        inv=engineering_plant(base)
        geometry=PowerOptLab._switching_envelope_geometry(inv)
        bounds=PowerOptLab._inverter_current_envelope(inv)
        for p in ENGINEERING_DOMAIN.p_W,q in ENGINEERING_DOMAIN.q_var,
                i2 in ENGINEERING_DOMAIN.negative_sequence_A,phi in ENGINEERING_DOMAIN.sequence_angle_rad,
                level in ENGINEERING_DOMAIN.modulation_fraction,theta in ENGINEERING_DOMAIN.carrier_line_angle_rad
            u=positive .+ ENGINEERING_DOMAIN.voltage_negative_sequence_fraction .* cis(phi) .* negative
            needed=maximum(sqrt(2)*abs(u[j]-u[k]) for (j,k) in ((1,2),(2,3),(3,1)))
            four && (needed=max(needed,sqrt(2)*maximum(abs,u)))
            u.*=level*inv.m_max*inv.v_dc/needed
            u1=sum(u.*conj.(positive))/3;u2=sum(u.*conj.(negative))/3
            neg=i2*cis(phi);zero=four ? ENGINEERING_DOMAIN.four_leg_zero_sequence_A*cis(phi+pi/6) : 0im
            # Correct I1 so the complete sequence combination has exactly the
            # stated total converter P/Q, including its negative-sequence power.
            pos=conj((complex(p,q)/3-u2*conj(neg))/u1)
            current=pos.*positive .+ neg.*negative .+ zero
            v=sqrt(2).*real.(u.*cis(theta));fund=sqrt(2).*real.(current.*cis(theta))
            legs=four ? [v;0.] : v
            duties=.5 .+ (legs .- (maximum(legs)+minimum(legs))/2) ./ inv.v_dc
            leg_current=four ? [fund;-sum(fund)] : fund
            events=PowerOptLab._carrier_event_ripple(duties,inv.v_dc,leg_current,
                geometry,1/inv.f_sw,inv.c_dc)
            dc=PowerOptLab._dc_switching_voltage_bound(inv,real.(current),imag.(current))
            R,X=PowerOptLab._series_filter_matrices(inv,3,1.;side=:converter)
            primitive=complex.(R,X)*(four ? [current;-sum(current)] : current)
            drop=four ? primitive[1:3].-primitive[4] : primitive
            report=inverter_physical_feasibility(inv;converter_voltage=u,converter_current=current,
                grid_voltage=u-drop,filter_voltage=u-drop,pwm_thermal_headroom=PowerOptLab._PWM_ALLOCATION_HEADROOM)
            pairs(x)=[[real(z),imag(z)] for z in x]
            push!(refs,Dict("id"=>length(refs)+1,"topology"=>String(inv.topology),
                "p_W"=>p,"q_var"=>q,"i2_A"=>i2,"sequence_angle_rad"=>phi,
                "modulation_fraction"=>level,"carrier_line_angle_rad"=>theta,
                "converter_voltage"=>pairs(u),"converter_current"=>pairs(current),
                "duties"=>duties,"frozen_leg_current_A"=>leg_current,
                "event_phase_rms_A"=>events.phase_rms,"event_phase_peak_A"=>events.phase_peak,
                "event_neutral_rms_A"=>events.neutral_rms,"event_neutral_peak_A"=>events.neutral_peak,
                "event_dc_min_V"=>events.dc_min,"event_dc_max_V"=>events.dc_max,
                "phase_rms_bound_A"=>bounds.converter_rms,"phase_peak_bound_A"=>bounds.phase_peak,
                "neutral_rms_bound_A"=>bounds.neutral_rms,"neutral_peak_bound_A"=>bounds.neutral_peak,
                "dc_voltage_bound_V"=>dc.total,"dc_frozen_bound_V"=>dc.frozen,
                "dc_ac_variation_bound_V"=>dc.ac_variation,
                "engineering_limits_passed"=>report.passed,
                "failed_engineering_checks"=>[c.name for c in report.checks if !c.passed]))
        end
    end
    repo=normpath(joinpath(ROOT,"..",".."));source=study_source_sha256()
    for file in ("switching_engineering.jl","thermal_energy.jl","four_leg.jl","voltage_aware.jl",
            "prepare_engineering_domain.jl")
        path=joinpath(@__DIR__,file);source[relpath(path,repo)]=bytes2hex(sha256(read(path)))
    end
    mkpath(dirname(output))
    writejson(output,Dict("domain"=>ENGINEERING_DOMAIN,"cases"=>refs,"source_sha256"=>source,
        "hardware"=>HARDWARE,"operating_limits"=>ENGINEERING_OPERATING_LIMITS,"sizing_changed"=>false,
        "julia_version"=>string(VERSION),"git_revision"=>strip(read(Cmd(["git","-C",repo,"rev-parse","HEAD"]),String)),
        "git_status"=>read(Cmd(["git","-C",repo,"status","--porcelain"]),String),
        "Project_sha256"=>bytes2hex(sha256(read(joinpath(repo,"Project.toml")))),
        "Manifest_sha256"=>bytes2hex(sha256(read(joinpath(repo,"Manifest.toml")))),
        "scope"=>"frozen-bus, frozen-reference centered carrier with inductive AC ripple; periodic mean-subtracted gate-current charge; not a transient or thermal-calibration certificate"))
    println("Prepared ",length(refs)," carrier-domain references; engineering failures retained: ",
        count(r->!r["engineering_limits_passed"],refs))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==1 || error("usage: prepare_engineering_domain.jl output.json")
    prepare_engineering_domain(ARGS[1])
end
