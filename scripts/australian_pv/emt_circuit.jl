"""Independent time-stepping L-filter/bridge/DC-capacitor circuit; no OPF imports."""
module IBREMT
using JSON3, SHA

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const SQRT2=sqrt(2.)
phasors(v)=ntuple(k->complex(v[k][1],v[k][2]),3)
wave(v,t,w)=ntuple(k->SQRT2*real(v[k]*cis(w*t)),3)
plus(x,y,a)=ntuple(k->x[k]+a*y[k],4)

function parameters(c)
    p=(four=c["topology"]=="FOUR_LEG",vg=phasors(c["grid_voltage"]),
        uc=phasors(c["converter_voltage"]),iref=phasors(c["current"]),
        r=Float64(c["phase_R_ohm"]),l=Float64(c["phase_L_H"]),
        rn=Float64(c["neutral_R_ohm"]),ln=Float64(c["neutral_L_H"]),
        w=2pi*c["frequency_Hz"],f=Float64(c["frequency_Hz"]),
        capacitance=Float64(c["dc_capacitance_F"]),vdc=Float64(c["dc_voltage_V"]),
        psource=Float64(c["p_dc_W"]),loss=Float64(c["mean_converter_loss_W"]),
        m=Float64(c["modulation_utilization"]))
    all(isfinite,(p.r,p.l,p.rn,p.ln,p.w,p.capacitance,p.vdc,p.psource,p.loss,p.m)) || error("nonfinite EMT parameters")
    p.l>0 && p.capacitance>0 && p.vdc>0 && p.r>=0 && p.rn>=0 && p.ln>=0 &&
        p.psource>=0 && p.loss>=0 && 0<p.m<=1 || error("invalid EMT parameters")
    p
end

# The physical neutral is eliminated by its KCL, not set to zero. Four-wire
# magnetic energy uses L*I + Ln*11'; three-wire currents stay on sum(i)=0.
function rhs(p,t,x,g,mode)
    voltage=wave(p.vg,t,p.w);total=x[1]+x[2]+x[3]
    bridge=mode===:averaged ? wave(p.uc,t,p.w) :
        ntuple(k->x[4]*(p.four ? g[k]-g[4] : g[k]-(g[1]+g[2]+g[3])/3),3)
    offset=p.four ? p.rn*total : -(voltage[1]+voltage[2]+voltage[3])/3
    u=ntuple(k->bridge[k]-voltage[k]-p.r*x[k]-offset,3)
    common=p.four ? p.ln/(p.l+3p.ln)*sum(u) : 0.
    di=ntuple(k->(u[k]-common)/p.l,3)
    pac=sum(bridge[k]*x[k] for k in 1:3)
    dv=(p.psource-p.loss-pac)/(p.capacitance*x[4])
    filterloss=p.r*sum(x[k]^2 for k in 1:3)+(p.four ? p.rn*total^2 : 0.)
    net=p.psource-p.loss-sum(voltage[k]*x[k] for k in 1:3)-filterloss
    ((di...,dv),net)
end

energy(p,x)=p.capacitance*x[4]^2/2+p.l*sum(x[k]^2 for k in 1:3)/2+
    (p.four ? p.ln*(x[1]+x[2]+x[3])^2/2 : 0.)

# Match the declared MEAN DC voltage in the ideal averaged orbit. Simply
# initializing Vdc + Re(D) instead fixes an energy-center voltage and shifts
# the actual mean by O(|D|²/Vdc), which matters at a binding ripple boundary.
function initial_dc_voltage(p,d)
    amplitude=2p.vdc*abs(d)
    amplitude<p.vdc^2 || error("DC energy ripple reaches zero voltage")
    lower=p.vdc^2;upper=lower+amplitude
    for _ in 1:55
        center=(lower+upper)/2
        mean=sum(sqrt(center+amplitude*cos(2pi*k/256)) for k in 0:255)/256
        mean<p.vdc ? (lower=center) : (upper=center)
    end
    sqrt((lower+upper)/2+2p.vdc*real(d))
end

function step(p,t,x,g,mode,h)
    a,pa=rhs(p,t,x,g,mode)
    xa=plus(x,a,h/2);b,pb=rhs(p,t+h/2,xa,g,mode)
    xb=plus(x,b,h/2);c,pc=rhs(p,t+h/2,xb,g,mode)
    xc=plus(x,c,h);d,pd=rhs(p,t+h,xc,g,mode)
    next=ntuple(k->x[k]+h/6*(a[k]+2b[k]+2c[k]+d[k]),4)
    quadrature=((t,x,h/6),(t+h/2,xa,h/3),(t+h/2,xb,h/3),(t+h,xc,h/6))
    next,quadrature,h/6*(pa+2pb+2pc+pd)
end

function duties(p,t,x,period)
    # Regular sampled feedforward; the carrier midpoint is the voltage sample.
    # Predict only the slow DC slope over half a carrier for compensation.
    tm=t+period/2;u=wave(p.uc,tm,p.w);i=wave(p.iref,tm,p.w)
    slope=(p.psource-p.loss-sum(u[k]*i[k] for k in 1:3))/(p.capacitance*x[4])
    v=x[4]+period/2*slope
    r=p.four ? (u[1]/v,u[2]/v,u[3]/v,0.) : (u[1]/v,u[2]/v,u[3]/v)
    common=(maximum(r)+minimum(r))/2
    d=map(z->.5+z-common,r)
    margin=p.m-(maximum(r)-minimum(r))
    all(z->0<=z<=1,d) || error("PWM duty saturation: no silent clipping")
    d,margin
end

function simulate(c;switching_Hz=10000.,cycles=10,discard_cycles=6,steps_per_carrier=64,mode=:switched)
    mode in (:averaged,:switched) || throw(ArgumentError("unsupported EMT mode"))
    0<=discard_cycles<cycles && steps_per_carrier>=4 || throw(ArgumentError("invalid time window/mesh"))
    p=parameters(c);ratio=switching_Hz/p.f
    isinteger(ratio) || error("integer carrier/line ratio required for declared Fourier window")
    period=1/switching_Hz;n=round(Int,cycles*ratio);discard=round(Int,discard_cycles*ratio)
    window=(cycles-discard_cycles)/p.f
    s2=sum(p.uc[k]*p.iref[k] for k in 1:3);expected_d=im*s2/(2p.w*p.capacitance*p.vdc)
    x=(ntuple(k->SQRT2*real(p.iref[k]),3)...,initial_dc_voltage(p,expected_d))
    initial=NaN;initial_state=x;energy_input=0.;fund=zeros(ComplexF64,3);isq=zeros(3);ipeak=zeros(3)
    vmean=0.;v2=0im;icmean=0.;ic2=0im;icsq=0.;icpeak=0.;nsq=0.;pgrid=0.;minmargin=Inf
    trace=NamedTuple[];carrier=NamedTuple[];voltage_min=Inf;voltage_max=-Inf
    for index in 0:n-1
        start=index*period;active=index>=discard
        index==discard && (initial=energy(p,x);initial_state=x)
        d,margin=duties(p,start,x,period);minmargin=min(minmargin,margin)
        # Exact gate event locations for centered PWM. Every gate discontinuity
        # is an integration boundary; substeps cannot straddle a switch edge.
        edges=mode===:switched ? sort!([(period*(1-d[k])/2,k,1.) for k in eachindex(d)] ∪
            [(period*(1+d[k])/2,k,0.) for k in eachindex(d)];by=first) : Tuple{Float64,Int,Float64}[]
        push!(edges,(period,0,0.));g=(0.,0.,0.,0.);left=0.
        active && push!(carrier,(time_s=start,v_dc_V=x[4],i_a_A=x[1],i_b_A=x[2],i_c_A=x[3]))
        for (right,leg,value) in edges
            length=right-left
            count=max(1,ceil(Int,length/(period/steps_per_carrier)))
            h=length/count
            for sub in 1:count
                h>0 || continue
                t=start+left+(sub-1)*h
                next,quadrature,balance=step(p,t,x,g,mode,h)
                all(isfinite,next) && next[4]>0 || error("nonphysical EMT state")
                if active
                  for (tm,mid,weight) in quadrature
                    ii=(mid[1],mid[2],mid[3]);vg=wave(p.vg,tm,p.w)
                    bridge=mode===:averaged ? wave(p.uc,tm,p.w) : ntuple(k->mid[4]*(p.four ? g[k]-g[4] : g[k]-(g[1]+g[2]+g[3])/3),3)
                    ic=(p.psource-p.loss-sum(bridge[k]*ii[k] for k in 1:3))/mid[4]
                    for k in 1:3
                        fund[k]+=weight*SQRT2*ii[k]*cis(-p.w*tm)
                        isq[k]+=weight*ii[k]^2;ipeak[k]=max(ipeak[k],abs(ii[k]))
                    end
                    vmean+=weight*mid[4];v2+=2weight*mid[4]*cis(-2p.w*tm)
                    icmean+=weight*ic;ic2+=2weight*ic*cis(-2p.w*tm);icsq+=weight*ic^2
                    icpeak=max(icpeak,abs(ic))
                    nsq+=weight*sum(ii)^2;pgrid+=weight*sum(vg[k]*ii[k] for k in 1:3)
                    voltage_min=min(voltage_min,mid[4]);voltage_max=max(voltage_max,mid[4])
                    index>=n-2 && push!(trace,(time_s=tm,v_dc_V=mid[4],capacitor_current_A=ic,
                        i_a_A=ii[1],i_b_A=ii[2],i_c_A=ii[3],neutral_return_A=-sum(ii)))
                  end
                  energy_input+=balance
                end
                x=next
            end
            leg==0 || (g=ntuple(k->k==leg ? value : g[k],4))
            left=right
        end
    end
    fund./=window;actual_d=v2/window;actual_ic2=ic2/window
    expected_ic2=im*2p.w*p.capacitance*expected_d
    engineering=Dict("phase_rms"=>all(sqrt.(isq./window).<=c["phase_current_rating_A"]+1e-5),
        "neutral_rms"=>!p.four || sqrt(nsq/window)<=c["neutral_current_rating_A"]+1e-5,
        "capacitor_total_rms"=>sqrt(icsq/window)<=c["capacitor_current_rating_A"]+1e-5,
        "dc_2omega_peak"=>abs(actual_d)<=c["dc_2omega_peak_rating_V"]+1e-5,
        "mean_dc_source"=>p.psource<=c["p_available_W"]+.01,"modulation"=>minmargin>=-1e-8)
    metric=Dict("case_id"=>c["case_id"],"mode"=>String(mode),"steps_per_carrier"=>steps_per_carrier,
        "current_fundamental"=>[[real(i),imag(i)] for i in fund],"current_total_rms_A"=>sqrt.(isq./window),
        "current_peak_A"=>ipeak,"neutral_total_rms_A"=>sqrt(nsq/window),
        "capacitor_total_rms_A"=>sqrt(icsq/window),"capacitor_peak_A"=>icpeak,"capacitor_mean_A"=>icmean/window,
        "capacitor_2omega_peak"=>[real(actual_ic2),imag(actual_ic2)],
        "capacitor_2omega_error_A"=>abs(actual_ic2-expected_ic2),
        "v_dc_mean_V"=>vmean/window,"v_dc_min_V"=>voltage_min,"v_dc_max_V"=>voltage_max,
        "dc_2omega_peak"=>[real(actual_d),imag(actual_d)],"expected_dc_2omega_peak"=>[real(expected_d),imag(expected_d)],
        "fundamental_current_max_error_A"=>maximum(abs.(fund.-collect(p.iref))),
        "dc_2omega_error_V"=>abs(actual_d-expected_d),"p_grid_mean_W"=>pgrid/window,
        "window_initial_state"=>collect(initial_state),"window_final_state"=>collect(x),
        "window_s"=>window,"window_net_energy_J"=>energy_input,
        "energy_residual_J"=>energy(p,x)-initial-energy_input,
        "energy_relative_residual"=>abs(energy(p,x)-initial-energy_input)/max(1.,abs(p.psource)*window),
        "minimum_modulation_margin"=>minmargin,"engineering_checks"=>engineering,
        "all_engineering_checks_passed"=>all(values(engineering)))
    (metrics=metric,carrier=carrier,trace=trace)
end

function writecsv(path,rows)
    isempty(rows) && return
    open(path,"w") do io
        println(io,join(string.(keys(first(rows))),","))
        for row in rows;println(io,join(values(row),","));end
    end
end

function run(input,output)
    isdir(output) && !isempty(readdir(output)) && error("EMT output must be empty")
    env=JSON3.read(read(input,String),Dict{String,Any})
    for (path,digest) in env["source_sha256"]
        bytes2hex(sha256(read(joinpath(ROOT,path))))==digest || error("stale EMT input: $path")
    end
    mkpath(output);timing=env["timing"];bands=env["validation_bands"]
    results=Dict{String,Any}();all_metrics=Dict{String,Any}[]
    for c in env["cases"]
        solves=Any[]
        for mode in (:averaged,:switched),mesh in timing["steps_per_carrier"]
            result=simulate(c;switching_Hz=timing["switching_Hz"],cycles=timing["cycles"],
                discard_cycles=timing["discard_cycles"],steps_per_carrier=mesh,mode)
            push!(solves,result);push!(all_metrics,result.metrics)
            if mesh==last(timing["steps_per_carrier"])
                prefix=c["case_id"]*"_"*String(mode)
                writecsv(joinpath(output,prefix*"_carrier.csv"),result.carrier)
                writecsv(joinpath(output,prefix*"_switching.csv"),result.trace)
            end
        end
        accepted=Bool[]
        details=Dict{String,Any}()
        for mode in ("averaged","switched")
            rows=[r.metrics for r in solves if r.metrics["mode"]==mode]
            fine=last(rows);coarse=rows[end-1]
            iref=phasors(c["current"]);d=complex(fine["expected_dc_2omega_peak"]...)
            physical= fine["fundamental_current_max_error_A"]<=bands["fundamental_absolute_A"]+bands["fundamental_relative"]*maximum(abs,iref) &&
                fine["dc_2omega_error_V"]<=bands["dc_2omega_absolute_V"]+bands["dc_2omega_relative"]*abs(d) &&
                fine["capacitor_2omega_error_A"]<=2pi*2c["frequency_Hz"]*c["dc_capacitance_F"]*(bands["dc_2omega_absolute_V"]+bands["dc_2omega_relative"]*abs(d)) &&
                fine["energy_relative_residual"]<=bands["energy_relative"]
            mesh_i=maximum(abs.(collect(phasors(fine["current_fundamental"])).-collect(phasors(coarse["current_fundamental"]))))
            mesh_d=abs(complex(fine["dc_2omega_peak"]...)-complex(coarse["dc_2omega_peak"]...))
            mesh_cap=abs(fine["capacitor_total_rms_A"]-coarse["capacitor_total_rms_A"])
            converged=mesh_i<=bands["mesh_current_A"] && mesh_d<=bands["mesh_dc_2omega_V"] && mesh_cap<=bands["mesh_current_A"]
            details[mode]=Dict("mechanisms_accepted"=>physical && converged,"mesh_current_difference_A"=>mesh_i,
                "mesh_dc_difference_V"=>mesh_d,"mesh_capacitor_rms_difference_A"=>mesh_cap,"finest"=>fine)
            push!(accepted,physical && converged)
        end
        results[c["case_id"]]=Dict("mechanisms_accepted"=>all(accepted),"modes"=>details)
        open(io->JSON3.pretty(io,Dict("cases"=>results,"completed"=>length(results),
            "all_mechanisms_accepted"=>length(results)==length(env["cases"]) && all(r["mechanisms_accepted"] for r in values(results)),
            "all_switched_hardware_accepted"=>length(results)==length(env["cases"]) && all(r["modes"]["switched"]["finest"]["all_engineering_checks_passed"] for r in values(results)),
            "planned"=>length(env["cases"]),"input_sha256"=>bytes2hex(sha256(read(input))),
            "waveform_sha256"=>Dict(file=>bytes2hex(sha256(read(joinpath(output,file)))) for file in readdir(output) if endswith(file,".csv")),
            "circuit_sha256"=>bytes2hex(sha256(read(@__FILE__))),"metrics"=>all_metrics)),joinpath(output,"results.json"),"w")
        println(c["case_id"]," mechanisms accepted=",all(accepted)," switched capacitor RMS=",last(solves).metrics["capacitor_total_rms_A"]);flush(stdout)
    end
    results
end

end # module

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==2 || error("usage: emt_circuit.jl input.json empty_output_directory")
    IBREMT.run(ARGS[1],ARGS[2])
end
