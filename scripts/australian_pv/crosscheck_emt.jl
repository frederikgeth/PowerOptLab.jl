# Map measured non-2ω RMS stress to the plant's existing fixed i_sw diagnostic.
# This is a point-specific rejection check, never a fleet-wide reserve estimate.
include("prepare_emt.jl")
using Test

function crosscheck_emt(input,folder)
    env=JSON3.read(read(input,String),Dict{String,Any})
    result=JSON3.read(read(joinpath(folder,"results.json"),String),Dict{String,Any})
    result["input_sha256"]==bytes2hex(sha256(read(input))) || error("EMT input mismatch")
    result["circuit_sha256"]==bytes2hex(sha256(read(joinpath(@__DIR__,"emt_circuit.jl")))) || error("stale circuit")
    result["all_mechanisms_accepted"] || error("mechanism audit required before stress cross-check")
    phasors(v)=ComplexF64[complex(z[1],z[2]) for z in v]
    rows=Dict{String,Any}()
    @testset "EMT stress maps to the existing capacitor engineering gate" begin
        for c in env["cases"]
            m=result["cases"][c["case_id"]]["modes"]["switched"]["finest"]
            lf_rms=abs(complex(m["capacitor_2omega_peak"]...))/sqrt(2)
            reserve=sqrt(max(0.,m["capacitor_total_rms_A"]^2-lf_rms^2))
            hw=HARDWARE
            function plant_with_reserve(value)
                AdvancedInverter(id="emt",bus="poc",topology=Symbol(c["topology"]),
                    phase_terminals=["a","b","c"],neutral=c["topology"]=="FOUR_LEG" ? "n" : nothing,
                    s_max=hw.s_max,power_limit_location=:poc,i_max=hw.i_max,In_max=c["neutral_current_rating_A"],
                    v_dc=hw.v_dc,c_dc=hw.c_dc,m_max=hw.m_max,dv2_max=hw.dv2_max,i_cap_max=hw.i_cap_max,
                    i_sw=value,r_filter=c["phase_R_ohm"],x_filter=c["phase_L_H"]*2pi*c["frequency_Hz"],
                    r_filter_neutral=c["neutral_R_ohm"],x_filter_neutral=c["neutral_L_H"]*2pi*c["frequency_Hz"],
                    p_loss_fixed=hw.p_loss_fixed,a_loss=hw.a_loss,c_loss=hw.c_loss)
            end
            function report(inv)
                inverter_physical_feasibility(inv;converter_voltage=phasors(c["converter_voltage"]),
                    converter_current=phasors(c["current"]),grid_voltage=phasors(c["grid_voltage"]),
                    grid_current=phasors(c["current"]),dc_min=0.,dc_max=5000.)
            end
            base=report(plant_with_reserve(0.))
            base_rms=only(ch.value for ch in base.checks if ch.name=="capacitor_current")
            # A reserve already larger than the entire budget is rejected by
            # configuration validation before any operating point is stamped.
            rejection=try
                validate_device(plant_with_reserve(reserve));nothing
            catch err
                err
            end
            thermal=hypot(base_rms,reserve)
            @test base.passed
            @test rejection isa ArgumentError
            @test occursin("weighted i_sw must be smaller",sprint(showerror,rejection))
            @test thermal ≈ m["capacitor_total_rms_A"] atol=.001
            @test reserve>hw.i_cap_max
            rows[c["case_id"]]=Dict("fundamental_reference_accepted"=>base.passed,
                "measured_total_capacitor_rms_A"=>m["capacitor_total_rms_A"],"point_specific_non_2omega_reserve_A"=>reserve,
                "phasor_gate_with_reserve_rms_A"=>thermal,"phasor_gate_with_reserve_accepted"=>false,
                "configuration_rejection"=>sprint(showerror,rejection))
        end
    end
    writejson(joinpath(folder,"stress_crosscheck.json"),Dict("cases"=>rows,
        "scope"=>"point-specific measured non-2omega stress includes other harmonics and any DC drift; not a conservative reserve over dispatch or a PWM-aware allocator",
        "input_sha256"=>result["input_sha256"],"results_sha256"=>bytes2hex(sha256(read(joinpath(folder,"results.json")))),
        "crosscheck_code_sha256"=>bytes2hex(sha256(read(@__FILE__)))))
    rows
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==2 || error("usage: crosscheck_emt.jl input.json results_folder")
    crosscheck_emt(ARGS[1],ARGS[2])
end
