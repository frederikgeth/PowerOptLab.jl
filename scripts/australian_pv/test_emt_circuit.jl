# Independent circuit oracles. This file does not load PowerOptLab or its lifts.
include("emt_circuit.jl")
using Test

function emt_fixture(;four=false)
    a=cis(2pi/3);i=4 .* ComplexF64[1,a^2,a];v=230 .* ComplexF64[1,a^2,a]
    pairs(v)=[[real(z),imag(z)] for z in v]
    Dict{String,Any}("case_id"=>"analytic","topology"=>four ? "FOUR_LEG" : "THREE_LEG",
        "grid_voltage"=>pairs(v),"converter_voltage"=>pairs(v.+(.2+im*2pi*50*.004).*i),"current"=>pairs(i),
        "phase_R_ohm"=>.2,"phase_L_H"=>.004,"neutral_R_ohm"=>.1,"neutral_L_H"=>.002,
        "frequency_Hz"=>50.,"dc_voltage_V"=>750.,"dc_capacitance_F"=>300e-6,
        "p_dc_W"=>3*230*4+3*.2*4^2,"mean_converter_loss_W"=>0.,"p_available_W"=>5000.,
        "phase_current_rating_A"=>10.,"neutral_current_rating_A"=>10.,
        "capacitor_current_rating_A"=>1.,"dc_2omega_peak_rating_V"=>7.5,"modulation_utilization"=>.96)
end

@testset "Independent EMT circuit: energy and neutral inductance" begin
    for four in (false,true)
        p=IBREMT.parameters(emt_fixture(;four));x=(4.,-1.,four ? -2. : -3.,750.);g=(1.,0.,1.,0.)
        dx,power=IBREMT.rhs(p,.001,x,g,:switched)
        grad=(ntuple(k->p.l*x[k]+(four ? p.ln*sum(x[1:3]) : 0.),3)...,p.capacitance*x[4])
        @test sum(grad[k]*dx[k] for k in 1:4) ≈ power atol=1e-9
        four || @test abs(sum(dx[1:3]))<1e-9
    end
    c=emt_fixture(four=true)
    c["grid_voltage"]=fill([0.,0.],3);c["p_dc_W"]=0.
    p=IBREMT.parameters(c)
    dx,_=IBREMT.rhs(p,0.,(0.,0.,0.,750.),(0.,0.,0.,1.),:switched)
    @test collect(dx[1:3]) ≈ fill(-750/(.004+3*.002),3) atol=1e-9
    @test dx[4]==0.
    p=IBREMT.parameters(emt_fixture())
    d=7.5+0im;initial=IBREMT.initial_dc_voltage(p,d)
    center=initial^2-2p.vdc*real(d)
    average=sum(sqrt(center+2p.vdc*abs(d)*cos(2pi*k/4096)) for k in 0:4095)/4096
    @test average ≈ p.vdc atol=1e-9
    @test_throws ErrorException IBREMT.initial_dc_voltage(p,p.vdc+0im)
end

@testset "Independent EMT circuit: analytic balanced averaged orbit" begin
    for four in (false,true)
        r=IBREMT.simulate(emt_fixture(;four);cycles=2,discard_cycles=1,steps_per_carrier=32,mode=:averaged)
        @test r.metrics["fundamental_current_max_error_A"]<1e-6
        @test r.metrics["dc_2omega_error_V"]<1e-6
        @test r.metrics["capacitor_total_rms_A"]<1e-6
        @test r.metrics["energy_relative_residual"]<1e-10
        @test r.metrics["v_dc_mean_V"] ≈ 750. atol=1e-6
        @test r.metrics["p_grid_mean_W"] ≈ 2760. atol=1e-5
        @test r.metrics["all_engineering_checks_passed"]
    end
    c=emt_fixture();c["dc_capacitance_F"]=-1.
    @test_throws ErrorException IBREMT.parameters(c)
    @test_throws ArgumentError IBREMT.simulate(emt_fixture();mode=:unknown)
    @test_throws ArgumentError IBREMT.simulate(emt_fixture();cycles=2,discard_cycles=2)
end
