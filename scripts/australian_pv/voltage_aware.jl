# Study network and solver settings shared by the matched PF/OPF benchmark.
include("run.jl")

const VA_OPTIONS=("tol"=>1e-9,"max_iter"=>1000,"bound_relax_factor"=>0.,
    "mu_strategy"=>"adaptive","nlp_scaling_method"=>"none")
const VA_BASE=1e5

function voltage_aware_network(base,case)
    net=study_network(base,case.pv_fraction,case.source_pu;
        source_vuf=case.source_vuf,phase_shares=case.phase_shares)
    # Rotate U2 relative to U1 without altering |U1| or adding U0.
    a=cis(2pi/3)
    for source in values(net["voltage_source"])
        v=source["v_magnitude"].*cis.(source["v_angle"])
        u1=(v[1]+a*v[2]+a^2*v[3])/3
        u2=case.source_vuf*u1*cis(case.source_negative_angle_rad)
        v=u1.*[1,a^2,a].+u2.*[1,a,a^2]
        source["v_magnitude"]=abs.(v);source["v_angle"]=angle.(v)
    end
    for load in values(net["load"])
        load["p_nom"].*=case.load_scale;load["q_nom"].*=case.load_scale
    end
    net
end

# A different injection formulation in the same network engine. This is a
# fixed-PQ cross-check, not an independent circuit or external simulator.
