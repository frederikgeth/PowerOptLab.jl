using Test, PowerOptLab, BMOPFTools, JSON3, JSONSchema

# Minimized high-ratio witness: raw secondary ohms alone cannot classify the
# finite-admittance domain (1e-7 * 100^2 = 1e-3 ohm on the HV coil basis).
function _adapter_transformer_net(; ratio=100.0, zto=1e-7, tap=1.0)
    Dict{String,Any}(
        "bus"=>Dict{String,Any}(
            b=>Dict{String,Any}("terminal_names"=>["a","n"],
                               "perfectly_grounded_terminals"=>["n"]) for b in ("hv","lv")),
        "voltage_source"=>Dict("s"=>Dict{String,Any}("bus"=>"hv", "terminal_map"=>["a"],
            "v_magnitude"=>[230.0*ratio], "v_angle"=>[0.0])),
        "transformer"=>Dict("single_phase"=>Dict("t"=>Dict{String,Any}(
            "bus_from"=>"hv", "bus_to"=>"lv", "terminal_map_from"=>["a","n"],
            "terminal_map_to"=>["a","n"], "v_nom_from"=>230.0*ratio, "v_nom_to"=>230.0,
            "r_series_from"=>0.0, "x_series_from"=>0.0,
            "r_series_to"=>zto, "x_series_to"=>0.0, "tap"=>tap))))
end

@testset "upstream finite-transformer domain in SE and operability" begin
    cases = [
        ("referred leakage", _adapter_transformer_net(), true),
        ("ideal", _adapter_transformer_net(zto=0.0), false),
        ("below cutoff", _adapter_transformer_net(ratio=1.0,zto=prevfloat(1e-6)), false),
        ("at cutoff", _adapter_transformer_net(ratio=1.0,zto=1e-6), false),
        ("above cutoff", _adapter_transformer_net(ratio=1.0,zto=nextfloat(1e-6)), true),
        ("tap-referral below cutoff", _adapter_transformer_net(ratio=1.0,zto=2e-6,tap=0.5), false),
        ("tap-referral above cutoff", _adapter_transformer_net(ratio=1.0,zto=5e-7,tap=2.0), true),
    ]
    for (label, net, supported) in cases
        @testset "$label" begin
            original = deepcopy(net)
            # Exercise actual JSON transport as well as the in-memory witness.
            for input in (net, JSON3.read(JSON3.write(net), Dict{String,Any}))
                scope = operability_scope_audit(input)
                @test scope["status"] == (supported ? :supported : :not_applicable)
                @test JSON3.read(JSON3.write(scope))["status"] ==
                    (supported ? "supported" : "not_applicable")
                report = state_estimator_preflight(input)
                @test report.supported == supported
                if supported
                    @test compile_state_estimator(input).preflight.supported
                else
                    @test any(f->f.code===:invalid_transformer &&
                        f.element=="transformer/single_phase/t", report.findings)
                    @test_throws SEUnsupportedNetwork compile_state_estimator(input)
                    ms = [Measurement(kind=:vmag,bus="lv",terminal="a",value=230.0,sigma=1.0)]
                    @test_throws SEUnsupportedNetwork solve_state_estimation(input, ms)
                    spec = OperabilitySpec(scaling_policy=BMOPFTools.SIUnitsScaling())
                    # No solved state is needed to return scope refusal.
                    @test check_opf_operability(input, Dict(); spec).status == :not_applicable
                    @test continue_opf_operability(input, Dict(); spec).status == :not_applicable
                    @test continue_opf_operability_pseudo_arclength(input, Dict(); spec).status == :not_applicable
                    @test locate_opf_operability_fold(input, Dict(); spec).status == :not_applicable
                end
            end
            @test net == original
        end
    end
    finite = _adapter_transformer_net()
    s = compile_state_estimator(finite)
    x = initial_state_estimator(s, SEParameters(s))
    @test x[s.free_state_map[("lv","a")]] ≈ 230.0
    @test helm_series(_adapter_transformer_net(zto=0.0)).converged

    # Partially singular n-winding leakage requires an augmented current even
    # though one of its leakage modes is finite.
    nw = _adapter_transformer_net()
    nw["bus"]["third"] = deepcopy(nw["bus"]["lv"])
    nw["transformer"] = Dict("n_winding"=>Dict("t"=>Dict{String,Any}(
        "windings"=>[Dict{String,Any}("bus"=>b,"terminal_map"=>["a","n"],
            "configuration"=>"WYE","v_nom"=>v,"r_winding"=>0.0)
            for (b,v) in (("hv",23000.0),("lv",230.0),("third",230.0))],
        "x_sc"=>Dict("1_2"=>0.0,"1_3"=>0.02,"2_3"=>0.02), "s_rating"=>100e3)))
    @test operability_scope_audit(nw)["status"] == :not_applicable
    @test !state_estimator_preflight(nw).supported
    @test_throws SEUnsupportedNetwork compile_state_estimator(nw)
    malformed = _adapter_transformer_net()
    malformed["transformer"]["single_phase"]["t"]["v_nom_from"] = 0.0
    @test operability_scope_audit(malformed)["status"] == :not_applicable
    @test !state_estimator_preflight(malformed).supported
end

@testset "Kron coordinate serialization policy uses raw schema validation" begin
    schema = JSONSchema.Schema(JSON3.read(read(joinpath(dirname(pathof(BMOPFTools)),
        "validation", "schemas", "draft_bmopf_schema.json"), String)))
    # Cover both writer paths: real neutral elimination and no-op early return.
    for neutral in (true, false), coordinates in (true, false)
        terminals = neutral ? ["a","n"] : ["a"]
        net = Dict{String,Any}(
            "bus"=>Dict("b"=>Dict{String,Any}("terminal_names"=>terminals)),
            "voltage_source"=>Dict("s"=>Dict{String,Any}("bus"=>"b",
                "terminal_map"=>["a"],"v_magnitude"=>[230.0],"v_angle"=>[0.0])))
        neutral && (net["bus"]["b"]["perfectly_grounded_terminals"] = ["n"])
        if coordinates
            net["bus"]["b"]["longitude"] = 151.0
            net["bus"]["b"]["latitude"] = -27.0
        end
        original = deepcopy(net)
        for input in (net, JSON3.write(net))
            dict = kron_reduce_bmopf(input)
            @test haskey(dict["bus"]["b"],"longitude") == coordinates
            strict = JSON3.read(kron_reduce_bmopf(input; as_json=true), Dict{String,Any})
            @test JSONSchema.validate(schema, strict) === nothing
            @test !haskey(strict["bus"]["b"],"longitude")
            @test !haskey(strict["bus"]["b"],"latitude")
            extended = JSON3.read(kron_reduce_bmopf(input; as_json=true, include_coordinates=true), Dict{String,Any})
            @test (JSONSchema.validate(schema, extended) === nothing) == !coordinates
            if coordinates
                @test extended["bus"]["b"]["longitude"] == 151.0
                @test extended["bus"]["b"]["latitude"] == -27.0
                delete!(extended["bus"]["b"],"longitude")
                delete!(extended["bus"]["b"],"latitude")
            end
            @test extended == strict # identical electrical data and provenance
        end
        @test net == original
    end
end
