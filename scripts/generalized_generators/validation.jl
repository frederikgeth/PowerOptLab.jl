# Independent complex-circuit and equality-rank checks. No timing claim is made
# against another OPF implementation; reported times describe this process only.
# Run: julia --project=. scripts/generalized_generators/validation.jl
include(joinpath(@__DIR__,"tradeoffs.jl"))

function generator_equality_audit(model)
    vars=JuMP.all_variables(model); index=Dict(v=>k for (k,v) in enumerate(vars))
    rows=Vector{Float64}[]
    for ref in JuMP.all_constraints(model;include_variable_in_set_constraints=true)
        c=JuMP.constraint_object(ref)
        c.set isa JuMP.MOI.EqualTo || continue
        f=c.func; row=zeros(length(vars))
        if f isa JuMP.VariableRef
            row[index[f]]=1
        else
            for (coef,v) in JuMP.linear_terms(f); row[index[v]]+=coef; end
            if f isa JuMP.QuadExpr
                for (coef,v,w) in JuMP.quad_terms(f)
                    row[index[v]]+=coef*JuMP.value(w)
                    row[index[w]]+=coef*JuMP.value(v)
                end
            end
        end
        nr=norm(row);nr>0 && (row./=nr)
        push!(rows,row)
    end
    J=isempty(rows) ? zeros(0,length(vars)) : reduce(vcat,permutedims.(rows))
    s=svdvals(J); tol=1e-9
    (variables=length(vars),equalities=size(J,1),rank=count(>(tol),s),
     smallest=isempty(s) ? NaN : minimum(s),jacobian=J)
end

function generator_validation_study()
    net=generator_tradeoff_network();empty!(net["load"])
    grid=vcat([240,225,232].*cis.([0.03,-2.12,2.10]),0.0)
    L=Diagonal(ComplexF64[0.1+0.3im,0.1+0.3im,0.1+0.3im,0.2+0.1im])
    ephase=[235,227,233].*cis.([0.02,-2.11,2.08])
    Cw=vcat(Matrix{Float64}(I,3,3),-ones(1,3))
    Cd=[1.0 0 -1;-1 1 0;0 -1 1;0 0 0]
    options=Pair{String,Any}["tol"=>1e-9,"constr_viol_tol"=>1e-9,"bound_relax_factor"=>0.0,"max_iter"=>1200]
    cases=NamedTuple[]
    for (connection,C) in ((:wye,Cw),(:delta,Cd))
        e=connection==:wye ? ephase : Cd[1:3,:]'*ephase
        matrices=["zero"=>zeros(ComplexF64,3,3),
            "diagonal"=>Matrix(Diagonal(ComplexF64[0.4+0.6im,0.5+0.7im,0.3+0.8im])),
            "mutual"=>ComplexF64[0.4+0.6im 0.02+0.03im 0.01;0.02+0.03im 0.5+0.7im 0.04im;0.01 0.04im 0.3+0.8im],
            "singular"=>(0.4+0.8im)*(Cd'*Cd)]
        for (name,Z) in matrices
            A=Z+C'*L*C
            ideal_loop=connection==:delta && all(iszero,sum(Z;dims=1))
            # The singular oracle explicitly selects zero circulation. The
            # optimization gets the same declared zero-sequence current limit.
            i=ideal_loop ? pinv(A)*(e-C'*grid) : A\(e-C'*grid)
            cap=ideal_loop ? GeneratorCapability(sequence=GeneratorSequenceLimits(current_max=[0.0,1000.0,1000.0])) : GeneratorCapability()
            ports=connection==:wye ? [("a","n"),("b","n"),("c","n")] : [("a","b"),("b","c"),("c","a")]
            d=GeneralizedGenerator(id="g",bus="poc",connections=ports,impedance=Z,
                voltage=GeneratorVoltageLaw(:fixed_phasor;phasor=e),capability=cap)
            j=C*i;v=grid+L*j
            push!(cases,(label="$(connection)/$name",device=d,i=i,j=connection==:wye ? j : j[1:3],
                v=connection==:wye ? v : v[1:3]))
        end
    end
    Z=Matrix(Diagonal(ComplexF64[0.4+0.6im,0.5+0.7im,0.3+0.8im,0.2+0.1im]))
    Z[1,4]=Z[4,1]=0.02+0.01im
    for (label,zg) in (("open",nothing),("finite",0.8+0.1im),("ideal",0.0))
        M=L+Z
        phase=M[1:3,:]-ones(3)*transpose(M[4,:])
        ground=zg===nothing ? ones(1,4) : transpose(M[4,:])+zg*ones(1,4)
        A=vcat(phase,ground)
        j=A\vcat(ephase-grid[1:3],0.0)
        d=SourceGenerator(id="g",bus="poc",impedance=Z,grounding=zg,
            voltage=GeneratorVoltageLaw(:fixed_phasor;phasor=ephase))
        push!(cases,(label="source/$label",device=d,i=j[1:3],j=j,v=grid+L*j))
    end
    for connection in (:wye,:delta,:source),mode in (:none,:fixed_phasor,:fixed_magnitudes,:common_magnitude,:phase_magnitudes)
        delta=connection==:delta
        e=delta ? Cd[1:3,:]'*ephase : ephase
        law=mode in (:fixed_phasor,:fixed_magnitudes) ? GeneratorVoltageLaw(mode;phasor=e) :
            mode==:none ? GeneratorVoltageLaw() : GeneratorVoltageLaw(mode;
                magnitude_min=delta ? 320.0 : 180.0,magnitude_max=delta ? 470.0 : 280.0)
        ctl=mode==:fixed_phasor ? GeneratorControl() : mode==:fixed_magnitudes ? GeneratorControl(p=1800.0) :
            mode==:none || mode==:phase_magnitudes ? GeneratorControl(p=1800.0,q=[50.0,100.0,150.0]) :
            GeneratorControl(p=1800.0,q=300.0)
        args=(id="g",bus="poc",voltage=law,control=ctl)
        d=connection==:source ? SourceGenerator(;args...,impedance=fill(0.8+1.2im,4),grounding=0.8+0.1im) :
            GeneralizedGenerator(;args...,impedance=0.8+1.2im,
                connections=delta ? [("a","b"),("b","c"),("c","a")] : [("a","n"),("b","n"),("c","n")])
        push!(cases,(label="$(connection)/$(mode)",device=d,i=nothing,j=nothing,v=nothing))
    end
    rows=NamedTuple[]
    # Warm one ordinary case; other specializations may still compile.
    warm=build_generator_model(net,[cases[1].device];per_unit=true,s_base=1e4,optimizer=JuMP.optimizer_with_attributes(Ipopt.Optimizer,options...))
    BMOPFTools.enforce_kcl!(warm.context);JuMP.optimize!(BMOPFTools.opf_model(warm.context))
    for case in cases,(pu,base) in ((false,1e4),(true,1e4),(true,1e8))
        started=time_ns()
        built=build_generator_model(net,[case.device];per_unit=pu,s_base=base,optimizer=JuMP.optimizer_with_attributes(Ipopt.Optimizer,options...))
        BMOPFTools.enforce_kcl!(built.context)
        build_ms=(time_ns()-started)/1e6
        model=BMOPFTools.opf_model(built.context)
        started=time_ns();JuMP.optimize!(model);solve_ms=(time_ns()-started)/1e6
        status=SolveStatus(PowerOptLab._solve_outcome(model))
        @assert status.publishable "$(case.label): $status"
        g=extract_device(case.device,built.handles["g"],status)
        verr=case.v===nothing ? NaN : maximum(abs.(g.terminal_voltage-case.v))
        ierr=case.i===nothing ? NaN : max(maximum(abs.(g.port_current-case.i)),maximum(abs.(g.terminal_current-case.j)))
        case.v===nothing || @assert verr<1e-5 && ierr<1e-5
        @assert abs(g.power_balance_error)<1e-4
        audit=generator_equality_audit(model)
        @assert audit.rank==audit.equalities "dependent equality rows in $(case.label)"
        @assert built.handles["g"].internal_voltage_variables==0
        push!(rows,(case=case.label,per_unit=pu,base=base,voltage_error=verr,current_error=ierr,
            rank=audit.rank,equalities=audit.equalities,variables=audit.variables,
            smallest=audit.smallest,build_ms=build_ms,solve_ms=solve_ms))
    end
    rows
end

function print_generator_validation(rows)
    @printf("%-26s %3s %8s %10s %10s %8s %7s %10s %9s %9s\n",
        "Circuit","pu","Sbase","V error","I error","rank/eq","vars","min sv","build ms","solve ms")
    for r in rows
        @printf("%-26s %3s %8.0e %10.2e %10.2e %3d/%-4d %7d %10.2e %9.2f %9.2f\n",
            r.case,r.per_unit,r.base,r.voltage_error,r.current_error,r.rank,r.equalities,r.variables,r.smallest,r.build_ms,r.solve_ms)
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    print_generator_validation(generator_validation_study())
end
