# Coupling of AC switching ripple into the DC-link current. The optimizer uses
# an all-duty bound; the event integral remains independent of the NLP.

function _carrier_ac_bridge_variance(d,K)
    # On the first half of a centered carrier, h_j=-d_j*t+(t-a_j)_+.
    # i=Vdc*T*K*h and b=g'i. On the second half b is odd, hence has zero
    # mean and is orthogonal to the even frozen-fundamental bridge current.
    a=(1 .-d)./2
    events=sort([zero(eltype(d));a;one(eltype(d))/2])
    result=zero(eltype(d))
    for j in 1:length(events)-1
        left,right=events[j:j+1];width=right-left
        gate=(left+right)/2 .>= a
        h(t)=[-d[k]*t+(gate[k] ? t-a[k] : zero(t)) for k in eachindex(d)]
        x=sum(gate.*(K*h(left)));y=sum(gate.*(K*h(right)))
        result+=2width*(x*x+x*y+y*y)/3
    end
    result
end

function _carrier_centered_duties(x,n)
    r=x[1];low=(1-r)/2
    n==3 ? [low+r,low+r*x[2],low] :
        [low+r,low+r*x[2],low+r*x[2]*x[3],low]
end

function _carrier_permutations(v)
    isempty(v) && return [Int[]]
    [[k;p] for k in v for p in _carrier_permutations(filter(!=(k),v))]
end

# A polynomial of degree <=5 in each coordinate is integrated exactly above.
# Interpolate at six rational nodes, then change to tensor Bernstein form.
# Rational arithmetic makes this construction independent of interpolation
# conditioning. Conversion and subdivision round every upper coefficient up.
function _carrier_bernstein_transform()
    q=Rational{BigInt};nodes=q.(0:5)./5
    V=[x^k for x in nodes,k in 0:5]
    B=[k<=i ? q(binomial(i,k),binomial(5,k)) : zero(q) for i in 0:5,k in 0:5]
    B*inv(V)
end

function _carrier_tensor_transform(A,T)
    B=similar(A)
    for axis in 1:ndims(A)
        for index in CartesianIndices(A)
            B[index]=sum(T[index[axis],k]*A[CartesianIndex(ntuple(
                j->j==axis ? k : index[j],ndims(A)))] for k in 1:6)
        end
        A,B=B,A
    end
    A
end

function _carrier_float_upper(x::Rational)
    y=Float64(x)
    Rational{BigInt}(y)<x ? nextfloat(y) : y
end

function _carrier_bernstein_split(A,axis)
    left=similar(A);right=similar(A)
    for index in CartesianIndices(A)
        index[axis]==1 || continue
        at(k)=CartesianIndex(ntuple(j->j==axis ? k : index[j],ndims(A)))
        v=[A[at(k)] for k in 1:6]
        left[at(1)]=v[1];right[at(6)]=v[6]
        for level in 1:5
            v=[nextfloat((v[k]+v[k+1])/2) for k in 1:length(v)-1]
            left[at(level+1)]=first(v);right[at(6-level)]=last(v)
        end
    end
    left,right
end

function _carrier_ac_variance_bound(K,m_max;max_boxes=20000)
    n=size(K,1);n in (3,4) || throw(ArgumentError("three/four physical legs required"))
    q=Rational{BigInt};T=_carrier_bernstein_transform();shape=ntuple(_->6,n-1)
    boxes=Tuple{Array{Float64},Int}[];sample_max=zero(q)
    for order in _carrier_permutations(collect(1:n))
        ordered=q.(K[order,order])
        A=Array{q}(undef,shape)
        for index in CartesianIndices(A)
            x=q.(collect(Tuple(index)).-1)./5;x[1]*=q(m_max)
            d=_carrier_centered_duties(x,n)
            A[index]=_carrier_ac_bridge_variance(d,ordered)
        end
        sample_max=max(sample_max,maximum(A))
        push!(boxes,(_carrier_float_upper.(_carrier_tensor_transform(A,T)),0))
    end
    # Accuracy affects conservatism only: every intermediate maximum is a
    # valid bound. A budget stop therefore never promotes a sampled maximum.
    target=_carrier_float_upper(5sample_max/4)
    subdivisions=0
    while subdivisions<max_boxes
        upper,index=findmax([maximum(b[1]) for b in boxes])
        upper<=target && return (variance=upper,subdivisions=subdivisions)
        A,depth=boxes[index];axis=mod(depth,n-1)+1
        left,right=_carrier_bernstein_split(A,axis)
        boxes[index]=(left,depth+1);push!(boxes,(right,depth+1));subdivisions+=1
    end
    (variance=maximum(maximum(b[1]) for b in boxes),subdivisions=subdivisions)
end

const _CARRIER_THERMAL_CACHE=Dict{Any,Any}()
const _CARRIER_THERMAL_LOCK=ReentrantLock()

function _carrier_thermal_bound(inv)
    geometry=_switching_envelope_geometry(inv)
    K=inv.topology===:THREE_LEG ? geometry.gate_map :
        [geometry.gate_map;geometry.neutral_gate_map']
    key=(Tuple(vec(K)),inv.m_max,inv.f_sw,geometry.damping_rate)
    lock(_CARRIER_THERMAL_LOCK) do
        get!(_CARRIER_THERMAL_CACHE,key) do
            proof=_carrier_ac_variance_bound(K,inv.m_max)
            # For each periodic passive R/L output, the difference from its
            # lossless primitive is <= alpha/omega_sw times that primitive
            # in L2. Since sum(i_leg)=0, g'i=(g-1/2)'i. Mean subtraction is
            # an L2 contraction, giving this separate damping residual.
            residual=geometry.damping_rate/(2pi*inv.f_sw)*sum(abs,K)/(16inv.f_sw)
            (variance_per_V2=proof.variance/inv.f_sw^2,
                damping_rms_per_V=residual,damping_ratio=geometry.damping_rate/(2pi*inv.f_sw),
                subdivisions=proof.subdivisions)
        end
    end
end

function _allocation_coupled_pwm_variance(inv,frozen_variance,voltage_upper)
    inv.pwm_strategy===:NONE && return frozen_variance
    bound=_carrier_thermal_bound(inv)
    lossless=frozen_variance+bound.variance_per_V2*voltage_upper^2
    iszero(bound.damping_ratio) && return lossless
    # Young's inequality avoids a square root of a temporarily negative
    # frozen-current expression at infeasible NLP iterates. No clipping or
    # hidden smoothing floor is introduced into the physical constraint.
    eta=bound.damping_ratio
    (1+eta)*lossless+(1+1/eta)*bound.damping_rms_per_V^2*voltage_upper^2
end
