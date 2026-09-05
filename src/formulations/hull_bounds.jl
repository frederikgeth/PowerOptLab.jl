# Rational arithmetic treats the supplied Float64 curve/domain data exactly.
# Returned gap bounds are rounded outwards; envelope/witness values are displays.
_hull_rational(x) = Rational{BigInt}(x)
function _hull_endpoint(f,x)
    x<=first(f.breakpoints) && return _hull_rational(first(f.values))
    x>=last(f.breakpoints) && return _hull_rational(last(f.values))
    j=findlast(k->k<=x,f.breakpoints)
    a,b=_hull_rational.(f.breakpoints[j:j+1])
    u,v=_hull_rational.(f.values[j:j+1])
    u+(v-u)*(_hull_rational(x)-a)/(b-a)
end
function _hull_chain(xs,ys,upper)
    chain=Int[]
    for j in eachindex(xs)
        while length(chain)>=2
            a,b=chain[end-1:end]
            left=(ys[b]-ys[a])*(xs[j]-xs[b])
            right=(ys[j]-ys[b])*(xs[b]-xs[a])
            (upper ? left<=right : left>=right) || break
            pop!(chain)
        end
        push!(chain,j)
    end
    chain
end
function _hull_max_gap(xs,ys,chain,upper)
    segment=1
    largest=zero(first(ys)); witness=1; envelope=first(ys)
    for j in eachindex(xs)
        while segment<length(chain)-1 && xs[j]>xs[chain[segment+1]]
            segment+=1
        end
        a,b=chain[segment],chain[min(segment+1,length(chain))]
        value=a==b ? ys[a] : ys[a]+(ys[b]-ys[a])*(xs[j]-xs[a])/(xs[b]-xs[a])
        gap=upper ? value-ys[j] : ys[j]-value
        if gap>largest
            largest,witness,envelope=gap,j,value
        end
    end
    bound=Float64(largest)
    isfinite(bound) || throw(ArgumentError("Hull gap exceeds Float64 range; rescale output units"))
    _hull_rational(bound)<largest && (bound=nextfloat(bound))
    isfinite(bound) || throw(ArgumentError("Hull gap exceeds Float64 range; rescale output units"))
    return bound,(input=Float64(xs[witness]),canonical=Float64(ys[witness]),
        envelope=Float64(envelope))
end

"""
    hull_gap_bound(curve, domain=(first(curve.breakpoints), last(curve.breakpoints)))

Compute maximum vertical deviations of the bounded graph hull from its canonical
PWL curve. `upper_gap` is the concave-envelope excess, `lower_gap` the deficit to
the convex envelope; any hull point obeys `-lower_gap <= y-f(x) <= upper_gap`.
Returns physical units, envelope vertices and attaining witness inputs. A supplied
finite domain can cut segments or extend into the canonical flat tails.

Monotone-chain scans use O(n) rational operations on the exact supplied Float64
coordinates; maximum gaps occur at restricted curve knots. Gap bounds are rounded
outwards to Float64; displayed envelope/witness values use nearest rounding.
These bound vertical scalar graph discrepancy, not a network objective gap,
equilibrium error or DOE safety. Constraint-coefficient rounding and solver
residuals remain additional. No solver or JuMP model is constructed.
"""
function hull_gap_bound(f::PWLFunction,domain=(first(f.breakpoints),last(f.breakpoints)))
    lo,hi=_relation_domain(domain)
    interior=[j for j in eachindex(f.breakpoints) if lo<f.breakpoints[j]<hi]
    xs=_hull_rational.([lo;[f.breakpoints[j] for j in interior];hi])
    ys=[_hull_endpoint(f,lo);[_hull_rational(f.values[j]) for j in interior];_hull_endpoint(f,hi)]
    lo==hi && (pop!(xs);pop!(ys))
    upper,lower=_hull_chain(xs,ys,true),_hull_chain(xs,ys,false)
    ug,uw=_hull_max_gap(xs,ys,upper,true)
    lg,lw=_hull_max_gap(xs,ys,lower,false)
    vertices(chain)=Tuple((Float64(xs[j]),Float64(ys[j])) for j in chain)
    return (domain=(lo,hi),input_unit=f.input_unit,output_unit=f.output_unit,
        upper_gap=ug,lower_gap=lg,error_lower=-lg,error_upper=ug,
        upper_witness=uw,lower_witness=lw,
        upper_envelope=vertices(upper),lower_envelope=vertices(lower))
end
