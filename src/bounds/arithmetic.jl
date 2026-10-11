# Outward binary64 interval arithmetic for the small propagation kernel.
# Rational comparisons make elementary operations rigorous without changing
# process-global MPFR precision/rounding or adding a solver dependency.
const _BRational = Rational{BigInt}
_brat(x::Real) = _BRational(x)

struct _BInterval
    lo::Float64
    hi::Float64
end
_bpoint(x::Real) = _BInterval(Float64(x), Float64(x))
const _BENTIRE = _BInterval(-Inf, Inf)
const _BNONNEG = _BInterval(0.0, Inf)

function _bround(q::_BRational, up::Bool)
    iszero(q) && return 0.0
    q < 0 && return -_bround(-q,!up)
    # Do not seed from Float64(q): Julia's Rational{BigInt} conversion can use
    # the caller's current BigFloat precision. Integer divrem directly computes
    # the adjacent binary64 significands, including subnormals and overflow.
    n=numerator(q); d=denominator(q)
    exponent=ndigits(n;base=2)-ndigits(d;base=2)
    below=exponent>=0 ? n < (d<<exponent) : (n<<(-exponent)) < d
    below && (exponent-=1)
    exponent>1023 && return up ? Inf : floatmax(Float64)
    step=max(exponent-52,-1074)
    significand,remainder=step>=0 ? divrem(n,d<<step) : divrem(n<<(-step),d)
    up && !iszero(remainder) && (significand+=1)
    ldexp(Float64(significand),step)
end

function _baddend(a::Float64, b::Float64, up::Bool)
    isfinite(a) && isfinite(b) && return _bround(_brat(a) + _brat(b), up)
    x = a + b
    isnan(x) ? (up ? Inf : -Inf) : x
end
_badd(a::_BInterval, b::_BInterval) =
    _BInterval(_baddend(a.lo,b.lo,false), _baddend(a.hi,b.hi,true))
_bneg(a::_BInterval) = _BInterval(-a.hi,-a.lo)
_bsub(a::_BInterval, b::_BInterval) = _badd(a,_bneg(b))

function _bmulend(a::Float64, b::Float64, up::Bool)
    (iszero(a) || iszero(b)) && return 0.0
    isfinite(a) && isfinite(b) && return _bround(_brat(a)*_brat(b), up)
    return a*b
end
function _bmul(a::_BInterval, b::_BInterval)
    pairs = ((a.lo,b.lo),(a.lo,b.hi),(a.hi,b.lo),(a.hi,b.hi))
    _BInterval(minimum(_bmulend(x,y,false) for (x,y) in pairs),
               maximum(_bmulend(x,y,true) for (x,y) in pairs))
end
function _bscale(c::_BRational, a::_BInterval)
    iszero(c) && return _bpoint(0)
    c < 0 && return _bneg(_bscale(-c,a))
    _BInterval(isfinite(a.lo) ? _bround(c*_brat(a.lo),false) : a.lo,
               isfinite(a.hi) ? _bround(c*_brat(a.hi),true) : a.hi)
end
function _bdiv(a::_BInterval, b::_BInterval)
    b.lo <= 0 <= b.hi && return _BENTIRE
    b.hi < 0 && return _bdiv(_bneg(a),_bneg(b))
    divide(x,y,up) = isfinite(x) && isfinite(y) ? _bround(_brat(x)/_brat(y),up) :
        isfinite(x) ? 0.0 : x
    _BInterval(divide(a.lo,a.lo>=0 ? b.hi : b.lo,false),
               divide(a.hi,a.hi<=0 ? b.hi : b.lo,true))
end
_bminabs(a::_BInterval) = a.lo <= 0 <= a.hi ? 0.0 : min(abs(a.lo),abs(a.hi))
_bmaxabs(a::_BInterval) = max(abs(a.lo),abs(a.hi))
_bsquare(a::_BInterval) = _BInterval(
    _bmulend(_bminabs(a),_bminabs(a),false),
    _bmulend(_bmaxabs(a),_bmaxabs(a),true))
function _bsqrtend(x::Float64, up::Bool)
    isinf(x) && return x
    y = sqrt(x)
    yy = _brat(y)^2
    up && yy < _brat(x) && return nextfloat(y)
    !up && yy > _brat(x) && return prevfloat(y)
    y
end
_bsqrt(a::_BInterval) = _BInterval(_bsqrtend(max(0.0,a.lo),false),_bsqrtend(max(0.0,a.hi),true))
_bhypot(a::_BInterval,b::_BInterval) = _bsqrt(_badd(_bsquare(a),_bsquare(b)))
_bsum(xs) = foldl(_badd,xs;init=_bpoint(0))
