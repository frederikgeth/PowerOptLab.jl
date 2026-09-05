module PowerOptLabPiecewiseLinearOptExt
using PowerOptLab
using PiecewiseLinearOpt
using JuMP
function PowerOptLab._exact_pwl_graph!(m::JuMP.Model,x,xs,ys,r::PowerOptLab.ExactPWLGraph)
    return PiecewiseLinearOpt.piecewiselinear(m,x,xs,ys; method=r.method)
end
end
