# Export values from the actual Julia primitive implementation for plotting.
using PowerOptLab
path=isempty(ARGS) ? joinpath(@__DIR__,"figure_data.csv") : only(ARGS)
soft=SoftplusFormulation(.05/log(2))
local_c2=LocalC2Formulation(.05/(3/16))
open(path,"w") do io
    println(io,"z,exact,soft,c2,soft_error,c2_error,soft_curvature,c2_curvature")
    for z in range(-.8,.8;length=801)
        exact=max(z,0.)
        s,c=hinge_value(z,soft),hinge_value(z,local_c2)
        println(io,join((z,exact,s,c,s-exact,c-exact,
            hinge_derivatives(z,soft)[2],hinge_derivatives(z,local_c2)[2]),","))
    end
end
