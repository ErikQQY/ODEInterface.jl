using ODEInterface
@ODEInterface.import_huge

# We look at the boundary value problem
#
#                       sin⁴(x)
#    ε⋅y'(x) = sin²(x)-λ──────  ;   y(-π/2)=y(π/2)=1
#                         y
# 
# Here ε is given (e.g. ε=0.1) and λ is an unknown parameter: λ'=0
#
# We use colsys/colnew and a "homotopy": starting with ε=1.0, using this
# solution as start-guess for ε=0.5, then using this as start guess
# for ε=0.2 and then ε=0.1

a, b = 0.0, 1.0
orders = [1, 1, 1]
ζ = [a, a, b]

global ε = nothing 
global ε_old = nothing
global sol_old = nothing

function rhs(x, z, y, f)
    e = 2.7
    f[1] = (1 + z[2] - sin(x)) * y[1] + cos(x)
    f[2] = cos(x)
    f[3] = y[1]
    f[4] = (z[1] - sin(x)) *(y[1] - e^x)
end

function Drhs(x, z, y, df)
    e = 2.7
    df[:] .= 0.0
    df[1,2] = y[1]
    df[1,4] = 1 + z[2] - sin(x)
    df[3,4] = 1.0
    df[4,1] = y[1] - e^x
    df[4,4] = z[1] - sin(x)
end

function bc(i, z, bc)
    if i == 1
        bc[1] = z[1]
    elseif i == 2
        bc[1] = z[3] - 1.0
    else
        bc[1] = z[2] - sin(1.0)
    end
end

function Dbc(i, z, dbc)
    if i == 1
        dbc[1] = 1.0
        dbc[2] = 0.0
        dbc[3] = 0.0
    elseif i == 2
        dbc[1] = 0.0
        dbc[2] = 0.0
        dbc[3] = 1.0
    else
        dbc[1] = 0.0
        dbc[2] = 1.0
        dbc[3] = 0.0
    end
end

opt = OptionsODE("example 8",
      OPT_BVPCLASS => 2, OPT_COLLOCATIONPTS => 7,
      OPT_RTOL => [1e-4, 1e-4, 1e-4], OPT_MAXSUBINTERVALS => 200)
      
sol, retcode, stats = coldae([a,b], orders, 1, ζ, rhs, Drhs, bc, Dbc, nothing ,opt);
@assert retcode>0
      
zz, yy = evalSolution(sol, xx)
xx = collect(LinRange(a, b, 100))
println("sol z: ", zz)
println("sol y: ", yy)

# vim:syn=julia:cc=79:fdm=indent:
