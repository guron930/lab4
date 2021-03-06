using DifferentialEquations
using Gadfly
using Plots

function portret(w, g, x0, y0)

    function SDU(du,u,p,t)
        du[1] = u[2]
        du[2] = -w*w*u[1]-g*u[2]-f(t)
    end

    u0 = [x0, y0] 
    tspan = (0.0, 25)

    prob = ODEProblem(SDU, u0, tspan) 
    sol = solve(prob, RK4(),reltol=1e-6, timeseries_steps = 0.05)

    N = length(sol.u)
    J = length(sol.u[1])

    U = zeros(N, J)

    for i in 1:N, j in 1:J
        U[i,j] = sol.u[i][j] 
    end
    U
end

f(t) = 0
ans1 = portret(7, 0,1, 1);


Plots.plot(ans1)

set_default_plot_size(30cm, 20cm)
 Gadfly.plot(x = ans1[:,1], y = ans1[:,2],
        Guide.title("Колебания без затухания без действия внешней силы"))

f(t) = 0
ans2 = portret(6, 2, -1, 1)
set_default_plot_size(40cm, 20cm)
 Gadfly.plot(x = ans2[:,1], y = ans2[:,2],
        Guide.title("Колебания c затуханием и без действия внешней силы"))

Plots.plot(ans2)

f(t) = cos(3t)
ans3 = portret(1, 5, -1, 1)
set_default_plot_size(40cm, 20cm)
 Gadfly.plot(x = ans3[:,1], y = ans3[:,2],
        Guide.title("Колебания c затуханием и под действием внешней силы"))

Plots.plot(ans3)


