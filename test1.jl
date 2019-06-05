using PyPlot
include("mpc.jl")

const Ms = 1e-5
const Area = 0.754 * Ms^2
const m = 1000.0 * Ms # Vehicle Property - vehicle quality, not given in paper
const g = 9.81 * Ms # Gravity
const K = 1.11 # Vehicle Property
const dt = 1.0
const R = 6371000 * Ms
const Cd0 = 0.8 # Vehicle Property
const K = 1.11 # Vehicle Property
const rho0 = 1.225 * Ms / Ms^3 # Sea level air dens - kg/m^3
const beta = 1 / (8000 * Ms) # Scale height - 1 / m
const eps = 1e-7

function veh_Alin(t, x, u)
    r = R + x[3] 
    Cd = get_Cd(u)
    rho = dens(x[3])

    a11 = -Cd * rho * Area * x[1] / m
    a12 = -g * cos(x[2]) 
    a21 = u[1] * rho * Area / (2 * m) + ((1 / r) + (g / (x[1] + sqrt(eps))^2)) * cos(x[2])
    a22 = ((g / (x[1] + eps)) - (x[1] / r)) * sin(x[2])
    a23 = -x[1] * cos(x[2]) / (r^2)
    a31 = sin(x[2])
    a32 = x[1] * cos(x[2])

    Amat = [a11 a12 0;
            a21 a22 a23;
            a31 a32 0]
    
    Amat = eye(3) + dt * Amat


    return Amat
end

function veh_Alin3(t, x, u)
  u = u[1]
  x1 = x[1]
  x2 = x[2]
  x3 = x[3]
  Amat = zeros(3, 3)

  Amat[1,1] = (1 - (Area * dt * rho0 * (1.1100000000000002e+0 * u^2 + 3.9999999999999997e-1) * x1 * exp(-beta * x3)) / m)
  Amat[1,2] = -g * cos(x2)
  Amat[1,3] = (((Area * beta * dt * rho0 * (1.1100000000000002e+0 * u^2 + 3.9999999999999997e-1) * x1^2 * exp(-beta * x3)) / m) / 2.e+0)
  Amat[2,1] = (((Area * dt * rho0 * u * exp(-beta * x3)) / m) / 2.e+0 + cos(x2) * (1 / (x3 + R) + g / x1^2))
  Amat[2,2] = 1 - sin(x2) * (x1 / (x3 + R) - g / x1)
  Amat[2,3] = (-((Area * beta * dt * rho0 * u * x1 * exp(-beta * x3)) / m) / 2.e+0) - (x1 * cos(x2)) / (x3 + R)^2
  Amat[3,1] = dt * sin(x2)
  Amat[3,2] = dt * x1 * cos(x2)
  Amat[3,3] = 1

  return Amat
end


function veh_Blin(t, x, u)
    rho = dens(x[3])
    B = [(K * rho * x[1]^2 * Area * u[1] / m);
         (rho * x[1] * Area / (2 * m));
         0]
    
    return dt .* B
end

function dens(h)
   return rho0 * exp(-beta * h) 
end

function get_Cd(Cl)
   return Cd0 + K * Cl[1] ^ 2
end

function vehf(t, x, u)

   # x = [v, th, h]


    xn = copy(x)
    rho = dens(x[3])
    Cd = get_Cd(u)
    r = R + x[3] 

    xn[1] += dt * -(Cd * rho * x[1] ^ 2 * Area) / (2 * m) - g * sin(x[2])
    xn[2] += dt * (u[1] * rho * x[1] * Area) / (2 * m) + (x[1] / r - g / (x[1] + eps)) * cos(x[2])
    xn[3] += dt * x[1] * sin(x[2])
    return xn
end

function veh_test()
    # Following the problem setup in the paper
    x0 = [11000.0 * Ms; -0.13962; 110000 * Ms];

    Q = 1.0 * speye(3)
    Q[2, 2] = 0.0 # Quadratic costs on velocity and height
    Q[1, 1] = 1.0 # Penalize high velocities
    Q[3, 3] = 0.0 # Penalize losing altitude
    R = 0.0 * speye(1)
    P = 1e1 * Q
    N = 2

    Clmax = 1.5
    ub = [[-Clmax], [Clmax]]
    xb = nothing # we may need to include some terminal velocity constraints in this

    (Xplan, Uplan) = scpMPC(vehf, veh_Alin, veh_Blin, Q, R, P, x0, N, ub=ub, xb=xb)
    #@show(Xplan, Uplan)

    x = x0
    Xactual = [x]
    Uactual = []
    Xold = nothing
    Uold = nothing
    for i in 1:4000
        #(X, U) = scpMPC(vehf, veh_Alin, veh_Blin, Q, R, P, x, N, ub=ub, xb=xb,
        #                Xguess=Xold, Uguess=Uold)
        #(X, U) = scpMPC(f, Alin, Blin, Q, R, P, x, N, ub=ub)
        #u = [U[1] + 0.1 * randn()]
        
        u = 1.5
        push!(Uactual, u)
        x = vehf(i, x, u) #+ 0.01 * randn(3)
        push!(Xactual, x)
        #println(i)
    
        #Xold = [X[4:end]; zeros(3)]
        #Uold = [U[2:end]; zeros(1)]

        if x[3] <= 5000.0 * Ms
            break
        end
    end
    
    Xactual = vcat(Xactual...)
    Uactual = vcat(Uactual...)
   
    figure(1)
    clf()
    plot(Xactual[1:3:end], label="\$x_{1,a}\$", color="red")
    plot(Xactual[2:3:end], label="\$x_{2,a}\$", color="blue")
    plot(Xactual[3:3:end], label="\$x_{3,a}\$", color="green")
    plot(Uactual, label="\$u_a\$", color="black")
    
    #plot(Xplan[3:3:30*3], label="\$x_{3,a}\$", color="green", linestyle="--")
    #plot(Xplan[2:2:30*2], label="\$x_{2,a}\$", color="blue", linestyle="--")
    #plot(Uplan[1:30], label="\$u_a\$", color="black", linestyle="--")
    legend()
    return
end


