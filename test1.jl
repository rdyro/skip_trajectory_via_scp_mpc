using PyPlot
include("mpc.jl")

const rho0 = 1.225 # Sea level air dens - kg/m^3
const beta = 1 / 8000 # Scale height - 1 / m
const Cd0 = 0.4 # Vehicle Property
const K = 1.11 # Vehicle Property
const A = 0.754 # Vehicle Property - Surface Area
const m = 1.0 # Vehicle Property - vehicle quality, not given in paper
const g = 9.81 # Gravity
const dt = 1.0
const R = 6371000

function f(t, x, u)
  xn = copy(x)
  xn[1] += x[2] - 0.6 / 2.0 * x[2]^2
  xn[2] += 0.3 * x[2] + u[]
  return xn
end

function Alin(t, x, u)
  return [1.0 1.0 - 0.6 * x[2]; 0.0 1.3]
  #return [1.1 1.0; 0.0 1.1]
end

function Blin(t, x, u)
  return [0.0; 1.0]
end

function veh_Alin(t, x, u)
    A = 0.754 # Vehicle Property - Surface Area
    m = 1.0 # Vehicle Property - vehicle quality, not given in paper
    g = 9.81 # Gravity
    R = 6371000
    dt = 1.0

    r = R + x[3] 
    Cd = get_Cd(u)
    rho = dens(x[3])
    a11 = (-Cd * rho * A * x[1] / m)
    a12 = (-g * cos(x[2])) 
    a21 = (u[1] * rho * A / (2 * m) + ((1 / r) + (g / x[2]^2)) * cos(x[2]))
    a23 = -x[1] * cos(x[2]) / r ^ 2
    a22 = (((g / x[1]) - (x[1] / r)) * cos(x[2]))
    a31 = (sin(x[2]))
    a32 = (x[1] * cos(x[2]))

    Amat = [a11 a12 0;
            a21 a22 a23;
            a31 a32 0]
    
    Amat = eye(3) + dt * Amat


    return Amat
end

function veh_Alin2(t, x, u)
  A = zeros(3, 3)

  A[1, 1] = 1 - (A * dt * rho0 * (1.1100000000000002e+0 * u^2 + 3.9999999999999997e-1) * x1 * exp(-beta * x3)) / m
  A[1, 2] = -g * cos(x2)
  A[1, 3] = ((A * beta * dt * rho0 * (1.1100000000000002e+0 * u^2 + 3.9999999999999997e-1) * x1^2 * exp(-beta * x3)) / m) / 2.e+0
  A[2, 1] = ((A * dt * rho0 * u * exp(-beta * x3)) / m) / 2.e+0 + cos(x2) * (1 / (x3 + R) + g / x1^2)
  A[2, 2] = 1 - sin(x2) * (x1 / (x3 + R) - g / x1)
  A[2, 3] = (-((A * beta * dt * rho0 * u * x1 * exp(-beta * x3)) / m) / 2.e+0) - (x1 * cos(x2)) / (x3 + R)^2
  A[3, 1] = dt * sin(x2)
  A[3, 2] = dt * x1 * cos(x2)
  A[3, 3] = 1

  return A
end

function veh_Blin(t, x, u)
    A = 0.754 # Vehicle Property - Surface Area
    r = 0.4899 # Vehicle Property - vehicle radius, not given in paper
    m = 1.0 # Vehicle Property - vehicle quality, not given in paper
    g = 9.81 # Gravity
    K = 1.11 # Vehicle Property
    dt = 1.0

    rho = dens(x[3])
    B = [(K * rho * x[1]^2 * A * u[1] / m);
         (rho * x[1] * A / (2 * m));
         0]
    
    return dt * B
end

function dens(h)
    rho0 = 1.225 # Sea level air dens - kg/m^3
    beta = 1 / 8000 # Scale height - 1 / m
    return rho0 * exp(-beta * h) 
end

function get_Cd(Cl)
    Cd0 = 0.4 # Vehicle Property
    K = 1.11 # Vehicle Property
    return Cd0 + K * Cl[1] ^ 2
end

function vehf(t, x, u)
    rho0 = 1.225 # Sea level air dens - kg/m^3
    beta = 1 / 8000 # Scale height - 1 / m
    Cd0 = 0.4 # Vehicle Property
    K = 1.11 # Vehicle Property
    A = 0.754 # Vehicle Property - Surface Area
    m = 1.0 # Vehicle Property - vehicle quality, not given in paper
    g = 9.81 # Gravity
    dt = 1.0
    R = 6371000
    # x = [v, th, h]

    xn = copy(x)
    rho = dens(x[3])
    Cd = get_Cd(u)
    r = R + x[3] 

    @show(Cd)

    xn[1] += dt * -(Cd * rho * x[1] ^ 2 * A) / 2 * m - g * sin(x[2])
    xn[2] += dt * (u[1] * rho * x[1] * A) / 2 * m + (x[1] / r - g / x[1]) * cos(x[2])
    xn[3] += dt * x[1] * sin(x[2])
    return xn
end

function veh_test()
    # Following the problem setup in the paper
    x0 = [11000.0; -0.13962; 110000];

    Q = 1.0 * speye(3)
    #Q[2, 2] = 0.0 # Quadratic costs on velocity and height
    #Q[1, 1] = 10.0 # Penalize high velocities
    #Q[3, 3] = 1.0 # Penalize losing altitude
    R = 0.0 * speye(1)
    P = 1e1 * Q
    N = 100

    Clmax = 1.5
    ub = [[-Clmax], [Clmax]]
    xb = nothing # we may need to include some terminal velocity constraints in this

    #(Xplan, Uplan) = scpMPC(vehf, veh_Alin, veh_Blin, Q, R, P, x0, N, ub=ub, xb=xb)

    x = x0
    Xactual = [x]
    Uactual = []
    Xold = nothing
    Uold = nothing
    for i in 1:30
        #(X, U) = scpMPC(vehf, veh_Alin, veh_Blin, Q, R, P, x, N, ub=ub, xb=xb,
        #                Xguess=Xold, Uguess=Uold)
        #(X, U) = scpMPC(f, Alin, Blin, Q, R, P, x, N, ub=ub)
        #u = [U[1] + 0.1 * randn()]
        u = 0.5
        push!(Uactual, u)
        x = vehf(i, x, u) #+ 0.01 * randn(3)
        push!(Xactual, x)
        println(i)
    
        #Xold = [X[4:end]; zeros(3)]
        #Uold = [U[2:end]; zeros(1)]
    end
    
    Xactual = vcat(Xactual...)
    Uactual = vcat(Uactual...)
   
    figure(1)
    clf()
    plot(Xactual[1:3:end], label="\$x_{1,a}\$", color="red")
    plot(Xactual[2:3:end], label="\$x_{2,a}\$", color="blue")
    plot(Xactual[3:3:end], label="\$x_{3,a}\$", color="green")
    plot(Uactual, label="\$u_a\$", color="black")
    
    #plot(Xplan[1:2:30*2], label="\$x_{1,a}\$", color="red", linestyle="--")
    #plot(Xplan[2:2:30*2], label="\$x_{2,a}\$", color="blue", linestyle="--")
    #plot(Uplan[1:30], label="\$u_a\$", color="black", linestyle="--")
    legend()
    return
end


function solver_test()
  x0 = [3.0; 3.0]

  A = [1.0 1.0;
       0.0 1.0]
  B = [0.0;
       1.0][:, :]
  Q = 1.0 * speye(2)
  R = 1.0 * speye(1)
  P = 1e1 * Q
  N = 30
  ub = [[-0.7], [0.7]]
  xb = [[-10.0, ], [10.0]]
  #ub = nothing

  (Xlin, Ulin) = linMPC(A, B, Q, R, P, x0, N, ub=ub)
  #@btime $fa()
  #@time fa()

  #fb() = scpMPC2(f, Alin, Blin, Q, R, P, x0, N, ub=ub)
  #@btime $fb()
  #@time fb()

  f2(t, x, u) = A * x + B * u
  Alin2(t, x, u) = A
  Blin2(t, x, u) = B[:, :]

  #@time (X, U) = scpMPC2(f, Alin, Blin, Q, R, P, x0, N, ub=ub)
  (X, U) = scpMPC2(f2, Alin2, Blin2, Q, R, P, x0, N, ub=ub)
  #X += 0.1 * randn(size(X))
  #U += 0.1 * randn(size(U))
  fc() = scpMPC2(f, Alin, Blin, Q, R, P, x0, N, ub=ub, Xguess=X, Uguess=U)
  #@btime $fc()
  @time fc()
  #(X, U) = fc()

  COLORS = ["red", "blue", "green"]
  figure(2)
  clf()
  plot(X[1:2:end], label="x_1", color=COLORS[1])
  plot(X[2:2:end], label="x_2", color=COLORS[2])
  plot(U, label="u", color=COLORS[3])
  plot(Xlin[1:2:end], label="lin x_1", color=COLORS[1], linestyle="--")
  plot(Xlin[2:2:end], label="lin x_2", color=COLORS[2], linestyle="--")
  plot(Ulin, label="lin u", color=COLORS[3], linestyle="--")
  legend()

  simulate(x0)
  return
end

function simulate(x0)
  Q = 1.0 * eye(2)
  R = 1.0 * eye(1)
  P = 1e1 * Q
  N = 30
  ub = [[-1.0], [1.0]]
  xb = [[-50.0, -50.0], [50.0, 50.0]]

  (Xplan, Uplan) = scpMPC2(f, Alin, Blin, Q, R, P, x0, N, ub=ub, xb=xb)

  x = x0
  Xactual = [x]
  Uactual = []
  Xold = nothing
  Uold = nothing
  for i in 1:20
    println(i)
    (X, U) = scpMPC2(f, Alin, Blin, Q, R, P, x, N, ub=ub, xb=xb, Xguess=Xold,
                    Uguess=Uold)
    #(X, U) = scpMPC(f, Alin, Blin, Q, R, P, x, N, ub=ub)
    u = [U[1] + 0.0 * randn()]
    push!(Uactual, u)
    x = f(i, x, u) + 0.0 * randn(2)
    push!(Xactual, x)

    Xold = [X[3:end]; zeros(2)]
    Uold = [U[2:end]; zeros(1)]
  end

  Xactual = vcat(Xactual...)
  Uactual = vcat(Uactual...)

  figure(1)
  clf()
  plot(Xactual[1:2:end], label="\$x_{1,a}\$", color="red")
  plot(Xactual[2:2:end], label="\$x_{2,a}\$", color="blue")
  
  plot(Uactual, label="\$u_a\$", color="black")

  plot(Xplan[1:2:30*2], label="\$x_{1,p}\$", color="red", linestyle="--")
  plot(Xplan[2:2:30*2], label="\$x_{2,p}\$", color="blue", linestyle="--")
  plot(Uplan[1:30], label="\$u_p\$", color="black", linestyle="--")
  legend()
  return
end
