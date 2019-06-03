using PyPlot
include("mpc.jl")

function f(t, x, u)
  xn = copy(x)
  xn[1] += x[2] - 0.3 * x[2]^2
  xn[2] += 0.1 * x[2] + u[]
  return xn
end

function Alin(t, x)
  return [1.0 1.0 - 0.6 * x[2]; 0.0 1.1]
  #return [1.1 1.0; 0.0 1.1]
end

function Blin(t, x, u)
  return [0.0; 1.0]
end

function veh_Alin(t, x, u)
    A = 0.754 # Vehicle Property - Surface Area
    r = 0.4899 # Vehicle Property - vehicle radius, not given in paper
    m = 1.0 # Vehicle Property - vehicle quality, not given in paper
    g = 9.81 # Gravity

    Cd = get_Cd(u)
    rho = dens(x[3])
    a11 = (-Cd * rho * A * x[1] / m)
    a12 = (-g * cos(x[2])) 
    a21 = (u[1] * rho * A / (2 * m) + ((1 / r) + (g / x[2]^2)) * cos(x[2]))
    a22 = (((g / x[1]) - (x[1] / r)) * cos(x[2]))
    a31 = (sin(x[2]))
    a32 = (x[1] * cos(x[2]))

    Amat = [a11 a12 0;
            a21 a22 0;
            a31 a32 0]

    return Amat
end

function veh_Blin(t, x, u)
    A = 0.754 # Vehicle Property - Surface Area
    r = 0.4899 # Vehicle Property - vehicle radius, not given in paper
    m = 1.0 # Vehicle Property - vehicle quality, not given in paper
    g = 9.81 # Gravity
    K = 1.11 # Vehicle Property

    rho = dens(x[3])
    B = [(K * rho * x[1]^2 * A * u[1] / m);
         (rho * x[1] * A / (2 * m));
         0]
    
    return B
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
    A = 0.754 # Vehicle Property - Surface Area
    r = 0.4899 # Vehicle Property - vehicle radius, not given in paper
    m = 1.0 # Vehicle Property - vehicle quality, not given in paper
    g = 9.81 # Gravity
    # x = [v, th, h]

    xn = copy(x)
    rho = dens(x[3])
    Cd = get_Cd(u)

    xn[1] += -(Cd * rho * x[1] ^ 2 * A) / 2 * m - g * sin(x[2])
    xn[2] += (u[1] * rho * x[1] * A) / 2 * m + (x[1] / r - g / x[1]) * cos(x[2])
    xn[3] += x[1] * sin(x[2])
    return xn
end

function veh_test()
    # Following the problem setup in the paper
    x0 = [11000.0; -0.13962; 110000];

    Q = 1.0 * speye(3)
    Q[2, 2] = 0.0 # Quadratic costs on velocity and height
    Q[1, 1] = 10.0 # Penalize high velocities
    Q[3, 3] = 1.0 # Penalize losing altitude
    R = 0.0 * speye(1)
    P = 1e1 * Q
    N = 100

    Clmax = 1.5
    ub = [[-Clmax], [Clmax]]
    xb = nothing # we may need to include some terminal velocity constraints in this

    (Xplan, Uplan) = scpMPC(vehf, veh_Alin, veh_Blin, Q, R, P, x0, N, ub=ub, xb=xb)

    x = x0
    Xactual = [x]
    Uactual = []
    Xold = nothing
    Uold = nothing
    for i in 1:30
        (X, U) = scpMPC(vehf, veh_Alin, veh_Blin, Q, R, P, x, N, ub=ub, xb=xb,
                        Xguess=Xold, Uguess=Uold)
        #(X, U) = scpMPC(f, Alin, Blin, Q, R, P, x, N, ub=ub)
        u = [U[1] + 0.1 * randn()]
        push!(Uactual, u)
        x = vehf(i, x, u) + 0.1 * randn(3)
        push!(Xactual, x)
        println(i)
    
        Xold = [X[4:end]; zeros(3)]
        Uold = [U[2:end]; zeros(1)]
    end
    
    Xactual = vcat(Xactual...)
    Uactual = vcat(Uactual...)
   
    figure(1)
    clf()
    plot(Xactual[1:2:end], label="\$x_{1,a}\$", color="red")
    plot(Xactual[2:2:end], label="\$x_{2,a}\$", color="blue")
    plot(Uactual, label="\$u_a\$", color="black")
    
    plot(Xplan[1:2:30*2], label="\$x_{1,a}\$", color="red", linestyle="--")
    plot(Xplan[2:2:30*2], label="\$x_{2,a}\$", color="blue", linestyle="--")
    plot(Uplan[1:30], label="\$u_a\$", color="black", linestyle="--")
    legend()
    return
end


function solver_test()
  x0 = [3.0; 3.0]

  A = [1.0 1.0;
       0.0 1.0]
  B = [0.0;
       1.0]
  Q = 1.0 * speye(2)
  R = 1.0 * speye(1)
  P = 1e1 * Q
  N = 100
  ub = [[-0.5], [0.5]]
  xb = [[-10.0, ], [10.0]]
  #ub = nothing

  fa() = linMPC(A, B, Q, R, P, x0, N, ub=ub)
  #@btime $fa()
  @time fa()

  fb() = scpMPC(f, Alin, Blin, Q, R, P, x0, N, ub=ub)
  #@btime $fb()
  @time fb()

  (X, U) = scpMPC(f, Alin, Blin, Q, R, P, x0, N, ub=ub)
  X += 0.1 * randn(size(X))
  U += 0.1 * randn(size(U))
  fc() = scpMPC(f, Alin, Blin, Q, R, P, x0, N, ub=ub, Xguess=X, Uguess=U)
  #@btime $fc()
  @time fc()
  (X, U) = fc()

  clf()
  plot(X[1:2:end])
  plot(X[2:2:end])
  plot(U)

  simulate(x0)
  return
end

function simulate(x0)
  Q = 1.0 * eye(2)
  R = 1.0 * eye(1)
  P = 1e1 * Q
  N = 100
  ub = [[-0.5], [0.5]]
  xb = [[-50.0, -50.0], [50.0, 50.0]]

  (Xplan, Uplan) = scpMPC(f, Alin, Blin, Q, R, P, x0, N, ub=ub, xb=xb)

  x = x0
  Xactual = [x]
  Uactual = []
  Xold = nothing
  Uold = nothing
  for i in 1:30
    (X, U) = scpMPC(f, Alin, Blin, Q, R, P, x, N, ub=ub, xb=xb, Xguess=Xold,
                    Uguess=Uold)
    #(X, U) = scpMPC(f, Alin, Blin, Q, R, P, x, N, ub=ub)
    u = [U[1] + 0.1 * randn()]
    push!(Uactual, u)
    x = f(i, x, u) + 0.1 * randn(2)
    push!(Xactual, x)
    println(i)

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

  plot(Xplan[1:2:30*2], label="\$x_{1,a}\$", color="red", linestyle="--")
  plot(Xplan[2:2:30*2], label="\$x_{2,a}\$", color="blue", linestyle="--")
  plot(Uplan[1:30], label="\$u_a\$", color="black", linestyle="--")
  legend()
  return
end
