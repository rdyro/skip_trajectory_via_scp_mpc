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
