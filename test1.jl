using PyPlot
include("mpc.jl")

function Alin(t, x)
  return [1.0 1.0 - 0.6 * x[2]; 0.0 1.1]
  #return [1.1 1.0; 0.0 1.1]
end

function Blin(t, x, u)
  return [0.0; 1.0]
end

function main()
  x0 = [3.0; 3.0]
  A = [1.0 1.0;
       0.0 1.0]
  B = [0.0;
       1.0]
  Q = 1.0 * speye(2)
  R = 1.0 * speye(1)
  P = 1e1 * Q
  N = 100
  ub = [[-2.0], [2.0]]
  xb = [[-10.0, ], [10.0]]
  #ub = nothing

  fa() = linMPC(A, B, Q, R, P, x0, N, ub=ub)
  #@btime $fa()
  @time fa()

  fb() = scpMPC(Alin, Blin, Q, R, P, x0, N, ub=ub)
  #@btime $fb()
  @time fb()

  (X, U) = scpMPC(Alin, Blin, Q, R, P, x0, N, ub=ub)
  X += 0.1 * randn(size(X))
  U += 0.1 * randn(size(U))
  fc() = scpMPC(Alin, Blin, Q, R, P, x0, N, ub=ub, Xguess=X, Uguess=U)
  #@btime $fc()
  @time fc()
  (X, U) = fc()

  clf()
  plot(X[1:2:end])
  plot(X[2:2:end])
  plot(U)
  return
end
