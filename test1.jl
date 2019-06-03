using PyPlot

function Alin(t, x)
  return [1.0 1.0 - 0.4 * x[2]; 0.0 1.1]
  #return [1.1 1.0; 0.0 1.1]
end

function Blin(t, x, u)
  return [0.0; 1.0]
end

function simulate(Alin, Blin, A, B, Q, R, P, x0, N, ub)
    # This is the first step no guess
    (X, U) = scpMPC(Alin, Blin, Q, R, P, x0, N, ub=ub)

    fc() = scpMPC(Alin, Blin, Q, R, P, x0, N, ub=ub, Xguess=X, Uguess=U)

    # at each step:
    # x0 = current state
    # solve the scpMPC
    # propogate dynamics with u and noise
    # shift the X,U for next guess 
    # 

    return 
end

function main()
  x0 = [3.0; 3.0]
  A = [1.0 0.1;
       0.0 1.0]
  B = [0.0;
       1.0]
  Q = 1.0 * speye(2)
  R = 1.0 * speye(1)
  P = 1e1 * Q
  N = 100
  ub = [[-0.05], [0.05]]
  #ub = nothing

  fa() = linMPC(A, B, Q, R, P, x0, N, ub=ub)
  #@btime $fa()
  @time fa()
  (X, U) = linMPC(A, B, Q, R, P, x0, N, ub=ub)
  clf()
  plot(X[1:2:end])
  plot(X[2:2:end])
  plot(U)
  return

  fb() = scpMPC(Alin, Blin, Q, R, P, x0, N, ub=ub)
  #@btime $fb()
  @time fb()

  (X, U) = scpMPC(Alin, Blin, Q, R, P, x0, N, ub=ub)
  fc() = scpMPC(Alin, Blin, Q, R, P, x0, N, ub=ub, Xguess=X, Uguess=U)
  #@btime $fc()
  @time fc()



  clf()
  plot(X[1:2:end])
  plot(X[2:2:end])
  plot(U)
  return
end
