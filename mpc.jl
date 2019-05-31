using Convex
#using Mosek, Gurobi
using ECOS
using LinearAlgebra, SparseArrays
using BenchmarkTools, Suppressor

eye(k) = Matrix{Float64}(I, k, k)
speye(k) = sparse(I, k, k)
function linMPC(A, B, Q, R, P, x0, N; xb=nothing, ub=nothing, Xref=nothing,
                Uref=nothing, Xguess=nothing, Uguess=nothing)
  #solver = MosekSolver()
  #solver = GurobiSolver()
  solver = ECOSSolver(maxit=1e3, eps=1e-7, verbose=false)
  # create variables and reference trajectories ###############################
  xdim = size(Q, 1)
  udim = size(R, 1)

  Xref = Xref == nothing ? [x0; zeros(xdim * N)] : Xref
  Uref = Uref == nothing ? zeros(udim * N) : Uref
  X = Variable((N + 1) * xdim)
  U = Variable(N * udim)
  X.value = Xguess == nothing ? Xref : Xguess
  U.value = Uguess == nothing ? Uref : Uguess

  # enforce system dynamics ###################################################
  # building a large block diagonal matrix is faster than mapping over cstrs
  Ap = Variable(xdim * N, xdim * N)
  Bp = Variable(xdim * N, udim * N)
  Aa = spzeros(xdim * N, xdim * N)
  Ba = spzeros(xdim * N, udim * N)
  for i in 1:N
    Aa[(xdim*(i-1)+1):(xdim*i), (xdim*(i-1)+1):(xdim*i)] = A[:, :]
    Ba[(xdim*(i-1)+1):(xdim*i), (udim*(i-1)+1):(udim*i)] = B[:, :]
  end
  Ap.value = Aa; fix!(Ap)
  Bp.value = Ba; fix!(Bp)
  cstr = [X[(xdim + 1):end] == Ap * X[1:(end - xdim)] + Bp * U,
          X[1:xdim] == x0]

  # enforce control and state limits ##########################################
  if xb != nothing
    cstr = [cstr; reshape(X, xdim, N) >= xb[1]; reshape(X, xdim, N) <= xb[2]]
  end
  if ub != nothing
    cstr = [cstr; reshape(U, udim, N) >= ub[1]; reshape(U, udim, N) <= ub[2]]
  end

  # formulate the objective ###################################################
  # vector-wise quadform is faster than block diagonal quad form
  # vector-wise quadform
  obj = (quadform(reshape(X[1:(end - xdim)] - Xref[1:(end - xdim)], xdim, N), 
                  Q) + 
         quadform(reshape(U - Uref, udim, N), R) + 
         quadform(X[(end - xdim + 1):end] - Xref[(end - xdim + 1):end], P))

  # build the problem and solve ###############################################
  prob = minimize(obj, cstr)
  # redirect printing during execution; the only way to enforce verbose=false
  @suppress solve!(prob, solver, warmstart=true)
  return (X.value, U.value)
end

function scpMPC(Alin, Blin, Q, R, P, x0, N; xb=nothing, ub=nothing,
                  Xref=nothing, Uref=nothing, Xguess=nothing, Uguess=nothing)
  #solver = GurobiSolver()
  #solver = MosekSolver()
  solver = ECOSSolver(maxit=10^3, eps=1e-9, verbose=false)
  # create variables and reference trajectories ###############################
  rho = Variable(); rho_init_val = 1e-1; rho.value = rho_init_val; fix!(rho)
  xdim = size(Q, 1)
  udim = size(R, 1)

  #Xref = Xref == nothing ? [x0; zeros(xdim * N)] : Xref
  Xref = Xref == nothing ? [x0; 10 * zeros(xdim * N)] : Xref
  Uref = Uref == nothing ? zeros(udim * N) : Uref
  X = Variable((N + 1) * xdim)
  U = Variable(N * udim)
  if Xguess != nothing && Uguess != nothing
    rho.value = 1e-5
  end
  Xprev = Variable(xdim * (N + 1))
  Uprev = Variable(udim * N)
  Xprev.value = Xguess == nothing ? Xref : Xguess; fix!(Xprev)
  Uprev.value = Xguess == nothing ? Uref : Uguess; fix!(Uprev)

  # enforce control and state limits ##########################################
  Ap = Variable(xdim * N, xdim * N)
  Bp = Variable(xdim * N, udim * N)
  Aa = blockdiag(map(i -> sparse(Alin(i, 
    Xprev.value[(xdim*(i-1)+1):(xdim*i)])[:, :]), 1:N)...)
  Ba = blockdiag(map(i -> sparse(Blin(i, 
    Xprev.value[(xdim*(i-1)+1):(xdim*i)], 
    Uprev.value[(udim*(i-1)+1):(udim*i)])[:, :]), 1:N)...)
  Ap.value = Aa; fix!(Ap)
  Bp.value = Ba; fix!(Bp)

  cstr = [abs(X[(xdim + 1):end] - Ap * X[1:(end - xdim)] + Bp * U) <= rho]
  cstr = [cstr; X[1:xdim] == x0]
  #cstr = [cstr; abs(X - Xprev) <= 1e0 * repeat(1:(N+1), inner=xdim) * rho]
  #cstr = [cstr; abs(X - Xprev) <= 1000.0]
  if xb != nothing
    cstr = [cstr; reshape(X, xdim, N) >= xb[1]; 
               reshape(X, xdim, N) <= xb[2]]
  end
  if ub != nothing
    cstr = [cstr; reshape(U, udim, N) >= ub[1];
               reshape(U, udim, N) <= ub[2]]
  end

  # formulate the objective ###################################################
  # penalty and trust region
  obj = (quadform(reshape(X[1:(end - xdim)] - Xref[1:(end - xdim)], xdim,
                             N), Q) + quadform(reshape(U - Uref, udim, N), R) +
            quadform(X[(end - xdim + 1):end] - Xref[(end - xdim + 1):end], P))
  obj += 1e2 * rho * sum(abs(X - Xprev))

  # build the problem and solve ###############################################
  prob = minimize(obj, cstr)

  residual = Inf
  # penalty and trust region solution 
  for i in 1:100
    for j in 1:N
      Aa[(xdim*(j-1)+1):(xdim*j), (xdim*(j-1)+1):(xdim*j)] = Alin(j,
          Xprev.value[(xdim*(j-1)+1):(xdim*j)])
      Ba[(xdim*(j-1)+1):(xdim*j), (udim*(j-1)+1):(udim*j)] = Blin(j,
          Xprev.value[(xdim*(j-1)+1):(xdim*j)],
          Uprev.value[(udim*(j-1)+1):(udim*j)])
    end
    Ap.value = Aa
    Bp.value = Ba

    X.value = Xprev.value
    U.value = Uprev.value
    @suppress solve!(prob, solver, warmstart=true)
    if prob.status == :Optimal || prob.status == :Suboptimal
      residual = norm(X.value - Xprev.value) / ((N + 1) * xdim)
      Xprev.value = X.value
      Uprev.value = U.value

      rho.value = max(0.5 * rho.value, 1e-5)
      #println("residual = ", residual, " rho = ", rho.value)
      if rho.value[] <= 1e-5 && residual < 1e-4
        break
      end
    else
      X.value = Xprev.value; U.value = Uprev.value
      rho.value = min(2.0 * rho.value, 1e1)
    end
  end
  #println("RHO VALUE = ", rho.value)
  if rho.value > 2e-5
    @warn "Solution not found , the solution returned is meaningless"
    return (X.value, U.value)
    return (fill(NaN, size(X.value)), fill(NaN, size(U.value)))
  else
    return (X.value, U.value)
  end
end
