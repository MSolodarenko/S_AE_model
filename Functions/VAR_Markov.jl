using LinearAlgebra
using SpecialFunctions
using Statistics

#
function phi(z)
    #
    #  Standard statistical normal distribution cdf
    #
    p = erfc( -z/sqrt(2) )/2
    return p
end
#
function phinv(w)
    #
    #  Standard statistical inverse normal distribution
    #
    z = -sqrt(2)*erfcinv( 2*w )
    return z
end

function fn_var_to_markov(A0,A1,A2,SIGMA,N,random_draws,method)

    #
    ###########################################################################
    ###########################################################################
    ###########################################################################
    #
    # GENERALIZED MARKOV APPROXIMATIONS TO VAR PROCESSES
    #
    # This function converts a VAR to a discretized Markov process,
    # generalizing the approach in Tauchen (1986) by allowing for more general
    # var./cov. structure in the error term.  The required multivariate normal
    # probabilities are calculated using a Monte Carlo-type technique
    # implemented in the function qscmvnv.m, developed by and available on the
    # website of Alan Genz: http://www.math.wsu.edu/faculty/genz/homepage.
    #
    # Original VAR: A0*Z(t) = A1 + A2*Z(t-1) + e(t), e(t) ~ N(0,SIGMA)
    #
    # INPUTS:
    # 1 - A0, A1, A2 are the VAR coefficients, as indicated above, with
    #     A0 assumed non-singular
    # 2 - N is n x 1, where n = # of vars, N[i] = # grid points for ith var.
    # 3 - SIGMA is the arbitrary positive semi-definite error var./cov. matrix
    # 4 - random_draws is the number of random draws used in the required
    #     Monte Carlo-type integration of the multivariate normal
    # 5 - method switch determines the grid selection method
    #       - method = 1 uses a uniformly spaced grid covering a fixed number
    #         of std. dev. of the relevant component variables.  This is the
    #         grid spacing strategy proposed in Tauchen (1986).
    #       - method = 2 selects grid points based on approximately equal
    #         weighting from the UNIVARIATE normal cdf.  This method is adapted
    #         from code written by Jonathan Willis.  (Note that method = 2
    #         requires the use of the MATLAB statistics toolbox.)
    #
    # OUTPUTS:
    # 1 - Pr_mat is the Prod(N) x Prod(N) computed transition probability matrix
    #     for the discretized Markov chain
    # 2 - Pr_mat_key is n x Prod(N) matrix s.t. if Z* is the discretized Markov
    #     approximation to the VAR Z, then Z*(state i) = Pr_mat_key[:,i]
    # 3 - zbar is the n x max(N) matrix s.t. zbar(i,1:N[i]) is the univariate
    #     grid for the ith component of Z*
    ###########################################################################
    ###########################################################################
    ###########################################################################

    n = size(N,1)  #number of variables in VAR

    #compute reduced form parameters & steady-state mean
    A1bar = inv(A0)*A1
    A2bar = inv(A0)*A2
    SIGMAbar = inv(A0)*SIGMA*(inv(A0)')

    sstate_mean = inv(Matrix{Float64}(I, n, n)-A2bar)*A1bar

    m = 3  #number std deviations of the VAR covered by grid

    #iterate to obtain var./cov. structure of the PROCESS (not error term)
    SIGMAprocess = SIGMAbar
    SIGMAprocess_last = SIGMAprocess
    dif = 1
    while dif>0.00000001
        SIGMAprocess = A2bar*SIGMAprocess_last*(A2bar') + SIGMAbar
        dif = maximum(SIGMAprocess-SIGMAprocess_last)
        SIGMAprocess_last = SIGMAprocess
    end
    #SIGMAprocess

    #This block equally spaces grid points bounded by m*(std.deviation of
    #process) on either side of the unconditional mean of the process.  Any
    #more sophisticated spacing of the grid points could be implemented by
    #changing the definition of zbar below.
    zbar = zeros(n,maximum(N))
    grid_stdev = diag(SIGMAprocess).^0.5
    if method==1
        grid_increment = zeros(n,1)
        for i = 1:n
            #grid_increment[i] = 2*m*grid_stdev[i]/(N[i]-1)
            grid_increment[i] = 2*grid_stdev[i]*sqrt(N[i]-1)/(N[i]-1)
            #zbar[i,1] = -m*grid_stdev[i] + sstate_mean[i]
            zbar[i,1] = -grid_stdev[i]*sqrt(N[i]-1) + sstate_mean[i]
            for j = 1:N[i]-1
                zbar[i,j+1] = zbar[i,j] + grid_increment[i]
            end
        end
    #=elseif method==2
        d = zeros(n,max(N))
        b = -4:.005:4
        c = normcdf(b,0,1)
        for i = 1:n
            a = (1/(2*N[i])):(1/N[i]):1
            for j = 1:N[i]
                [d1,d[i,j]] = min((a[j]-c).^2)
            end
            zbar(i,1:N[i]) = grid_stdev[i]*b(d[i,:])+sstate_mean[i]
        end=#
    end

    #compute key matrix & pos matrix
    Pr_mat_key = zeros(length(N),prod(N))
    Pr_mat_key_pos = zeros(length(N),prod(N))

    A = zbar[length(N), 1:N[length(N)]]
    counts = [1 prod(N)/N[length(N)]]
    Pr_mat_key[length(N),:] = repeat(A, Int64(counts[1]), Int64(counts[2]))

    A = collect(1:N[length(N)])
    counts = [1 prod(N)/N[length(N)]]
    Pr_mat_key_pos[length(N),:] = repeat(A, Int64(counts[1]), Int64(counts[2]))

    for i=length(N)-1:-1:1
        A = kron(zbar[i,1:N[i]], ones(1,prod(N[i+1:length(N)]))')
        counts = [1 prod(N)/prod(N[i:length(N)])]
        Pr_mat_key[i,:] = repeat(A,Int64(counts[1]),Int64(counts[2]))

        A = kron(collect(1:N[i]), ones(1,prod(N[i+1:length(N)]))')
        counts = [1 prod(N)/prod(N[i:length(N)])]
        Pr_mat_key_pos[i,:] = repeat(A,Int64(counts[1]),Int64(counts[2]))
    end

    nstate = prod(N)
    Pr_mat_intervals = zeros(n,nstate,2)   #this will store the unadjusted limits of integration for each variable in each state, for input into the Genz code
    if method==1
        for i = 1:nstate  #number of states
            for j = 1:n    #number of variables
                if Int64(Pr_mat_key_pos[j,i])==1
                    Pr_mat_intervals[j,i,1] = -Inf
                    Pr_mat_intervals[j,i,2] = zbar[j,Int64(Pr_mat_key_pos[j,i])] + (grid_increment[j]/2)
                elseif Int64(Pr_mat_key_pos[j,i])==N[j]
                    Pr_mat_intervals[j,i,1] = zbar[j,Int64(Pr_mat_key_pos[j,i])] - (grid_increment[j]/2)
                    Pr_mat_intervals[j,i,2] = Inf
                else
                    Pr_mat_intervals[j,i,1] = zbar[j,Int64(Pr_mat_key_pos[j,i])] - (grid_increment[j]/2)
                    Pr_mat_intervals[j,i,2] = zbar[j,Int64(Pr_mat_key_pos[j,i])] + (grid_increment[j]/2)
                end
            end
        end
    #=elseif method==2
        for i = 1:nstate  #number of states
            for j = 1:n    #number of variables
                if Int64(Pr_mat_key_pos[j,i])==1
                    Pr_mat_intervals[j,i,1] = -Inf
                    Pr_mat_intervals[j,i,2] = zbar[j,Int64(Pr_mat_key_pos[j,i])] + (zbar(j,Int64(Pr_mat_key_pos[j,i])+1)-zbar[j,Int64(Pr_mat_key_pos[j,i])])/2
                elseif Int64(Pr_mat_key_pos[j,i])==N[j]
                    Pr_mat_intervals[j,i,1] = zbar[j,Int64(Pr_mat_key_pos[j,i])] - (zbar[j,Int64(Pr_mat_key_pos[j,i])]-zbar(j,Int64(Pr_mat_key_pos[j,i])-1))/2
                    Pr_mat_intervals[j,i,2] = Inf
                else
                    Pr_mat_intervals[j,i,1] = zbar[j,Int64(Pr_mat_key_pos[j,i])] - (zbar[j,Int64(Pr_mat_key_pos[j,i])]-zbar(j,Int64(Pr_mat_key_pos[j,i])-1))/2
                    Pr_mat_intervals[j,i,2] = zbar[j,Int64(Pr_mat_key_pos[j,i])] + (zbar(j,Int64(Pr_mat_key_pos[j,i])+1)-zbar[j,Int64(Pr_mat_key_pos[j,i])])/2
                end
            end
        end=#
    end

    error_est = zeros(nstate,nstate)
    Pr_mat_intervals_adjusted = zeros(n,nstate,2)
    Pr_mat = zeros(nstate,nstate)
    for i = 1:nstate #rows of Pr_mat
        Pr_mat_intervals_adjusted[:,:,1] = Pr_mat_intervals[:,:,1] - repeat((A1bar + A2bar*Pr_mat_key[:,i]),1,nstate)
        Pr_mat_intervals_adjusted[:,:,2] = Pr_mat_intervals[:,:,2] - repeat((A1bar + A2bar*Pr_mat_key[:,i]),1,nstate)
        for j = 1:nstate   #columns of Pr_mat
            #Pr_mat[i,j] = P(state j|state i)
            #Pr_mat[i,j], error_est[i,j] = qscmvnv(random_draws,SIGMAbar,Pr_mat_intervals_adjusted[:,j,1],Matrix{Float64}(I, n, n),Pr_mat_intervals_adjusted[:,j,2])
            Pr_mat[i,j], error_est[i,j] = qscmvnv(copy(random_draws),copy(SIGMAbar),copy(Pr_mat_intervals_adjusted[:,j,1]),Matrix{Float64}(I, n, n),copy(Pr_mat_intervals_adjusted[:,j,2]))
            if isnan(Pr_mat[i,j])
                println([i,j])
            end
        end
    end

    #rounding error adjustment
    round_sum = sum(Pr_mat;dims=2)
    for i = 1:size(Pr_mat,2)
        Pr_mat[i,:] = Pr_mat[i,:]/round_sum[i]
    end

    return Pr_mat,Pr_mat_key,zbar
end

function qscmvnv( m, r, a, cn, b )
    #=m = random_draws
    r = SIGMAbar
    a = Pr_mat_intervals_adjusted[:,1,1]
    cn = Matrix{Float64}(I, n, n)
    b = Pr_mat_intervals_adjusted[:,1,2]=#
    #
    #  [ P E ] = QSCMVNV( M, R, A, CN, B )
    #    uses a randomized quasi-random rule with m points to estimate an
    #    MVN probability for positive semi-definite covariance matrix r,
    #    with constraints a < cn*x < b. If r is nxn and cn is kxn, then
    #    a and b must be column k-vectors.
    #   Probability p is output with error estimate e.
    #    Example usage:
    #     >> r = [ 4 3 2 1; 3 5 -1 1; 2 -1 4 2; 1 1 2 5 ];
    #     >> a = [ -inf 1 -5 ]'; b = [ 3 inf 4 ]';
    #     >> cn = [ 1 2 3 -2; 2 4 1 2; -2 3 4 1 ];
    #     >> [ p e ] = qscmvnv( 5000, r, a, cn, b ); disp([ p e ])
    #
    #  This function uses an algorithm given in the paper by Alan Genz:
    #   "Numerical Computation of Multivariate Normal Probabilities", in
    #     J. of Computational and Graphical Stat., 1(1992), 141-149.
    #  The primary references for the numerical integration are
    #   "On a Number-Theoretical Integration Method"
    #     H. Niederreiter, Aequationes Mathematicae, 8(1972), 304-11, and
    #   "Randomization of Number Theoretic Methods for Multiple Integration"
    #     R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13(1976), 904-14.
    #
    #   Alan Genz is the author of this function and following Matlab functions.
    #          Alan Genz, WSU Math, PO Box 643113, Pullman, WA 99164-3113
    #          Email : AlanGenz@wsu.edu
    #
    # Initialization
    #
    as, ch, bs, clg, n = chlsrt( copy(r), copy(a), copy(cn), copy(b) )
    ci = phi(as[1])
    dci = phi(bs[1]) - ci
    p = 0
    e = 0
    ns = 8
    nv = Int64(trunc( max(m/(2*ns),1.0) ))
    q = 2.0 .^( collect(1:n-1)'./n)
    #
    # Randomization loop for ns samples
    #
    xx = Array{Float64}(undef, n-1,nv)
    for i = 1:ns
      # periodizing transformation
      xx[:,1:nv] = abs.( 2*mod.( q*collect(1:nv)' + rand(n-1,1)*ones(1,nv), 1 ) .- 1 )
      #vp =   mvndnv( n, as, ch, bs, clg, ci, dci,   xx, nv )
      vp =   mvndnv( copy(n), copy(as), copy(ch), copy(bs), copy(clg), copy(ci), copy(dci),   copy(xx), copy(nv) )
      #vp = ( mvndnv( n, as, ch, bs, clg, ci, dci, 1.0.-xx, nv ) + vp )/2 # symmetrize
      vp = ( mvndnv( copy(n), copy(as), copy(ch), copy(bs), copy(clg), copy(ci), copy(dci), 1.0.-xx, copy(nv) ) + vp )/2 # symmetrize
      d = ( mean(vp) - p )/i
      p = p + d
      if abs(d) > 0
        e = abs(d)*sqrt( 1 + ( e/d )^2*( i - 2 )/i )
      else
        if i > 1
            e = e*sqrt( ( i - 2 )/i )
        end
      end
    end
    #
    e = 3*e # error estimate is 3 x standard error with ns samples.

    return p, e

end
#
nanmax(x,y) = isnan(x) ? y : (isnan(y) ? x : max(x,y))
nanmin(x,y) = isnan(x) ? y : (isnan(y) ? x : min(x,y))
function mvndnv( n, a, ch, b, clg, ci, dci, x, nv )
    #
    #  Transformed integrand for computation of MVN probabilities.
    #
    y = zeros(n-1,nv)
    on = ones(1,nv)
    c = ci*on
    dc = dci*on
    p = dc
    li = 2
    lf = 1
    for i = 2:n
       y[i-1,:] = phinv.( c + x[i-1,:]'.*dc )
       lf = lf + clg[i]
       if lf < li
           c = 0
           dc = 1
       else
           s = ch[li:lf,1:i-1]*y[1:i-1,:]
           #ai =           max.( maximum( a[li:lf]*on - s, dims = 1 ), -9 )
           ai = nanmax.(maximum( a[li:lf]*on - s, dims = 1 ), -9)
           #bi = max.( ai, min.( minimum( b[li:lf]*on - s, dims = 1 ),  9 ) )
           bi = max.(ai, nanmin.(minimum( b[li:lf]*on - s, dims = 1 ),  9) )
           c = phi.(ai)
           dc = phi.(bi) - c
           p = p.*dc
       end
       li = li + clg[i]
    end
    return p
end
#
function chlsrt( r, a, cn, b )

    #
    #  Computes permuted lower Cholesky factor ch for covariance r which
    #   may be singular, combined with contraints a < cn*x < b, to
    #   form revised lower triangular constraint set ap < ch*x < bp;
    #   clg contains information about structure of ch: clg(1) rows for
    #   ch with 1 nonzero, ..., clg(np) rows with np nonzeros.
    #
    ep = 1e-10 # singularity tolerance;
    #
    n,n = size(r)
    m,n = size(cn)
    ch = cn
    np = 0
    ap = a
    bp = b
    y = zeros(n,1)
    sqtp = sqrt(2*pi)
    c = r
    d = sqrt.(max.(diag(c),0))
    for i = 1:n
        di = d[i]
        if di > 0
            c[:,i] = c[:,i]/di
            c[i,:] = c[i,:]/di
            ch[:,i] = ch[:,i]*di
        end
    end
    #
    # determine (with pivoting) Cholesky factor for r
    #  and form revised constraint matrix ch
    #
    clg = Array{Int64,1}(undef,n)
    for i = 1:n
        #np = np + 1
      clg[i] = 0
      epi = ep * i^2
      j = i
      for l = i+1:n
          if c[l,l] > c[j,j]
              j = l
          end
      end
      if j > i
        t = c[i,i]
        c[i,i] = c[j,j]
        c[j,j] = t

        t = c[i,1:i-1]
        c[i,1:i-1] = c[j,1:i-1]
        c[j,1:i-1] = t

        t = c[i+1:j-1,i]
        c[i+1:j-1,i] = c[j,i+1:j-1]'
        c[j,i+1:j-1] = t'

        t = c[j+1:n,i]
        c[j+1:n,i] = c[j+1:n,j]
        c[j+1:n,j] = t

        t = ch[:,i]
        ch[:,i] = ch[:,j]
        ch[:,j] = t
      end
      if c[i,i] < epi
          break
      end
      cvd = sqrt( c[i,i] )
      c[i,i] = cvd
      for l = i+1:n
        c[l,i] = c[l,i]/cvd
        c[l,i+1:l] = c[l,i+1:l] - c[l,i]*c[i+1:l,i]'
      end
      ch[:,i] = ch[:,i:n]*c[i:n,i]
      np = np + 1
    end
    #
    # use right reflectors to reduce ch to lower triangular
    #
    for i = 1:min(np-1,m)
      epi = ep * i * i
      vm = 1
      lm = i
      #
      # permute rows so that smallest variance variables are first.
      #
      for l = i:m
        v = ch[l,1:np]
        s = v[1:i-1]' * y[1:i-1]
        ss = max( sqrt( sum( v[i:np]'.^2 ) ), epi )
        al = ( ap[l] - s )/ss
        bl = ( bp[l] - s )/ss
        dna = 0
        dsa = 0
        dnb = 0
        dsb = 1
        if al > -9
            dna = exp(-al*al/2)/sqtp
            dsa = phi(al)
        end
        if bl <  9
            dnb = exp(-bl*bl/2)/sqtp
            dsb = phi(bl)
        end
        if dsb - dsa > epi
          if      al <= -9
              mn =      -dnb
              vr =         -bl*dnb
          elseif  bl >=  9
              mn = dna
              vr = al*dna
          else
              mn = dna - dnb
              vr = al*dna - bl*dnb
          end
          mn = mn/( dsb - dsa )
          vr = 1 + vr/( dsb - dsa ) - mn^2
        else
          if     al <= -9
              mn = bl
          elseif bl >=  9
              mn = al
          else
              mn = ( al + bl )/2
          end
          vr = 0
        end
        if vr <= vm
            lm = l
            vm = vr
            y[i] = mn
        end
      end
      v = ch[lm,1:np]
      if lm > i
        ch[lm,1:np] = ch[i,1:np]
        ch[i,1:np] = v

        tl = ap[i]
        ap[i] = ap[lm]
        ap[lm] = tl

        tl = bp[i]
        bp[i] = bp[lm]
        bp[lm] = tl
      end
      ch[i,i+1:np] .= 0
      ss = sum( v[i+1:np].^2 )
      if    ss > epi
        ss = sqrt( ss + v[i]^2 )
        if v[i] < 0
            ss = -ss
        end
        ch[i,i] = -ss
        v[i] = v[i] + ss
        vt = collect(v[i:np]/( ss*v[i] ))
        ch[i+1:m,i:np] = ch[i+1:m,i:np] - ch[i+1:m,i:np] * vt * v[i:np]'
      end
    end
    #
    # scale and sort constraints
    #
    clm = Array{Int64}(undef,m)
    for i = 1:m
      v = ch[i,1:np]
      clm[i] = min(i,np)
      jm = 1
      for j = 1:clm[i]
          if abs(v[j]) > ep*j*j
              jm = j
          end
      end
      if jm < np
          v[jm+1:np] .= 0
      end
      clg[jm] = clg[jm] + 1
      at = ap[i]
      bt = bp[i]
      j = i
      for l = i-1:-1:1
        if jm >= clm[l]
            break
        end
        ch[l+1,1:np] = ch[l,1:np]
        j = l
        ap[l+1] = ap[l]
        bp[l+1] = bp[l]
        clm[l+1] = clm[l]
      end
      clm[j] = jm
      vjm = v[jm]
      ch[j,1:np] = v/vjm
      ap[j] = at/vjm
      bp[j] = bt/vjm
      if vjm < 0
          tl = ap[j]
          ap[j] = bp[j]
          bp[j] = tl
      end
    end
    j = 0
    for i = 1:np
        if clg[i] > 0
            j = i
        end
    end
    np = j
    #
    # combine constraints for first variable
    #
    if clg[1] > 1
      ap[1] = maximum( ap[1:clg[1]] )
      bp[1] = max( ap[1], minimum( bp[1:clg[1]] ) )
      ap[2:m-clg[1]+1] = ap[clg[1]+1:m]
      bp[2:m-clg[1]+1] = bp[clg[1]+1:m]
      ch[2:m-clg[1]+1,:] = ch[clg[1]+1:m,:]
      clg[1] = 1
    end

    return ap, ch, bp, clg, np
end

function fn_var_to_markov(A1,A2,SIGMA,N,random_draws,method)
    return fn_var_to_markov([1.0 0.0; 0.0 1.0],A1,A2,SIGMA,N,random_draws,method)
end

function fn_var_to_markov(A2,SIGMA,N,random_draws,method)
    return fn_var_to_markov([0.0; 0.0],A2,SIGMA,N,random_draws,method)
end

function fn_var_to_markov(A2,SIGMA,N)
    return fn_var_to_markov(A2,SIGMA,N,1000,1)
end
#=
number_u_m_nodes = 3#5#7
number_u_w_nodes = 3#5#7

sigma_eps_m = 1.145 #Variance innovation managerial skill
sigma_eps_w = 0.073 #Variance innovation working skill

rho_m = 0.788 #Autocorrelation managerial skill
rho_w = 0.96  #Autocorrelation of working-ability shocks

rho_eps_m_w = 0.335

A0 = [1.0 0.0; 0.0 1.0]
A1 = [0.0; 0.0]
A2 = [rho_w 0.0; 0.0 rho_m]
SIGMA = [sigma_eps_w rho_eps_m_w*sqrt(sigma_eps_w*sigma_eps_m); rho_eps_m_w*sqrt(sigma_eps_w*sigma_eps_m) sigma_eps_m]
N = [number_u_w_nodes; number_u_m_nodes]
random_draws = 1000
method = 1
P,P_key,nodes = fn_var_to_markov(A0,A1,A2,SIGMA,N,random_draws,method)
display(P)
display(P_key)
display(nodes)
=#

import QuantEcon: tauchen,std_norm_cdf,MarkovChain
function tauchen(N::Integer, ρ::T1, σ::T2, μ=zero(promote_type(T1, T2)), n_std::Float64=3) where {T1 <: Real, T2 <: Real}
    # Get discretized space
    a_bar = n_std * sqrt(σ^2 / (1 - ρ^2))
    y = range(-a_bar, stop=a_bar, length=N)
    d = y[2] - y[1]

    # Get transition probabilities
    Π = zeros(promote_type(T1, T2), N, N)
    for row = 1:N
        # Do end points first
        Π[row, 1] = std_norm_cdf((y[1] - ρ*y[row] + d/2) / σ)
        Π[row, N] = 1 - std_norm_cdf((y[N] - ρ*y[row] - d/2) / σ)

        # fill in the middle columns
        for col = 2:N-1
            Π[row, col] = (std_norm_cdf((y[col] - ρ*y[row] + d/2) / σ) -
                           std_norm_cdf((y[col] - ρ*y[row] - d/2) / σ))
        end
    end

    # NOTE: I need to shift this vector after finding probabilities
    #       because when finding the probabilities I use a function
    #       std_norm_cdf that assumes its input argument is distributed
    #       N(0, 1). After adding the mean E[y] is no longer 0, so
    #       I would be passing elements with the wrong distribution.
    #
    #       It is ok to do after the fact because adding this constant to each
    #       term effectively shifts the entire distribution. Because the
    #       normal distribution is symmetric and we just care about relative
    #       distances between points, the probabilities will be the same.
    #
    #       I could have shifted it before, but then I would need to evaluate
    #       the cdf with a function that allows the distribution of input
    #       arguments to be [μ/(1 - ρ), 1] instead of [0, 1]

    yy = y .+ μ / (1 - ρ) # center process around its mean (wbar / (1 - rho)) in new variable

    # renormalize. In some test cases the rows sum to something that is 2e-15
    # away from 1.0, which caused problems in the MarkovChain constructor
    Π = Π./sum(Π, dims = 2)

    MarkovChain(Π, yy)
end
