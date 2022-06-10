using SchumakerSpline
using Statistics

include("print_sameline.jl")
print_sameline("Loading functions for calculating profit and income")
include("AllubErosa.jl")
include("profit.jl")

function transitional_dynamics(lambda_s, ss_star, ss_starstar, global_params, global_approx_params, model_params, file_name, GUESS_RS, maxiters)
    # Initialisation of model's internal parameters
    T = length(lambda_s)

    # approximation object for skills nodes
    #               1               2       3       4       5       6           7
    #               z_m_nodes, z_w_nodes, P_zeta, P_u, stat_P_u, P_alpha, number_u_nodes
    #approx_object = build_skill_nodes(global_approx_params, model_params)
    approx_object = ss_star[4]

    distr_tol = global_params[3]
    val_tol = global_params[4]

    number_a_nodes = global_approx_params[1]
    number_zeta_nodes = global_approx_params[4]
    number_alpha_m_nodes = global_approx_params[5]
    number_alpha_w_nodes = global_approx_params[6]

    beta = model_params[2]
    p_alpha = model_params[13]

    z_m_nodes = approx_object[1]
    z_w_nodes = approx_object[2]
    P_zeta = approx_object[3]
    P_u = approx_object[4]
    stat_P_u = approx_object[5]
    P_alpha = approx_object[6]
    number_u_nodes = approx_object[7]

    delta = model_params[3]
    gamma = model_params[4]
    eta = model_params[5]
    theta = model_params[6]
    c_e = model_params[7]

    z_m_nodes = approx_object[1]
    z_w_nodes = approx_object[2]
    number_u_nodes = approx_object[7]

    ############################ use a_nodes of the latter economy (ss_starstar)
    a_min = 0.01
    #!!!!!!!!!!!!!!! No idea why it doesn't work with max between a_max's of different ss'
    a_max = ss_starstar[1][1]#max(ss_star[1][1], ss_starstar[1][1])#
    #!!!!!!!!!!!!!!!
    a_nodes = exp.(collect(range(log(a_min+1); stop=log(a_max+1), length=number_a_nodes))).-1

    policies_small_grid_T = zeros(number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    values_T = zeros(number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)

    number_asset_grid = number_a_nodes*10
    asset_grid = exp.(collect(range(log(a_min+1); stop=log(a_max+1), length=number_asset_grid))).-1
    distrs_1 = zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
        policyT = Schumaker(ss_starstar[1][3],ss_starstar[1][4][:,u_i,zeta_i,alpha_m_i,alpha_w_i]; extrapolation = (Linear,Linear))
        policies_small_grid_T[:,u_i,zeta_i,alpha_m_i,alpha_w_i] .= evaluate.(policyT, a_nodes)

        valueT = Schumaker( exp.(collect(range(log(a_min+1); stop=log(ss_starstar[1][1]+1), length=number_a_nodes))).-1, ss_starstar[1][33][:,u_i,zeta_i,alpha_m_i,alpha_w_i]; extrapolation = (Linear,Linear))
        values_T[:,u_i,zeta_i,alpha_m_i,alpha_w_i] .= evaluate.(valueT, a_nodes)

        distr1 = Schumaker(ss_star[1][3], cumsum(ss_star[1][5][:,u_i,zeta_i,alpha_m_i,alpha_w_i]); extrapolation=(Linear,Linear))
        distrs_1[:,u_i,zeta_i,alpha_m_i,alpha_w_i] .= diff([0.0;evaluate.(distr1, asset_grid)])
        distrs_1[:,u_i,zeta_i,alpha_m_i,alpha_w_i] .*= sum(ss_star[1][5][:,u_i,zeta_i,alpha_m_i,alpha_w_i])/sum(distrs_1[:,u_i,zeta_i,alpha_m_i,alpha_w_i])
    end

    ###########################
    function plot_iter_results(r_s,w_s, new_r_s,new_w_s, best_r_s,best_w_s, len_r_s,len_w_s, best_len_r_s,best_len_w_s)
        p5 = Plots.plot(r_s.-best_r_s, legend=false)
        hline!(p5,[0.0])
        p6 = Plots.plot(w_s.-best_w_s, legend=false)
        hline!(p6,[0.0])

        p1 = Plots.plot([r_s, ones(T).*ss_star[2],ones(T).*ss_starstar[2], new_r_s, best_r_s], label=["Current" "" "" "New" "Best"], legend=false)
        p2 = Plots.plot([w_s, ones(T).*ss_star[3],ones(T).*ss_starstar[3], new_w_s, best_w_s], label=["Current" "" "" "New" "Best"], legend=:bottomright)
        #display(Plots.plot(p1,p2, layout=(1,2)))
        p3 = Plots.plot([len_r_s, best_len_r_s], label=["Current" "Best"], legend=false)
        hline!(p3,[0.0])
        p4 = Plots.plot([len_w_s, best_len_w_s], label=["Current" "Best"], legend=false)
        hline!(p4,[0.0])
        #display(Plots.plot(p1,p2,p3,p4, layout=(2,2)))
        #display(Plots.plot(p1,p2,p3,p4,p5,p6, layout=(3,2)))
        display(Plots.plot(p5,p6,p1,p2,p3,p4, layout=(3,2)))
    end

    function calculate_policies(r_s, w_s)
        policies_small_grid = zeros(T, number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        policies_small_grid[T,:,:,:,:,:] .= policies_small_grid_T

        policy_function = Array{Any}(undef,T,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)

        values = zeros(T, number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        values[T,:,:,:,:,:] .= values_T

        policies = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        a1_indices = Array{Int64}(undef,T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        lottery_prob_1 = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        lottery_prob_2 = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)

        Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            policy_function[T,u_i,zeta_i,alpha_m_i,alpha_w_i] = Schumaker(a_nodes,policies_small_grid[T,:,u_i,zeta_i,alpha_m_i,alpha_w_i]; extrapolation = (Linear,Linear))
            for a_i in 1:number_asset_grid
                try
                    a1 = evaluate(policy_function[T,u_i,zeta_i,alpha_m_i,alpha_w_i], asset_grid[a_i])
                    policies[T, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = a1
                    j_1 = sum(a1 .>= asset_grid)
                    j = j_1 + 1

                    a1_indices[T, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = j_1
                    if j <= number_asset_grid && j_1 >= 1
                        lottery_prob_1[T, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (asset_grid[j] - a1             )/(asset_grid[j]-asset_grid[j_1])
                        lottery_prob_2[T, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (a1            - asset_grid[j_1])/(asset_grid[j]-asset_grid[j_1])
                    elseif j_1 == number_asset_grid
                        lottery_prob_1[T, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 1.0
                        lottery_prob_2[T, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 0.0
                    end
                catch e
                    println_sameline("Trouble is in loop for increasing grid for policies")
                    throw(e)
                end
            end
        end

        for t = T-1:-1:1

            income, earnings = compute_income_and_earnings(a_nodes,number_a_nodes,r_s[t],w_s[t], number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, lambda_s[t], delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)

            # compute future payoffs if no change in fixed effects (no death)
            value_tran = zeros(number_a_nodes,number_u_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
            for u_i in 1:number_u_nodes
                for alpha_m_i in 1:number_alpha_m_nodes
                    for alpha_w_i in 1:number_alpha_w_nodes
                        for zeta_prime_i in 1:number_zeta_nodes
                            value_tran[:,u_i,alpha_m_i,alpha_w_i] .+= P_zeta[zeta_prime_i].*values[t+1,:,u_i,zeta_prime_i,alpha_m_i,alpha_w_i]
                        end
                    end
                end
            end

            value_tran_rhs = Array{Any}(undef,number_u_nodes,number_alpha_m_nodes,number_alpha_w_nodes)

            Threads.@threads for (u_prime_i,(alpha_m_i,alpha_w_i)) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes)))
                value_tran_rhs[u_prime_i,alpha_m_i,alpha_w_i] = Schumaker(a_nodes,value_tran[:,u_prime_i,alpha_m_i,alpha_w_i]; extrapolation = (Linear,Linear))
            end

            expectation_value_death = zeros(number_a_nodes)
            for aprime_i in 1:number_a_nodes
                for alpha_m_prime_i in 1:number_alpha_m_nodes
                    for alpha_w_prime_i in 1:number_alpha_w_nodes
                        expectation_value_death[aprime_i] += beta*p_alpha*P_alpha[alpha_m_prime_i,alpha_w_prime_i] * sum(value_tran[aprime_i,:,alpha_m_prime_i,alpha_w_prime_i].*stat_P_u)
                    end
                end
            end

            expectation_value_death_rhs = Schumaker(a_nodes,expectation_value_death; extrapolation = (Linear,Linear))

            # loop for finding new_value and new_aprime_indices
            Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))

                values[t,:,u_i,zeta_i,alpha_m_i,alpha_w_i], policies_small_grid[t,:,u_i,zeta_i,alpha_m_i,alpha_w_i] = new_val_and_a1(values[t+1,:,u_i,zeta_i,alpha_m_i,alpha_w_i],policies_small_grid[t+1,:,u_i,zeta_i,alpha_m_i,alpha_w_i], income[:,u_i,zeta_i,alpha_m_i,alpha_w_i], value_tran_rhs[:,alpha_m_i,alpha_w_i],expectation_value_death_rhs, u_i,alpha_m_i,alpha_w_i ,a_min,a_max,a_nodes, 0.0,0.0, number_a_nodes, beta, p_alpha, P_u, number_u_nodes, CRRA)

                policy_function[t,u_i,zeta_i,alpha_m_i,alpha_w_i] = Schumaker(a_nodes,policies_small_grid[t,:,u_i,zeta_i,alpha_m_i,alpha_w_i]; extrapolation = (Linear,Linear))

                for a_i in 1:number_asset_grid
                    try
                        a1 = evaluate(policy_function[t,u_i,zeta_i,alpha_m_i,alpha_w_i], asset_grid[a_i])
                        policies[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = a1
                        j_1 = sum(a1 .>= asset_grid)
                        j = j_1 + 1

                        a1_indices[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = j_1
                        if j <= number_asset_grid && j_1 >= 1
                            lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (asset_grid[j] - a1             )/(asset_grid[j]-asset_grid[j_1])
                            lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (a1            - asset_grid[j_1])/(asset_grid[j]-asset_grid[j_1])
                        elseif j_1 == number_asset_grid
                            lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 1.0
                            lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 0.0
                        end
                    catch e
                        println_sameline("Trouble is in loop for increasing grid for policies")
                        throw(e)
                    end
                end
            end

            print_sameline("#$(rw_iters) - Policy iteration: $(t)/$(T)")
        end

        #=
        Threads.@threads for (t,(u_i,(zeta_i,(alpha_m_i,alpha_w_i)))) in collect(Iterators.product(1:T,Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes)))))
            try
                policy_function = Schumaker(a_nodes,policies_small_grid[t,:,u_i,zeta_i,alpha_m_i,alpha_w_i]; extrapolation = (Linear,Linear))
                for a_i in 1:number_asset_grid
                    a1 = evaluate(policy_function, asset_grid[a_i])
                    policies[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = a1
                    j_1 = sum(a1 .>= asset_grid)
                    j = j_1 + 1

                    a1_indices[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = j_1
                    if j <= number_asset_grid && j_1 >= 1
                        lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (asset_grid[j] - a1             )/(asset_grid[j]-asset_grid[j_1])
                        lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (a1            - asset_grid[j_1])/(asset_grid[j]-asset_grid[j_1])
                    elseif j_1 == number_asset_grid
                        lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 1.0
                        lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 0.0
                    end
                end
            catch e
                println_sameline("Trouble is in loop for increasing grid for policies")
                throw(e)
            end
        end
        =#

        return policies, a1_indices, lottery_prob_1, lottery_prob_2
    end

    function calculate_distrs(policies, a1_indices, lottery_prob_1, lottery_prob_2, r_s, w_s)
        distrs = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        distrs[1,:,:,:,:,:] = distrs_1

        for t=1:T-1

            distr_a1_z0 = copy(distrs[t,:,:,:,:,:])
            distr_a1_z0 .= 0.0
            Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
                try
                    non_zero_asset_grid_iters = findall(!iszero,distrs[t,:,u_i,zeta_i,alpha_m_i,alpha_w_i])
                    for a_i in non_zero_asset_grid_iters#1:number_asset_grid#
                        j_1 = a1_indices[t,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        distr_a1_z0[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += distrs[t,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[t,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        if j_1 != number_asset_grid
                            distr_a1_z0[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += distrs[t,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[t,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        end
                    end

                catch e
                    println_sameline(("Trouble is in loop 2",e))
                    throw(e)
                end
            end

            new_distr = copy(distrs[t,:,:,:,:,:])
            new_distr .= 0.0
            Threads.@threads for (alpha_m_i,alpha_w_i) in collect(Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))
                try
                    # zeta is the transitory shock,
                    # so add over all levels of zeta
                    # and then draw new u_prime and new zeta_prime
                    temp_distr_sum_zeta = sum(distr_a1_z0[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
                    #println_sameline(temp_distr_sum_zeta[:,1])
                    #throw(error)
                    temp_distr_sum_zeta = temp_distr_sum_zeta*P_u
                    #println_sameline(temp_distr_sum_zeta[:,1])
                    #throw(error)
                    for zeta_prime_i in 1:number_zeta_nodes
                        new_distr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_distr_sum_zeta
                    end
                catch e
                    println_sameline("Trouble is in loop 3")
                    throw(e)
                end
            end

            new_distr2 = copy(distrs[t,:,:,:,:,:])
            # second calculate transition
            # for people that change alpha_m and alpha_w
            #   and therefore draw new alpha_m_prime, alpha_w_prime
            #       and zeta_prime, u_prime
            #           from stationary distributions for this shocks

            # calculate sum of capital of all people who change skills
            distr_marginal_assets = sum(sum(sum(sum(distr_a1_z0,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
            new_distr2 .= 0.0

            Threads.@threads for (u_prime_i,(zeta_prime_i,(alpha_m_prime_i,alpha_w_prime_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
                try
                    new_distr2[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] = (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*distr_marginal_assets

                catch e
                    println_sameline("Trouble is in loop 4")
                    throw(e)
                end
            end

            new_distr .+= new_distr2

            distrs[t+1,:,:,:,:,:] .= new_distr

            print_sameline("#$(rw_iters) - Simulation of the history: t$(t)/$(T) - $(round(sum(distrs[t+1,:,:,:,:,:]);digits=6))")
        end

        return distrs
    end

    function partial_equilibrium_transition(r_s, w_s)
        # calculate path for policies from T to 1 (a1_indices and lottery_prob_1_2)
        policies, a1_indices, lottery_prob_1, lottery_prob_2 = calculate_policies(Rs, Ws)

        # calculate path for distributions from 1 to T
        distrs = calculate_distrs(policies, a1_indices, lottery_prob_1, lottery_prob_2, r_s, w_s)

        capital_d_s = Array{Float64}(undef,T,number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        labour_d_s = Array{Float64}(undef,T,number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        labour_s_s = Array{Float64}(undef,T,number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        Threads.@threads for (t,(a_i,(u_i,(zeta_i,(alpha_m_i,alpha_w_i))))) in collect(Iterators.product(1:T,Iterators.product(1:number_asset_grid,Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))))
            t1, t2, t3, t4, capital_d_s[t,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i], t6, t7, labour_d_s[t,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i], labour_s_s[t,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i], t8, t9, t10, t11 = compute_income_profile(asset_grid[a_i],z_m_nodes[u_i,alpha_m_i,zeta_i],z_w_nodes[u_i,alpha_m_i,alpha_w_i],r_s[t],w_s[t], lambda_s[t], delta, gamma, eta, theta, c_e)
        end

        # calculate K_demand
        K_d = sum(distrs.*capital_d_s, dims=2:6)[:,1,1,1,1,1]
        # calculate K_supply
        K_s = sum( sum(distrs; dims=3:6)[:,:,1,1,1,1] .* reshape(asset_grid,1,number_asset_grid), dims=2)[:,1]

        # calculate L_demand
        L_d = sum(distrs.*labour_d_s, dims=2:6)[:,1,1,1,1,1]
        # calculate L_supply
        L_s = sum(distrs.*labour_s_s, dims=2:6)[:,1,1,1,1,1]

        # calculate len_Rs and len_Ws (other lens)
        len_Rs = (K_d.-K_s)./(K_d.+K_s)
        len_Ws = (L_d.-L_s)./(L_d.+L_s)

        return [len_Rs, len_Ws, K_d,K_s,L_d,L_s, distrs,policies]
    end

    # initialisation of procedure
    println_sameline("Preparing initial variables")

    rw_iters = 1
    rw_maxiters = maxiters
    gen_tol_x = global_params[1]
    gen_tol_f = global_params[2]
    relax_r = 0.05
    relax_w = 0.1

    relax_R = relax_r
    relax_W = relax_w

    r_min = -model_params[3]
    r_max = 1/model_params[2]-1.0

    w_min = W_MIN
    w_max = W_MAX

    # guess Rs and Ws
    Rs = (log.(collect(range(exp(-2.0); stop=exp(10.0), length=T))) .- (-2.0)).*(ss_starstar[2] - ss_star[2])./(10 - (-2)) .+ ss_star[2]
    Ws = (log.(collect(range(exp(-2.0); stop=exp(10.0), length=T))) .- (-2.0)).*(ss_starstar[3] - ss_star[3])./(10 - (-2)) .+ ss_star[3]
    if GUESS_RS
        Rs = [ss_star[2]; (1.0./(collect(range(0.5; stop=2.0, length=T-1))) .- 1.0/0.5).*(ss_starstar[2] - 2*ss_starstar[2])./(1.0/(2.0) - 1.0/0.5).+2*ss_starstar[2] ]
    end

    println_sameline("Initiate the search for general equilibrium factor prices")
    # calculate len_Rs and len_Ws by running partial_equilibrium_transition
    #       1       2       3   4   5   6    7      8
    # res = len_Rs, len_Ws, K_d,K_s,L_d,L_s, distrs,policies
    res = partial_equilibrium_transition(Rs, Ws)

    len_Rs = copy(res[1])
    len_Ws = copy(res[2])
    total_len = sum(abs,len_Rs)+sum(abs,len_Ws)
    abs_len_r = maximum(abs,len_Rs)
    abs_len_w = maximum(abs,len_Ws)

    new_Rs = Rs.*(1.0.+len_Rs)./(1.0.-len_Rs) .+ 2.0.*(1.0-r_min).*len_Rs./(1.0.-len_Rs)
    new_Rs = min.(max.(r_min, new_Rs), r_max).*relax_r .+ (1-relax_r).*Rs
    #new_Rs[end] = ss_starstar[2]
    new_Ws = Ws.*(1.0.+len_Ws)./(1.0.-len_Ws)
    new_Ws = min.(max.(w_min, new_Ws), w_max).*relax_w .+ (1-relax_w).*Ws
    #new_Ws[end] = ss_starstar[3]

    best_Rs = copy(Rs)
    best_Ws = copy(Ws)
    best_res = copy(res)
    best_len_Rs = copy(len_Rs)
    best_len_Ws = copy(len_Ws)
    best_total_len = total_len
    best_abs_len_r = abs_len_r
    best_abs_len_w = abs_len_w


    println_sameline("#0 - total_len:$(round(total_len;digits=6)) - total_len_Rs:$(round(sum(abs,len_Rs);digits=6)) - total_len_Ws:$(round(sum(abs,len_Ws);digits=6)) - abs_len_r:$(round(abs_len_r;digits=6)) - abs_len_w:$(round(abs_len_w;digits=6))")
    plot_iter_results(Rs,Ws, new_Rs,new_Ws, best_Rs,best_Ws, len_Rs,len_Ws, best_len_Rs,best_len_Ws)

    println_sameline("#$(rw_iters) - New best Rs and Ws")

    #=
    past_Rs = []
    past_Ws = []
    past_len_Rs = []
    past_len_Ws = []
    push!(past_Rs,Rs)
    push!(past_Ws,Ws)
    push!(past_len_Rs, len_Rs)
    push!(past_len_Ws, len_Ws)
    =#
    old_old_res = copy(res)
    old_old_Rs = copy(Rs)
    old_old_len_Rs = copy(len_Rs)
    old_old_Ws = copy(Ws)
    old_old_len_Ws = copy(len_Ws)
    old_old_total_len = total_len
    old_old_abs_len_r = abs_len_r
    old_old_abs_len_w = abs_len_w
    old_old_new_Rs = new_Rs
    old_old_new_Ws = new_Ws

    old_res = copy(res)
    old_Rs = copy(Rs)
    old_len_Rs = copy(len_Rs)
    old_Ws = copy(Ws)
    old_len_Ws = copy(len_Ws)
    old_total_len = total_len
    old_abs_len_r = abs_len_r
    old_abs_len_w = abs_len_w
    old_new_Rs = new_Rs
    old_new_Ws = new_Ws

    while rw_iters <= rw_maxiters && total_len>gen_tol_f*2*T
        old_old_res = copy(old_res)
        old_old_Rs = copy(old_Rs)
        old_old_len_Rs = copy(old_len_Rs)
        old_old_Ws = copy(old_Ws)
        old_old_len_Ws = copy(old_len_Ws)
        old_old_total_len = old_total_len
        old_old_abs_len_r = old_abs_len_r
        old_old_abs_len_w = old_abs_len_w
        old_old_new_Rs = copy(old_new_Rs)
        old_old_new_Ws = copy(old_new_Ws)

        old_res = copy(res)
        old_Rs = copy(Rs)
        old_len_Rs = copy(len_Rs)
        old_Ws = copy(Ws)
        old_len_Ws = copy(len_Ws)
        old_total_len = total_len
        old_abs_len_r = abs_len_r
        old_abs_len_w = abs_len_w
        old_new_Rs = copy(new_Rs)
        old_new_Ws = copy(new_Ws)

        Rs = copy(new_Rs)
        Ws = copy(new_Ws)
        try
            #=
            if !maximum([Rs.-best_Rs Ws.-best_Ws] .> gen_tol_x)
                throw(error("same as best_Rs and best_Ws"))
            end
            if !maximum([Rs.-old_Rs Ws.-old_Ws] .> gen_tol_x)
                throw(error("same as old_Rs and old_Ws"))
            end
            if !maximum([Rs.-old_old_Rs Ws.-old_old_Ws] .> gen_tol_x)
                throw(error("same as old_old_Rs and old_old_Ws"))
            end
            =#
            # calculate len_Rs and len_Ws by running partial_equilibrium_transition
            #       1       2       3   4   5   6    7      8
            # res = len_Rs, len_Ws, K_d,K_s,L_d,L_s, distrs,policies
            res = partial_equilibrium_transition(Rs, Ws)
            len_Rs = copy(res[1])
            len_Ws = copy(res[2])
            total_len = sum(abs,len_Rs)+sum(abs,len_Ws)
            abs_len_r = maximum(abs,len_Rs)
            abs_len_w = maximum(abs,len_Ws)
            #=
            push!(past_Rs,Rs)
            push!(past_Ws,Ws)
            push!(past_len_Rs, len_Rs)
            push!(past_len_Ws, len_Ws)
            =#

            # get new_Rs and new_Ws from len_Rs and len_Ws

            ##_new_Rs =     Rs .- len_Rs.*(Rs.-old_Rs) ./(len_Rs.-old_len_Rs)
            #best_new_Rs = Rs .- len_Rs.*(Rs.-best_Rs)./(len_Rs.-best_len_Rs)
            ##new_Rs = (_new_Rs.+best_new_Rs)./2.0
            #new_Rs = best_new_Rs
            new_Rs = Rs.*(1.0.+len_Rs)./(1.0.-len_Rs) .+ 2.0.*(1.0-r_min).*len_Rs./(1.0.-len_Rs)
            #relax_R = relax_r*(1.0-rw_iters/(rw_maxiters+1))
            if old_total_len < total_len
                relax_R /= 2.0
            end
            if maximum(abs,new_Rs.-Rs)*relax_R < gen_tol_x
                relax_R = relax_r
                new_Rs = best_Rs.*(1.0.+best_len_Rs)./(1.0.-best_len_Rs) .+ 2.0.*(1.0-r_min).*best_len_Rs./(1.0.-best_len_Rs)
            end
            new_Rs = min.(max.(r_min, new_Rs), r_max).*relax_R .+ (1-relax_R).*Rs
            new_Rs[3:end-1] = min.(max.(ss_starstar[2], new_Rs[3:end-1]), best_Rs[2])
            #new_Rs[end] = ss_starstar[2]

            ##_new_Ws =     Ws .- len_Ws.*(Ws.-old_Ws) ./(len_Ws.-old_len_Ws)
            #best_new_Ws = Ws .- len_Ws.*(Ws.-best_Ws)./(len_Ws.-best_len_Ws)
            ##new_Ws = (_new_Ws.+best_new_Ws)./2.0
            #new_Ws = best_new_Ws
            new_Ws = Ws.*(1.0.+len_Ws)./(1.0.-len_Ws)
            #relax_W = relax_w*(1.0-rw_iters/(rw_maxiters+1))
            if old_total_len < total_len
                relax_W /= 2.0
            end
            if maximum(abs,new_Ws.-Ws)*relax_W < gen_tol_x
                relax_W = relax_w
                new_Ws = best_Ws.*(1.0.+best_len_Ws)./(1.0.-best_len_Ws)
            end
            new_Ws = min.(max.(w_min, new_Ws), w_max).*relax_W .+ (1-relax_W).*Ws
            new_Ws[2:end-1] = min.(max.(ss_star[3], new_Ws[2:end-1]), ss_starstar[3])
            #new_Ws[end] = ss_starstar[3]

            println_sameline("#$(rw_iters) - total_len:$(round(total_len;digits=6)) - total_len_Rs:$(round(sum(abs,len_Rs);digits=6)) - total_len_Ws:$(round(sum(abs,len_Ws);digits=6)) - abs_len_r:$(round(abs_len_r;digits=6)) - abs_len_w:$(round(abs_len_w;digits=6))")
            plot_iter_results(Rs,Ws, new_Rs,new_Ws, best_Rs,best_Ws, len_Rs,len_Ws, best_len_Rs,best_len_Ws)
            # save best candidates according to different metrics

            if total_len < best_total_len
                best_Rs = copy(Rs)
                best_Ws = copy(Ws)
                best_res = copy(res)
                best_len_Rs = copy(len_Rs)
                best_len_Ws = copy(len_Ws)
                best_total_len = total_len
                best_abs_len_r = abs_len_r
                best_abs_len_w = abs_len_w
                println_sameline("#$(rw_iters) - New best Rs and Ws")
            #end
            else
                is_best_updated = false
                cand_best_Rs = copy(best_Rs)
                cand_best_Ws = copy(best_Ws)
                cand_best_res = copy(best_res)
                cand_best_len_Rs = copy(best_len_Rs)
                cand_best_len_Ws = copy(best_len_Ws)
                cand_best_total_len = best_total_len
                cand_best_abs_len_r = best_abs_len_r
                cand_best_abs_len_w = best_abs_len_w
                Threads.@threads for i=1:T
                    if abs(best_len_Rs[i]) > abs(len_Rs[i])
                        cand_best_Rs[i] = Rs[i]
                        cand_best_Ws[i] = Ws[i]
                        is_best_updated = true
                    elseif abs(best_len_Ws[i]) > abs(len_Ws[i])
                        cand_best_Rs[i] = Rs[i]
                        cand_best_Ws[i] = Ws[i]
                        is_best_updated = true
                    end
                end
                if is_best_updated
                    #println_sameline("#$(rw_iters) - update the candidate for the best")
                    cand_best_res = partial_equilibrium_transition(cand_best_Rs, cand_best_Ws)
                    cand_best_len_Rs = copy(cand_best_res[1])
                    cand_best_len_Ws = copy(cand_best_res[2])
                    cand_best_total_len = sum(abs,cand_best_len_Rs)+sum(abs,cand_best_len_Ws)
                    cand_best_abs_len_r = maximum(abs,cand_best_len_Rs)
                    cand_best_abs_len_w = maximum(abs,cand_best_len_Ws)
                    #=
                    push!(past_Rs,cand_best_Rs)
                    push!(past_Ws,cand_best_Ws)
                    push!(past_len_Rs, cand_best_len_Rs)
                    push!(past_len_Ws, cand_best_len_Ws)
                    =#
                    println_sameline("#$(rw_iters) - cand_total_len:$(round(cand_best_total_len;digits=6)) - cand_total_len_Rs:$(round(sum(abs,cand_best_len_Rs);digits=6)) - cand_total_len_Ws:$(round(sum(abs,cand_best_len_Ws);digits=6)) - cand_abs_len_r:$(round(cand_best_abs_len_r;digits=6)) - cand_abs_len_w:$(round(cand_best_abs_len_w;digits=6))")
                    plot_iter_results(cand_best_Rs,cand_best_Ws, new_Rs,new_Ws, best_Rs,best_Ws, cand_best_len_Rs,cand_best_len_Ws, best_len_Rs,best_len_Ws)
                    if cand_best_total_len < total_len
                        new_Rs = best_Rs.*(1.0.+best_len_Rs)./(1.0.-best_len_Rs) .+ 2.0.*(1.0-r_min).*best_len_Rs./(1.0.-best_len_Rs)
                        if maximum(abs,new_Rs.-best_Rs)*relax_R < gen_tol_x
                            relax_R = relax_r
                            new_Rs = Rs.*(1.0.+len_Rs)./(1.0.-len_Rs) .+ 2.0.*(1.0-r_min).*len_Rs./(1.0.-len_Rs)
                        end
                        new_Rs = min.(max.(r_min, new_Rs), r_max).*relax_R .+ (1-relax_R).*best_Rs
                        new_Rs[3:end-1] = min.(max.(ss_starstar[2], new_Rs[3:end-1]), best_Rs[2])
                        #new_Rs[end] = ss_starstar[2]

                        new_Ws = best_Ws.*(1.0.+best_len_Ws)./(1.0.-best_len_Ws)
                        if maximum(abs,new_Ws.-Ws)*relax_W < gen_tol_x
                            relax_W = relax_w
                            new_Ws = Ws.*(1.0.+len_Ws)./(1.0.-len_Ws)
                        end
                        new_Ws = min.(max.(w_min, new_Ws), w_max).*relax_W .+ (1-relax_W).*best_Ws
                        new_Ws[2:end-1] = min.(max.(ss_star[3], new_Ws[2:end-1]), ss_starstar[3])
                        #new_Ws[end] = ss_starstar[3]
                    end
                    if cand_best_total_len < best_total_len
                        best_Rs = copy(cand_best_Rs)
                        best_Ws = copy(cand_best_Ws)
                        best_res = copy(cand_best_res)
                        best_len_Rs = copy(cand_best_len_Rs)
                        best_len_Ws = copy(cand_best_len_Ws)
                        best_total_len = cand_best_total_len
                        best_abs_len_r = cand_best_abs_len_r
                        best_abs_len_w = cand_best_abs_len_w
                        println_sameline("#$(rw_iters) - New best Rs and Ws")

                    end
                    is_best_updated = false
                end
            end
            # update Rs and Ws with the best candidate
        catch e
            println_sameline(("#$(rw_iters) - half the Rs and Ws", e))

            new_Rs = (old_Rs.+new_Rs.+best_Rs)./3.0
            new_Rs[end] = ss_starstar[2]
            new_Ws = (old_Ws.+new_Ws.+best_Ws)./3.0
            new_Ws[end] = ss_starstar[3]
            total_len = old_total_len
            println_sameline("#$(rw_iters) - total_len:-.------ - total_len_Rs:-.------ - total_len_Ws:-.------ - abs_len_r:-.------ - abs_len_w:-.------")

            Rs = copy(old_Rs)
            Ws = copy(old_Ws)
            old_Rs = copy(old_old_Rs)
            old_Ws = copy(old_old_Ws)
            old_len_Rs = copy(old_old_len_Rs)
            old_len_Ws = copy(old_old_len_Ws)

            plot_iter_results(Rs,Ws, new_Rs,new_Ws, best_Rs,best_Ws, len_Rs,len_Ws, best_len_Rs,best_len_Ws)

            rw_iters -= 1
        end

        if maximum([abs.(new_Rs.-Rs) abs.(new_Ws.-Ws)]) < gen_tol_x
            println_sameline("#$(rw_iters) - Paths has converged, but markets are still not clear, so restart from random deviation to paths")
            new_Rs = Rs .+ len_Rs
            new_Rs = min.(max.(r_min, new_Rs), r_max)
            new_Rs[end] = ss_starstar[2]
            
            new_Ws = Ws .+ len_Ws
            #new_Ws = min.(max.(w_min, new_Ws), w_max)
            new_Ws = min.(max.(min(ss_star[3],ss_starstar[3]), new_Ws), max(ss_star[3],ss_starstar[3]))
            new_Ws[end] = ss_starstar[3]

            plot_iter_results(Rs,Ws, new_Rs,new_Ws, best_Rs,best_Ws, len_Rs,len_Ws, best_len_Rs,best_len_Ws)
        end

        rw_iters += 1

    end

    # return best result
    Rs = copy(best_Rs)
    Ws = copy(best_Ws)
    res = copy(best_res)
    len_Rs = copy(best_len_Rs)
    len_Ws = copy(best_len_Ws)
    total_len = best_total_len
    abs_len_r = best_abs_len_r
    abs_len_w = best_abs_len_w

    println_sameline("#$(rw_iters) - total_len:$(round(total_len;digits=6)) - total_len_Rs:$(round(sum(abs,len_Rs);digits=6)) - total_len_Ws:$(round(sum(abs,len_Ws);digits=6)) - abs_len_r:$(round(abs_len_r;digits=6)) - abs_len_w:$(round(abs_len_w;digits=6))")
    plot_iter_results(Rs,Ws, Rs,Ws, best_Rs,best_Ws, len_Rs,len_Ws, best_len_Rs,best_len_Ws)

    #       1       2       3   4   5   6    7      8
    # res = len_Rs, len_Ws, K_d,K_s,L_d,L_s, distrs,policies

    #       1  2         3   4   5    6      7               8        9                  10
    return [T, lambda_s, Rs, Ws, res, a_max, number_a_nodes, a_nodes, number_asset_grid, asset_grid]

end
