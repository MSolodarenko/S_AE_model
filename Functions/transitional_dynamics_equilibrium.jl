include("print_sameline.jl")

function transitional_dynamics(lambda_s, ss_star, ss_starstar, global_params, global_approx_params, model_params, file_name, runway, GUESS_RS, smoothing, maxiters)

    println_sameline("Preparing initial variables")
    RUNWAY_HALF = Int64(round(runway/2; digits=0))

    lambda_s = [lambda_s; ones(runway).*lambda_s[end]]

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

    ############################
    a_min = 0.01#0.0
    a_max = ss_starstar[1][1]#2*entrep_optimal(maximum(z_m_nodes),maximum(z_w_nodes),r,w, delta, gamma, eta, theta, c_e)[1]

    a_nodes = exp.(collect(range(log(a_min+1); stop=log(a_max+1), length=number_a_nodes))).-1
    ###########################

    len = Inf
    len_T = Inf
    tol = global_params[2]/2.0#1e-4/2.0
    tot_tol = global_params[2]
    abs_tol = global_params[1]
    iters = 0
    relax = 0.05#0.1#
    # NEW ADDITION 25.11
    relax_r = 0.05#relax
    relax_w = 0.1#relax

    number_asset_grid = ss_star[1][2]
    asset_grid = exp.(collect(range(log(a_min+1); stop=log(a_max+1), length=number_asset_grid))).-1

    capital_s_distr_s = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)

    policy_s = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    policy_s[T,:,:,:,:,:] .= ss_starstar[1][4]

    occ_choice_s = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    occ_choice_s[T,:,:,:,:,:] .= ss_starstar[1][22]
    income_s = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    income_s[T,:,:,:,:,:] .= ss_starstar[1][23]
    earnings_s = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    earnings_s[T,:,:,:,:,:] .= ss_starstar[1][24]
    capital_excess_s = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    capital_excess_s[T,:,:,:,:,:] .= ss_starstar[1][25]
    capital_d_s = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    capital_d_s[T,:,:,:,:,:] .= ss_starstar[1][26]
    credit_s = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    credit_s[T,:,:,:,:,:] .= ss_starstar[1][27]
    labour_excess_s = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    labour_excess_s[T,:,:,:,:,:] .= ss_starstar[1][28]
    labour_d_s = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    labour_d_s[T,:,:,:,:,:] .= ss_starstar[1][29]
    labour_s_s = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    labour_s_s[T,:,:,:,:,:] .= ss_starstar[1][30]
    deposit_s = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    deposit_s[T,:,:,:,:,:] .= ss_starstar[1][31]
    output_s = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    output_s[T,:,:,:,:,:] .= ss_starstar[1][32]
    # Step 3
    #
    #r_s = [collect(range(ss_star[2]; stop=ss_starstar[2], length=T-runway)); ones(runway).*ss_starstar[2]]
    #w_s = [collect(range(ss_star[3]; stop=ss_starstar[3], length=T-runway)); ones(runway).*ss_starstar[3]]

    r_s = [(log.(collect(range(exp(-2.0); stop=exp(10.0), length=T-runway))) .- (-2.0)).*(ss_starstar[2] - ss_star[2])./(10 - (-2)) .+ ss_star[2]; ones(runway).*ss_starstar[2]]
    w_s = [(log.(collect(range(exp(-2.0); stop=exp(10.0), length=T-runway))) .- (-2.0)).*(ss_starstar[3] - ss_star[3])./(10 - (-2)) .+ ss_star[3]; ones(runway).*ss_starstar[3]]

    if GUESS_RS
        r_s = [ss_star[2]; (1.0./(collect(range(0.5; stop=2.0, length=T-runway-1))) .- 1.0/0.5).*(ss_starstar[2] - 2*ss_starstar[2])./(1.0/(2.0) - 1.0/0.5).+2*ss_starstar[2]; ones(runway).*ss_starstar[2]]
    end

    p1 = Plots.plot([r_s,ones(T).*ss_star[2],ones(T).*ss_starstar[2]], legend=false)
    p2 = Plots.plot([w_s,ones(T).*ss_star[3],ones(T).*ss_starstar[3]], legend=false)
    display(Plots.plot(p1,p2, layout=(1,2)))

    #throw(error)

    old_rs = copy(r_s)
    old_ws = copy(w_s)

    new_rs = copy(r_s)
    new_ws = copy(w_s)
    K_d = zeros(T)
    K_s = zeros(T)
    L_d = zeros(T)
    L_s = zeros(T)
    len_r = zeros(T)
    len_w = zeros(T)

    old_rs = copy(r_s)
    old_ws = copy(w_s)
    old_len_r = copy(len_r)
    old_len_w = copy(len_w)
    old_len = len
    old_new_rs = copy(new_rs)
    old_new_ws = copy(new_ws)

    old_capital_s_distr_s = copy(capital_s_distr_s)
    old_output_path = [sum(output_s[t,:,:,:,:,:] .* capital_s_distr_s[t,:,:,:,:,:]) for t=1:T]
    old_occ_W_path = [sum( capital_s_distr_s[t,:,:,:,:,:] .* Float64.(occ_choice_s[t,:,:,:,:,:].==1.0) ) for t=1:T]
    old_occ_SP_path = [sum( capital_s_distr_s[t,:,:,:,:,:] .* Float64.(occ_choice_s[t,:,:,:,:,:].==2.0) ) for t=1:T]
    old_occ_EMP_path = [sum( capital_s_distr_s[t,:,:,:,:,:] .* Float64.(occ_choice_s[t,:,:,:,:,:].==3.0) ) for t=1:T]

    #println_sameline("#i - T - len_r_T - len_w_T - len_r_max - len_w_max - len_r_sum - len_w_sum - len_sum - len_rs_sum - len_ws_sum                                               ")

    while len > tol && iters <= maxiters

        asprimes = zeros(T, number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            asprimeT = Schumaker(ss_starstar[1][3],ss_starstar[1][4][:,u_i,zeta_i,alpha_m_i,alpha_w_i]; extrapolation = (Linear,Linear))
            asprimes[T,:,u_i,zeta_i,alpha_m_i,alpha_w_i] .= evaluate.(asprimeT, a_nodes)
        end

        values = zeros(T, number_a_nodes,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        values[T,:,:,:,:,:] .= ss_starstar[1][33]

        for t = T-1:-1:1

            occ_choice, income, earnings, capital_excess, capital_d, credit, labour_excess, labour_d, labour_s, deposit, output = compute_income_profile(a_nodes,number_a_nodes,r_s[t],w_s[t], number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, lambda_s[t], delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)

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
                #print("$(u_i),$(zeta_i),$(alpha_m_i),$(alpha_w_i)\r")
                values[t,:,u_i,zeta_i,alpha_m_i,alpha_w_i], asprimes[t,:,u_i,zeta_i,alpha_m_i,alpha_w_i] = new_val_and_a1(values[t+1,:,u_i,zeta_i,alpha_m_i,alpha_w_i],asprimes[t+1,:,u_i,zeta_i,alpha_m_i,alpha_w_i], income[:,u_i,zeta_i,alpha_m_i,alpha_w_i], value_tran_rhs[:,alpha_m_i,alpha_w_i],expectation_value_death_rhs, u_i,alpha_m_i,alpha_w_i ,a_min,a_max,a_nodes, 0.0,0.0, number_a_nodes, beta, p_alpha, P_u, number_u_nodes, CRRA)
                #println_sameline(new_aprime_nodes[:,u_i,zeta_i,alpha_m_i,alpha_w_i])
                #throw(error)
            end

            print("#$(iters) - $(T) - $(round(len_T;digits = 4)) - $(round(len;digits = 4)) - Policy iteration: $(t)/$(T)                       \r")
        end

        # Step 5

        policy1 = copy(ss_star[1][4])
        density_distr1 = copy(ss_star[1][5])
        Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            asprime1 = Schumaker(ss_star[1][3],ss_star[1][4][:,u_i,zeta_i,alpha_m_i,alpha_w_i]; extrapolation = (Linear,Linear))
            policy1[:,u_i,zeta_i,alpha_m_i,alpha_w_i] .= evaluate.(asprime1, asset_grid)

            density_distr11 = Schumaker(ss_star[1][3],cumsum(ss_star[1][5][:,u_i,zeta_i,alpha_m_i,alpha_w_i]); extrapolation = (Linear,Linear))
            density_distr1[:,u_i,zeta_i,alpha_m_i,alpha_w_i] .= diff([0;evaluate.(density_distr11, asset_grid)])
            density_distr1[:,u_i,zeta_i,alpha_m_i,alpha_w_i] .*= sum(ss_star[1][5][:,u_i,zeta_i,alpha_m_i,alpha_w_i])/sum(density_distr1[:,u_i,zeta_i,alpha_m_i,alpha_w_i])
        end

        capital_s_distr_s = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        capital_s_distr_s[1,:,:,:,:,:] .= density_distr1

        occ_choice_s[1,:,:,:,:,:], income_s[1,:,:,:,:,:], earnings_s[1,:,:,:,:,:], capital_excess_s[1,:,:,:,:,:], capital_d_s[1,:,:,:,:,:], credit_s[1,:,:,:,:,:], labour_excess_s[1,:,:,:,:,:], labour_d_s[1,:,:,:,:,:], labour_s_s[1,:,:,:,:,:], deposit_s[1,:,:,:,:,:], output_s[1,:,:,:,:,:] = compute_income_profile(asset_grid,number_asset_grid,r_s[1],w_s[1], number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, lambda_s[1], delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)

        K_d[1], K_s[1], L_d[1], L_s[1], len_r[1], len_w[1] = find_aggregate_capital_labour_demand_supply(number_asset_grid,asset_grid,policy1,density_distr1,r_s[1],w_s[1], labour_excess_s[1,:,:,:,:,:], labour_d_s[1,:,:,:,:,:], labour_s_s[1,:,:,:,:,:], capital_excess_s[1,:,:,:,:,:], capital_d_s[1,:,:,:,:,:])
        K_d[1] = sum(density_distr1.*capital_d_s[1,:,:,:,:,:])
        # and K_supply for previous distribution
        # NEW ADDITION 29/11
        K_s[1] = sum( sum(density_distr1; dims=2:5)[:,1,1,1,1] .* asset_grid)
        #K_s[1] = sum(density_distr1.*policy1)

        len_r[1] = (K_d[1]-K_s[1])/(K_d[1]+K_s[1])
        # NEW ADDITION 25.11
        len_w[1] = (L_d[1]-L_s[1])/(L_d[1]+L_s[1])

        if true#iters == 0
            new_rs[1] = r_s[1]*(1.0+len_r[1])/(1.0-len_r[1]) + 2.0*(1.0-r_min)*len_r[1]/(1.0-len_r[1])
            new_rs[1] = min(max(r_min, new_rs[1]), r_max)*relax + r_s[1]*(1-relax)
            new_ws[1] = w_s[1]*(1.0+len_w[1])/(1.0-len_w[1])
            new_ws[1] = min(max(w_min, new_ws[1]), w_max)*relax + w_s[1]*(1-relax)
        else
            new_rs[1] = r_s[1] - len_r[1]*(r_s[1]-old_rs[1])/(len_r[1]-old_len_r[1])
            new_rs[1] = min(max(r_min, new_rs[1]), r_max)*relax + r_s[1]*(1-relax)
            new_ws[1] = w_s[1] - len_w[1]*(w_s[1]-old_ws[1])/(len_w[1]-old_len_w[1])
            new_ws[1] = min(max(w_min, new_ws[1]), w_max)*relax + w_s[1]*(1-relax)
        end
        # NEW ADDITION 25/11
        new_rs[1] = ss_star[2]
        new_ws[1] = ss_star[3]

        for t=1:T-1

            a1_indices = Array{Int64}(undef,number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
            lottery_prob_1 = zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
            lottery_prob_2 = zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
            Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
                try
                    policy_function = Schumaker(a_nodes,asprimes[t,:,u_i,zeta_i,alpha_m_i,alpha_w_i]; extrapolation = (Linear,Linear))
                    for a_i in 1:number_asset_grid
                        a1 = evaluate(policy_function, asset_grid[a_i])
                        policy_s[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = a1
                        j_1 = sum(a1 .>= asset_grid)
                        j = j_1 + 1

                        a1_indices[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = j_1
                        if j <= number_asset_grid && j_1 >= 1
                            lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (asset_grid[j] - a1             )/(asset_grid[j]-asset_grid[j_1])
                            lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (a1            - asset_grid[j_1])/(asset_grid[j]-asset_grid[j_1])
                        elseif j_1 == number_asset_grid
                            lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 1.0
                            lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 0.0
                        end
                    end
                catch e
                    println_sameline("Trouble is in loop 1")
                    throw(e)
                end
            end

            distr_a1_z0 = copy(capital_s_distr_s[t,:,:,:,:,:])
            distr_a1_z0 .= 0.0
            Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
                try
                    non_zero_asset_grid_iters = findall(!iszero,capital_s_distr_s[t,:,u_i,zeta_i,alpha_m_i,alpha_w_i])
                    for a_i in non_zero_asset_grid_iters#1:number_asset_grid#
                        j_1 = a1_indices[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        distr_a1_z0[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += capital_s_distr_s[t,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        if j_1 != number_asset_grid
                            distr_a1_z0[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += capital_s_distr_s[t,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        end
                    end

                catch e
                    println_sameline("Trouble is in loop 2")
                    throw(e)
                end
            end

            new_distr = copy(capital_s_distr_s[t,:,:,:,:,:])
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

            new_distr2 = copy(capital_s_distr_s[t,:,:,:,:,:])
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

            capital_s_distr_s[t+1,:,:,:,:,:] .= new_distr

            occ_choice_s[t+1,:,:,:,:,:], income_s[t+1,:,:,:,:,:], earnings_s[t+1,:,:,:,:,:], capital_excess_s[t+1,:,:,:,:,:], capital_d_s[t+1,:,:,:,:,:], credit_s[t+1,:,:,:,:,:], labour_excess_s[t+1,:,:,:,:,:], labour_d_s[t+1,:,:,:,:,:], labour_s_s[t+1,:,:,:,:,:], deposit_s[t+1,:,:,:,:,:], output_s[t+1,:,:,:,:,:] = compute_income_profile(asset_grid,number_asset_grid,r_s[t+1],w_s[t+1], number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, lambda_s[t+1], delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)

            K_d[t+1], K_s[t+1], L_d[t+1], L_s[t+1], len_r[t+1], len_w[t+1] = find_aggregate_capital_labour_demand_supply(number_asset_grid,asset_grid,policy_s[t,:,:,:,:,:],capital_s_distr_s[t+1,:,:,:,:,:],r_s[t+1],w_s[t+1], labour_excess_s[t+1,:,:,:,:,:], labour_d_s[t+1,:,:,:,:,:], labour_s_s[t+1,:,:,:,:,:], capital_excess_s[t+1,:,:,:,:,:], capital_d_s[t+1,:,:,:,:,:])
            # calculate K_demand for current distribution
            K_d[t+1] = sum(capital_s_distr_s[t+1,:,:,:,:,:].*capital_d_s[t+1,:,:,:,:,:])
            # NEW ADDITION 29/11
            K_s[t+1] = sum( sum(capital_s_distr_s[t+1,:,:,:,:,:]; dims=2:5)[:,1,1,1,1] .* asset_grid)
            # and K_supply for previous distribution
            #K_s[t+1] = sum(capital_s_distr_s[t,:,:,:,:,:].*policy_s[t,:,:,:,:,:])

            #=
            if t > 1
                K_d_dd = K_d[t+1]-2*K_d[t]+K_d[t-1]
                K_s_dd = K_s[t+1]-2*K_s[t]+K_s[t-1]
                if K_d_dd > 0.0
                    if K_s_dd < 0.0
                        K_s[t+1] = 2*K_s[t] - K_s[t-1]
                    end
                else#K_d_dd < 0.0
                    if K_s_dd > 0.0
                        K_s[t+1] = 2*K_s[t] - K_s[t-1]
                    end
                end
            else
                if K_d[t+1]-K_d[t] > 0.0
                    if K_s[t+1]-K_s[t] < 0.0
                        K_s[t+1] = K_s[t]
                    end
                else#K_d[t+1]-K_d[t] < 0.0
                    if K_s[t+1]-K_s[t] > 0.0
                        K_s[t+1] = K_s[t]
                    end
                end
            end
            =#

            len_r[t+1] = (K_d[t+1]-K_s[t+1])/(K_d[t+1]+K_s[t+1])
            # NEW ADDITION 25.11
            len_w[t+1] = (L_d[t+1]-L_s[t+1])/(L_d[t+1]+L_s[t+1])
            #
            # smoothing
            if (smoothing && t>1) || (!smoothing && t>T-runway)
                al_r = 0.05#0.1#0.15#1.0 - 0.05#t/T#
                al_w = 0.25#0.5#0.1#1.0 - 0.05#t/T#
                if t>T-runway
                    al_r = (1.0 - (t-T+runway)/runway)*al_r
                    al_w = (1.0 - (t-T+runway)/runway)*al_w
                end
                len_r[t+1] = al_r*len_r[t+1] + (1.0-al_r)*len_r[t]
                len_w[t+1] = al_w*len_w[t+1] + (1.0-al_w)*len_w[t]
            end

            if true#iters == 0
                new_rs[t+1] = r_s[t+1]*(1.0+len_r[t+1])/(1.0-len_r[t+1]) + 2.0*(1.0-r_min)*len_r[t+1]/(1.0-len_r[t+1])
                new_ws[t+1] = w_s[t+1]*(1.0+len_w[t+1])/(1.0-len_w[t+1])

                new_rs[t+1] = min(max(r_min, new_rs[t+1]), r_max)*relax_r + (1-relax_r)*r_s[t+1]
                new_ws[t+1] = min(max(w_min, new_ws[t+1]), w_max)*relax_w + (1-relax_w)*w_s[t+1]
            else #not working!
                relax = 1.0
                new_rs[t+1] = r_s[t+1] - len_r[t+1]*(r_s[t+1]-old_rs[t+1])/(len_r[t+1]-old_len_r[t+1])
                new_rs[t+1] = min(max(r_min, new_rs[t+1]), r_max)*relax + (1-relax)*r_s[t+1]
                new_ws[t+1] = w_s[t+1] - len_w[t+1]*(w_s[t+1]-old_ws[t+1])/(len_w[t+1]-old_len_w[t+1])
                new_ws[t+1] = min(max(w_min, new_ws[t+1]), w_max)*relax + (1-relax)*w_s[t+1]
            end

            p1 = Plots.plot([r_s[1:t],new_rs[1:t],ones(T)[1:t].*ss_star[2],ones(T)[1:t].*ss_starstar[2]],legend=false)
            p2 = Plots.plot([w_s[1:t],new_ws[1:t],ones(T)[1:t].*ss_star[3],ones(T)[1:t].*ss_starstar[3]],legend=false)

            p3 = Plots.plot([K_d[1:t],K_s[1:t], ones(T)[1:t].*(ss_star[1][6].+ss_star[1][7])./2.0, ones(T)[1:t].*(ss_starstar[1][6].+ss_starstar[1][7])./2.0],legend=false)
            p4 = Plots.plot([L_d[1:t],L_s[1:t], ones(T)[1:t].*(ss_star[1][8].+ss_star[1][9])./2.0, ones(T)[1:t].*(ss_starstar[1][8].+ss_starstar[1][9])./2.0],legend=false)

            display(Plots.plot(p1,p2,p3,p4,layout=(2,2)))

            #=
            distr_label = permutedims(["$(u),$(zeta),$(alpha_m),$(alpha_w)" for u in 1:number_u_nodes for zeta in 1:number_zeta_nodes for alpha_m in 1:number_alpha_m_nodes for alpha_w in 1:number_alpha_w_nodes])
            plot_distr = [capital_s_distr_s[t+1,:,u,zeta,alpha_m,alpha_w] for u in 1:number_u_nodes for zeta in 1:number_zeta_nodes for alpha_m in 1:number_alpha_m_nodes for alpha_w in 1:number_alpha_w_nodes]
            display(Plots.plot(plot_distr,title="$(iters)",label=distr_label,xrotation=-45,legend=:outertopleft,foreground_color_legend = nothing))
            =#

            print("#$(iters) - $(T) - $(round(len_T;digits = 4)) - $(round(len;digits = 4)) - Simulation of the history: t$(t)/$(T) - $(sum(capital_s_distr_s[t+1,:,:,:,:,:]))\r")
        end

        len_T = maximum(abs, [new_rs[T]-ss_starstar[2], new_ws[T]-ss_starstar[3]])

        len = maximum(abs, [new_rs-r_s; new_ws-w_s])
        #len = maximum(abs, [len_r; len_w])
        len_sum = sum(abs, [len_r; len_w])
        #i - T - len_r_T - len_w_T - len_r_max - len_w_max - len_r_sum - len_w_sum - len_sum - len_rs_sum - len_ws_sum
        #println_sameline("#$(iters) - $(T) - $(round(abs(new_rs[T]-ss_starstar[2]);digits = 6)) - $(round(abs(new_ws[T]-ss_starstar[3]);digits = 6)) - $(round(maximum(abs, new_rs-r_s);digits = 6)) - $(round(maximum(abs, new_ws-w_s);digits = 6)) - $(round(sum(abs, len_r);digits = 6)) - $(round(sum(abs, len_w);digits = 6)) - $(round(len_sum;digits = 6)) - $(round(sum(abs, new_rs.-r_s);digits = 6)) - $(round(sum(abs, new_ws.-w_s);digits = 6))                                               ")

        LENDIGITS = 6
        tot_cap_excess = round(sum(abs, len_r); digits = LENDIGITS)
        abs_cap_excess = round(maximum(abs, len_r); digits = LENDIGITS)
        tot_lab_excess = round(sum(abs, len_w); digits = LENDIGITS)
        abs_lab_excess = round(maximum(abs, len_w); digits = LENDIGITS)
        tot_excess = round(tot_cap_excess + tot_lab_excess; digits = LENDIGITS)
        abs_excess = round(maximum(abs, [abs_cap_excess; abs_lab_excess]); digits = LENDIGITS)
        println_sameline("#$(iters) - T$(T) - MarketExcess - $(tot_excess) ($(abs_excess)) - Cap - $(tot_cap_excess) ($(abs_cap_excess)) - Lab - $(tot_lab_excess) ($(abs_lab_excess))")
        write(file_name,"#$(iters) - T$(T) - MarketExcess - $(tot_excess) ($(abs_excess)) - Cap - $(tot_cap_excess) ($(abs_cap_excess)) - Lab - $(tot_lab_excess) ($(abs_lab_excess))\n")

        tot_r_len = round(sum(abs, new_rs-r_s); digits = LENDIGITS)
        abs_r_len = round(maximum(abs, new_rs-r_s); digits = LENDIGITS)
        tot_w_len = round(sum(abs, new_ws-w_s); digits = LENDIGITS)
        abs_w_len = round(maximum(abs, new_ws-w_s); digits = LENDIGITS)
        tot_rw_len = round(tot_r_len+tot_w_len; digits = LENDIGITS)
        abs_rw_len = round(maximum(abs, [abs_r_len; abs_w_len]); digits = LENDIGITS)
        println_sameline("          - RW - $(tot_rw_len) ($(abs_rw_len)) - Rs - $(tot_r_len) ($(abs_r_len)) - Ws - $(tot_w_len) ($(abs_w_len))")
        write(file_name,"          - RW - $(tot_rw_len) ($(abs_rw_len)) - Rs - $(tot_r_len) ($(abs_r_len)) - Ws - $(tot_w_len) ($(abs_w_len))\n")

        tot_distr_len = round(sum(abs, capital_s_distr_s-old_capital_s_distr_s); digits = LENDIGITS)
        abs_distr_len = round(maximum(abs, capital_s_distr_s-old_capital_s_distr_s); digits = LENDIGITS)
        println_sameline("          - Distr - $(tot_distr_len) ($(abs_distr_len))")
        write(file_name,"          - Distr - $(tot_distr_len) ($(abs_distr_len))\n")

        output_path = [sum(output_s[t,:,:,:,:,:] .* capital_s_distr_s[t,:,:,:,:,:]) for t=1:T]
        tot_output_len = round(sum(abs, output_path-old_output_path); digits = LENDIGITS)
        abs_output_len = round(maximum(abs, output_path-old_output_path); digits = LENDIGITS)
        println_sameline("          - Output - $(tot_output_len) ($(abs_output_len))")
        write(file_name,"          - Output - $(tot_output_len) ($(abs_output_len))\n")

        occ_W_path = [sum( capital_s_distr_s[t,:,:,:,:,:] .* Float64.(occ_choice_s[t,:,:,:,:,:].==1.0) ) for t=1:T]
        occ_SP_path = [sum( capital_s_distr_s[t,:,:,:,:,:] .* Float64.(occ_choice_s[t,:,:,:,:,:].==2.0) ) for t=1:T]
        occ_EMP_path = [sum( capital_s_distr_s[t,:,:,:,:,:] .* Float64.(occ_choice_s[t,:,:,:,:,:].==3.0) ) for t=1:T]
        tot_occ_W_len = round(sum(abs, occ_W_path-old_occ_W_path); digits = LENDIGITS)
        abs_occ_W_len = round(maximum(abs, occ_W_path-old_occ_W_path); digits = LENDIGITS)
        tot_occ_SP_len = round(sum(abs, occ_SP_path-old_occ_SP_path); digits = LENDIGITS)
        abs_occ_SP_len = round(maximum(abs, occ_SP_path-old_occ_SP_path); digits = LENDIGITS)
        tot_occ_EMP_len = round(sum(abs, occ_EMP_path-old_occ_EMP_path); digits = LENDIGITS)
        abs_occ_EMP_len = round(maximum(abs, occ_EMP_path-old_occ_EMP_path); digits = LENDIGITS)
        tot_occ_len = round(tot_occ_W_len+tot_occ_SP_len+tot_occ_EMP_len; digits = LENDIGITS)
        abs_occ_len = round(maximum(abs, [abs_occ_W_len; abs_occ_SP_len; abs_occ_EMP_len]); digits = LENDIGITS)
        println_sameline("          - Occ - $(tot_occ_len) ($(abs_occ_len)) - W - $(tot_occ_W_len) ($(abs_occ_W_len)) - SP - $(tot_occ_SP_len) ($(abs_occ_SP_len)) - EMP - $(tot_occ_EMP_len) ($(abs_occ_EMP_len))")
        write(file_name,"          - Occ - $(tot_occ_len) ($(abs_occ_len)) - W - $(tot_occ_W_len) ($(abs_occ_W_len)) - SP - $(tot_occ_SP_len) ($(abs_occ_SP_len)) - EMP - $(tot_occ_EMP_len) ($(abs_occ_EMP_len))\n")

        len = len_sum

        if iters < maxiters
            new_rs[end-RUNWAY_HALF:end] .= ss_starstar[2]
            new_ws[end-RUNWAY_HALF:end] .= ss_starstar[3]

            if tot_excess <= tot_tol
                println_sameline("Total market clearance converged at #$(iters)")
                write(file_name,"Total market clearance converged at #$(iters)\n")
            else
                if tot_cap_excess <= tot_tol
                    println_sameline("Total capital market clearance converged at #$(iters)")
                    write(file_name,"Total capital market clearance converged at #$(iters)\n")
                end
                if tot_lab_excess <= tot_tol
                    println_sameline("Total labour market clearance converged at #$(iters)")
                    write(file_name,"Total labour market clearance converged at #$(iters)\n")
                end
            end
            if abs_excess <= abs_tol
                println_sameline("Absolute market clearance converged at #$(iters)")
                write(file_name,"Absolute market clearance converged at #$(iters)\n")
            else
                if abs_cap_excess <= abs_tol
                    println_sameline("Absolute capital market clearance converged at #$(iters)")
                    write(file_name,"Absolute capital market clearance converged at #$(iters)\n")
                end
                if abs_lab_excess <= abs_tol
                    println_sameline("Absolute labour market clearance converged at #$(iters)")
                    write(file_name,"Absolute labour market clearance converged at #$(iters)\n")
                end
            end

            if tot_rw_len <= tot_tol
                println_sameline("Total factor prices paths converged at #$(iters)")
                write(file_name,"Total factor prices paths converged at #$(iters)\n")
            else
                if tot_r_len <= tot_tol
                    println_sameline("Total interest rate path converged at #$(iters)")
                    write(file_name,"Total interest rate path converged at #$(iters)\n")
                end
                if tot_w_len <= tot_tol
                    println_sameline("Total wage path converged at #$(iters)")
                    write(file_name,"Total wage path converged at #$(iters)\n")
                end
            end
            if abs_rw_len <= abs_tol
                println_sameline("Absolute factor prices paths converged at #$(iters)")
                write(file_name,"Absolute factor prices paths converged at #$(iters)\n")
            else
                if abs_r_len <= abs_tol
                    println_sameline("Absolute interest rate path converged at #$(iters)")
                    write(file_name,"Absolute interest rate path converged at #$(iters)\n")
                end
                if abs_w_len <= abs_tol
                    println_sameline("Absolute wage path converged at #$(iters)")
                    write(file_name,"Absolute wage path converged at #$(iters)\n")
                end
            end

            if tot_distr_len <= tot_tol
                println_sameline("Total distribution path convereged at #$(iters)")
                write(file_name,"Total distribution path convereged at #$(iters)\n")
            end
            if abs_distr_len <= abs_tol
                println_sameline("Total distribution path convereged at #$(iters)")
                write(file_name,"Total distribution path convereged at #$(iters)\n")
            end

            if tot_output_len <= tot_tol
                println_sameline("Total output path convereged at #$(iters)")
                write(file_name,"Total output path convereged at #$(iters)\n")
            end
            if abs_output_len <= abs_tol
                println_sameline("Total output path convereged at #$(iters)")
                write(file_name,"Total output path convereged at #$(iters)\n")
            end

            if tot_occ_len <= tot_tol
                println_sameline("Total occupation paths converged at #$(iters)")
                write(file_name,"Total occupation paths converged at #$(iters)\n")
            else
                if tot_occ_W_len <= tot_tol
                    println_sameline("Total Worker occupation path converged at #$(iters)")
                    write(file_name,"Total Worker occupation path converged at #$(iters)\n")
                end
                if tot_occ_SP_len <= tot_tol
                    println_sameline("Total Sole Proprietor occupation path converged at #$(iters)")
                    write(file_name,"Total Sole Proprietor occupation path converged at #$(iters)\n")
                end
                if tot_occ_EMP_len <= tot_tol
                    println_sameline("Total Employer occupation path converged at #$(iters)")
                    write(file_name,"Total Employer occupation path converged at #$(iters)\n")
                end
            end
            if abs_occ_len <= abs_tol
                println_sameline("Absolute occupation paths converged at #$(iters)")
                write(file_name,"Absolute occupation paths converged at #$(iters)\n")
            else
                if abs_occ_W_len <= abs_tol
                    println_sameline("Absolute Worker occupation path converged at #$(iters)")
                    write(file_name,"Absolute Worker occupation path converged at #$(iters)\n")
                end
                if abs_occ_SP_len <= abs_tol
                    println_sameline("Absolute Sole Proprietor occupation path converged at #$(iters)")
                    write(file_name,"Absolute Sole Proprietor occupation path converged at #$(iters)\n")
                end
                if abs_occ_EMP_len <= abs_tol
                    println_sameline("Absolute Employer occupation path converged at #$(iters)")
                    write(file_name,"Absolute Employer occupation path converged at #$(iters)\n")
                end
            end
        end

        old_rs = copy(r_s)
        old_ws = copy(w_s)
        old_len_r = copy(len_r)
        old_len_w = copy(len_w)
        old_len = len
        old_new_rs = copy(new_rs)
        old_new_ws = copy(new_ws)

        r_s = copy(new_rs)
        w_s = copy(new_ws)

        old_capital_s_distr_s = copy(capital_s_distr_s)
        old_output_path = copy(output_path)
        old_occ_W_path = copy(occ_W_path)
        old_occ_SP_path = copy(occ_SP_path)
        old_occ_EMP_path = copy(occ_EMP_path)

        p1 = Plots.plot([r_s,ones(T).*ss_star[2],ones(T).*ss_starstar[2]],legend=false)
        p2 = Plots.plot([w_s,ones(T).*ss_star[3],ones(T).*ss_starstar[3]],legend=false)
        display(Plots.plot(p1,p2,layout=(1,2)))

        iters += 1

    end

    #return r_s, w_s, asset_grid, capital_s_distr_s, a_nodes,asprimes
    #       1  2         3    4    5           6                  7         8             9         10          11                12           13        14               15          16          17         18        19   20   21   22   23     24
    return [T, lambda_s, r_s, w_s, asset_grid, capital_s_distr_s, policy_s, occ_choice_s, income_s, earnings_s, capital_excess_s, capital_d_s, credit_s, labour_excess_s, labour_d_s, labour_s_s, deposit_s, output_s, K_d, K_s, L_d, L_s, len_r, len_w]
end
