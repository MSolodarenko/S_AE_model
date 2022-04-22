include("AllubErosa.jl")

function find_policy_fixed_occ(a_min,a_max,a_nodes,r,w, income,earnings, val_tol, number_a_nodes, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, beta, p_alpha, P_zeta, P_u, stat_P_u, P_alpha, number_u_nodes, crra)

    #println_sameline("V $(value[:,end,end,end,end])")
    print_sameline("Initialise aprime_indices and aprime_nodes")
    # guess values and indices aprime
    aprime_nodes = Array{Any}(undef,3)
    for occ = 1:3
        aprime_nodes[occ] = 0.66.*income[occ]
    end

    print_sameline("Initialise value function")
    # guess the initial value for value_function
    value = Array{Any}(undef,3)
    for occ = 1:3
        #value = utility.(earnings, crra)./(1.0-beta)
        value[occ] = utility.(income[occ].-aprime_nodes[occ], crra)./(1.0-beta)
    end

    value_from_file_flag = true
    try
        #throw(error)
        path = "$(@__DIR__)/val_aprime/"
        if Sys.iswindows()
            path = "$(@__DIR__)\\val_aprime\\"
        end
        @load "$(path)val_aprime_$(number_a_nodes)_$(number_u_nodes)_$(number_zeta_nodes)_$(number_alpha_m_nodes)_$(number_alpha_w_nodes)_fixed_occ.jld2" local_value local_aprime_nodes
        value = copy(local_value)
        aprime_nodes = copy(local_aprime_nodes)
        print_sameline("Initialise value function from file")
    catch e
        value_from_file_flag = false
        print_sameline("Initialise value function from scratch")
    end

    val_Delta   = 0.5#0.05  # update parameter for value function iteration
    val_maxiters= 1000#250#500#
    if value_from_file_flag
        val_maxiters= 300
    end
    val_len     = Inf
    val_iters   = 0
    aprime_len = Inf
    print_sameline("Start of the main VFI loop")

    new_value = copy(value)
    new_aprime_nodes = copy(aprime_nodes)

    # future payoffs if no change in fixed effects (no death)
    value_tran = Array{Any}(undef,3)
    for occ = 1:3
        value_tran[occ] = zeros(number_a_nodes,number_u_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    end
    # future expected payoffs conditional on no death
    expectation_value_death = Array{Any}(undef,3)
    for occ = 1:3
        expectation_value_death[occ] = zeros(number_a_nodes)
    end

    old_val_len = copy(val_len)
    stable = true

    while val_len > val_tol && val_iters < val_maxiters

        # compute future payoffs if no change in fixed effects (no death)
        value_tran .*= 0.0
        for occ = 1:3
            for u_i in 1:number_u_nodes
                for alpha_m_i in 1:number_alpha_m_nodes
                    for alpha_w_i in 1:number_alpha_w_nodes
                        for zeta_prime_i in 1:number_zeta_nodes
                            value_tran[occ][:,u_i,alpha_m_i,alpha_w_i] .+= P_zeta[zeta_prime_i].*value[occ][:,u_i,zeta_prime_i,alpha_m_i,alpha_w_i]
                        end
                    end
                end
            end
        end

        value_tran_rhs = Array{Any}(undef,3)
        for occ = 1:3
            value_tran_rhs[occ] = Array{Any}(undef, number_u_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
            for u_prime_i in 1:number_u_nodes
                for alpha_m_i in 1:number_alpha_m_nodes
                    for alpha_w_i in 1:number_alpha_w_nodes
                        value_tran_rhs[occ][u_prime_i,alpha_m_i,alpha_w_i] = Schumaker(a_nodes,value_tran[occ][:,u_prime_i,alpha_m_i,alpha_w_i]; extrapolation = (Linear,Linear))
                    end
                end
            end
        end

        expectation_value_death .*= 0.0
        for occ = 1:3
            for aprime_i in 1:number_a_nodes
                for alpha_m_prime_i in 1:number_alpha_m_nodes
                    for alpha_w_prime_i in 1:number_alpha_w_nodes
                        expectation_value_death[occ][aprime_i] += beta*p_alpha*P_alpha[alpha_m_prime_i,alpha_w_prime_i] * sum(value_tran[occ][aprime_i,:,alpha_m_prime_i,alpha_w_prime_i].*stat_P_u)
                    end
                end
            end
        end

        expectation_value_death_rhs = Array{Any}(undef,3)
        for occ = 1:3
            expectation_value_death_rhs[occ] = Schumaker(a_nodes,expectation_value_death[occ]; extrapolation = (Linear,Linear))
        end

        # loop for finding new_value and new_aprime_indices
        Threads.@threads for (occ,(u_i,(zeta_i,(alpha_m_i,alpha_w_i)))) in collect(Iterators.product(1:3,Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes)))))
            #print_sameline("$(u_i),$(zeta_i),$(alpha_m_i),$(alpha_w_i)\r")
            new_value[occ][:,u_i,zeta_i,alpha_m_i,alpha_w_i], new_aprime_nodes[occ][:,u_i,zeta_i,alpha_m_i,alpha_w_i] = new_val_and_a1(value[occ][:,u_i,zeta_i,alpha_m_i,alpha_w_i],aprime_nodes[occ][:,u_i,zeta_i,alpha_m_i,alpha_w_i], income[occ][:,u_i,zeta_i,alpha_m_i,alpha_w_i], value_tran_rhs[occ][:,alpha_m_i,alpha_w_i],expectation_value_death_rhs[occ], u_i,alpha_m_i,alpha_w_i ,a_min,a_max,a_nodes, aprime_len,val_len, number_a_nodes, beta, p_alpha, P_u, number_u_nodes, crra)
            #println_sameline(new_aprime_nodes[:,u_i,zeta_i,alpha_m_i,alpha_w_i])
            #throw(error)
        end

        val_iters += 1
        b_lowerbar = zeros(3)
        b_upperbar = zeros(3)
        val_len = 0.0
        val_sumlen = 0.0
        aprime_len = 0.0
        aprime_sumlen = 0.0
        for occ = 1:3
            b_lowerbar[occ] = (beta/(1-beta))*minimum(new_value[occ] - value[occ])
            b_upperbar[occ] = (beta/(1-beta))*maximum(new_value[occ] - value[occ])

            val_len = max(val_len, maximum(abs,new_value[occ]-value[occ])/maximum(abs,new_value[occ]))

            val_sumlen = max(val_sumlen, sum(abs,new_value[occ]-value[occ])/maximum(abs,new_value[occ]))

            aprime_len = max(aprime_len, maximum(abs,new_aprime_nodes[occ]-aprime_nodes[occ])/maximum(abs,new_aprime_nodes[occ]))
            aprime_sumlen = max(aprime_len, sum(abs,new_aprime_nodes[occ]-aprime_nodes[occ])/maximum(abs,new_aprime_nodes[occ]))
        end

        if text_output
            print_sameline("VF#$(val_iters) - err: $(round(val_len;digits=12)) b_l:$(round.(b_lowerbar;digits=4)) b_u:$(round.(b_upperbar;digits=4)), sum_err:$(round(val_sumlen;digits=9)), a_err:$(round(aprime_len;digits=4)), a_sumerr:$(round(aprime_sumlen;digits=4))")
        end
        if val_len > val_tol*5 && old_val_len > val_len && stable
            for occ = 1:3
                if b_lowerbar[occ] > -10000000+1#=-Inf=# && b_upperbar[occ] < 10000000-1#=Inf=# #&& b_lowerbar < 0.0 && b_upperbar > 0.0
                    value[occ] .= new_value[occ] .+ (b_lowerbar[occ] + b_upperbar[occ])/2
                else
                    value[occ] .= new_value[occ]
                end
                aprime_nodes[occ] .= new_aprime_nodes[occ]
            end
        else
            stable = false
            for occ = 1:3
                value[occ] .= new_value[occ].*val_Delta .+ value[occ].*(1-val_Delta)
                aprime_nodes[occ] .= new_aprime_nodes[occ].*val_Delta .+ aprime_nodes[occ].*(1-val_Delta)
            end
        end

        old_val_len = copy(val_len)

        if fig_output
            #=
            # recover aprime_nodes
            for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            aprime_nodes[:,u_i,zeta_i,alpha_m_i,alpha_w_i] = a_nodes[aprime_indices[:,u_i,zeta_i,alpha_m_i,alpha_w_i]]
            end
            =#

            #=
            labels = permutedims(["$(u),$(zeta),$(alpha_m),$(alpha_w)" for u in 1:number_u_nodes for zeta in 1:number_zeta_nodes for alpha_m in 1:number_alpha_m_nodes for alpha_w in 1:number_alpha_w_nodes])
            p1 = plot(#=a_nodes,=#[value[:,u,zeta,alpha_m,alpha_w] for u in 1:number_u_nodes for zeta in 1:number_zeta_nodes for alpha_m in 1:number_alpha_m_nodes for alpha_w in 1:number_alpha_w_nodes], label=labels, legend=:outertopleft)
            p2 = plot(#=a_nodes,=#[aprime_nodes[:,u,zeta,alpha_m,alpha_w] for u in 1:number_u_nodes for zeta in 1:number_zeta_nodes for alpha_m in 1:number_alpha_m_nodes for alpha_w in 1:number_alpha_w_nodes], label=labels, legend=false)
            display(plot(p1,p2, title = "VF and a' - $(val_iters)"))
            #println_sameline(aprime_nodes[:,1,1,1,1])
            =#
        end

    end

    if text_output
        print_sameline("Calculation was finished on iteration - $(val_iters) with error: $(val_len)")
    end
    if isnan(val_len) || val_iters >= val_maxiters
        throw(error("Policy function is NaN"))
    end

    if true#!value_from_file_flag#false#
        local_value = copy(value)
        local_aprime_nodes = copy(aprime_nodes)
        path = "$(@__DIR__)/val_aprime/"
        if Sys.iswindows()
            path = "$(@__DIR__)\\val_aprime\\"
        end

        @save "$(path)val_aprime_$(number_a_nodes)_$(number_u_nodes)_$(number_zeta_nodes)_$(number_alpha_m_nodes)_$(number_alpha_w_nodes)_fixed_occ.jld2" local_value local_aprime_nodes
    end

    return aprime_nodes, value
end

function find_stationary_distribution_pdf(fixed_occ_shares, a1_nodes,a_min,a_max,a_nodes, distr_tol, number_a_nodes, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, p_alpha, P_zeta, P_u, stat_P_u, P_alpha, number_u_nodes)

    number_asset_grid = number_a_nodes*10#250#500

    asset_grid = exp.(collect(range(log(a_min+1); stop=log(a_max+1), length=number_asset_grid))).-1
    #asset_grid = collect(range(a_min; stop=a_max, length=number_asset_grid))

    policy = Array{Any}(undef,3)
    a1_indices = Array{Any}(undef,3)
    lottery_prob_1 = Array{Any}(undef,3)
    lottery_prob_2 = Array{Any}(undef,3)
    for occ = 1:3
        policy[occ] = zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        a1_indices[occ] = Array{Int64}(undef,number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        lottery_prob_1[occ] = zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        lottery_prob_2[occ] = zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    end
    Threads.@threads for (occ,(u_i,(zeta_i,(alpha_m_i,alpha_w_i)))) in collect(Iterators.product(1:3,Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes)))))
        policy_function = Schumaker(a_nodes,a1_nodes[occ][:,u_i,zeta_i,alpha_m_i,alpha_w_i]; extrapolation = (Linear,Linear))
        for a_i in 1:number_asset_grid
            a1 = evaluate(policy_function, asset_grid[a_i])
            policy[occ][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = a1
            j_1 = sum(a1 .>= asset_grid)
            j = j_1 + 1

            a1_indices[occ][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = max(1,min(Int64(round(j_1;digits=0)),number_asset_grid))
            if j <= number_asset_grid && j_1 >= 1
                lottery_prob_1[occ][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (asset_grid[j] - a1             )/(asset_grid[j]-asset_grid[j_1])
                lottery_prob_2[occ][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = (a1            - asset_grid[j_1])/(asset_grid[j]-asset_grid[j_1])
            elseif j_1 == number_asset_grid
                lottery_prob_1[occ][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 1.0
                lottery_prob_2[occ][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 0.0
            end
        end
    end
    #throw(error)
    distr = Array{Any}(undef,3)
    for occ = 1:3
        distr[occ] = zeros(number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
    end

    # initialise stationary distribution
    print_sameline("Initialise stationary distribution")
    Threads.@threads for (occ,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:3,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
        distr[occ][floor(Int,number_asset_grid/5),:,zeta_i,alpha_m_i,alpha_w_i] = stat_P_u.*(P_zeta[zeta_i]*P_alpha[alpha_m_i,alpha_w_i])
    end

    oldK_supply = 0.0
    for occ = 1:3
        distr[occ] = (distr[occ]./sum(distr[occ])).*fixed_occ_shares[occ]
        oldK_supply += sum(distr[occ].*policy[occ])
    end

    # main loop
    distr_Delta   = 1.0#0.5  # update parameter for iteration process
    distr_maxiters= 1600
    distr_len     = Inf
    distr_iters   = 0
    if text_output
        print_sameline("Start of the main stationary distribution iteration loop")
    end

    distr_a1_z0 = copy(distr)
    new_distr = copy(distr)
    new_distr2 = copy(distr)

    while distr_len > distr_tol && distr_iters < distr_maxiters

        # calculate distribution
        # over future capital and current level of skills
        distr_a1_z0 .*= 0.0
        Threads.@threads for (occ,(u_i,(zeta_i,(alpha_m_i,alpha_w_i)))) in collect(Iterators.product(1:3,Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes)))))
            non_zero_asset_grid_iters = findall(!iszero,distr[occ][:,u_i,zeta_i,alpha_m_i,alpha_w_i])
            for a_i in non_zero_asset_grid_iters#1:number_asset_grid#
                j_1 = a1_indices[occ][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                distr_a1_z0[occ][j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += distr[occ][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[occ][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                if j_1 != number_asset_grid
                    distr_a1_z0[occ][j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += distr[occ][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[occ][a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                end
            end
        end
        #println_sameline(distr_a1_z0[:,1,1,1,1])
        #display(sum(distr_a1_z0))
        #throw(error)

        # first calculate transition
        # for people that do not change alpha_m and alpha_w
        new_distr .*= 0.0
        Threads.@threads for (occ,(alpha_m_i,alpha_w_i)) in collect(Iterators.product(1:3,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes)))
            # zeta is the transitory shock,
            # so add over all levels of zeta
            # and then draw new u_prime and new zeta_prime
            temp_distr_sum_zeta = sum(distr_a1_z0[occ][:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            #println_sameline(temp_distr_sum_zeta[:,1])
            #throw(error)
            temp_distr_sum_zeta = temp_distr_sum_zeta*P_u
            #println_sameline(temp_distr_sum_zeta[:,1])
            #throw(error)
            for zeta_prime_i in 1:number_zeta_nodes
                new_distr[occ][:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_distr_sum_zeta
            end
        end
        #println_sameline(new_distr[:,1,1,1,1])
        #throw(error)

        # second calculate transition
        # for people that change alpha_m and alpha_w
        #   and therefore draw new alpha_m_prime, alpha_w_prime
        #       and zeta_prime, u_prime
        #           from stationary distributions for this shocks

        # calculate sum of capital of all people who change skills
        distr_marginal_assets = Array{Any}(undef,3)
        for occ = 1:3
            distr_marginal_assets[occ] = sum(sum(sum(sum(distr_a1_z0[occ],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
        end

        new_distr2 .*= 0.0
        Threads.@threads for (occ,(u_prime_i,(zeta_prime_i,(alpha_m_prime_i,alpha_w_prime_i)))) in collect(Iterators.product(1:3,Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes)))))
            new_distr2[occ][:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] = (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*distr_marginal_assets[occ]
        end

        for occ = 1:3
            new_distr[occ] .+= new_distr2[occ]
        end

        #println_sameline(new_distr[:,1,1,1,1])
        #throw(error)

        #new_distr .= new_distr./sum(new_distr)

        distr_len = 0.0
        for occ = 1:3
            distr_len = max(distr_len, maximum(abs,distr[occ]-new_distr[occ]))
            distr[occ] .= new_distr[occ]
        end
        distr_iters += 1

        newK_supply = 0.0
        for occ = 1:3
            newK_supply += sum(distr[occ].*policy[occ])
        end
        K_s_error = abs(newK_supply-oldK_supply)
        if text_output
            print_sameline("Distr#$(distr_iters) - err: $(distr_len) - sumdistr: $(round.([sum(distr[1]),sum(distr[2]),sum(distr[3])];digits=2)), K_s:$(round(newK_supply;digits=6)), K_s_err:$(K_s_error)")
        end
        distr_len = K_s_error

        oldK_supply = newK_supply
        if fig_output
            #=
            cum_distr = copy(distr)
            for a_i in 2:number_a_nodes
                cum_distr[a_i,:,:,:,:] .+= cum_distr[a_i-1,:,:,:,:]
            end
            distr_label = permutedims(["$(u),$(zeta),$(alpha_m),$(alpha_w)" for u in 1:number_u_nodes for zeta in 1:number_zeta_nodes for alpha_m in 1:number_alpha_m_nodes for alpha_w in 1:number_alpha_w_nodes])
            plot_distr = [cum_distr[:,u,zeta,alpha_m,alpha_w] for u in 1:number_u_nodes for zeta in 1:number_zeta_nodes for alpha_m in 1:number_alpha_m_nodes for alpha_w in 1:number_alpha_w_nodes]
            display(plot(#=asset_grid,=#plot_distr,label=distr_label,#=xticks = a_min:5.0:a_max,=#xrotation=-45,legend=:outertopleft,foreground_color_legend = nothing))
            =#
        end
    end
    if text_output
        print_sameline("Calculation was finished on iteration - $(distr_iters) with error: $(distr_len)")
    end

    if fig_output
        #=
        distr_label = permutedims(["$(u),$(zeta),$(alpha_m),$(alpha_w)" for u in 1:number_u_nodes for zeta in 1:number_zeta_nodes for alpha_m in 1:number_alpha_m_nodes for alpha_w in 1:number_alpha_w_nodes])
        plot_distr = [distr[:,u,zeta,alpha_m,alpha_w] for u in 1:number_u_nodes for zeta in 1:number_zeta_nodes for alpha_m in 1:number_alpha_m_nodes for alpha_w in 1:number_alpha_w_nodes]
        display(plot(#=asset_grid,=#plot_distr,label=distr_label,#=xticks = a_min:5.0:a_max,=#xrotation=-45,legend=:outertopleft,foreground_color_legend = nothing))

        sum_distr = zeros(number_asset_grid)
        for a_i = 1:number_asset_grid
            for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
                sum_distr[a_i] += distr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
            end
        end
        display(plot(#=asset_grid,=#sum_distr,#=xticks = a_min:5.0:a_max,=#xrotation=-45,legend=:outertopleft,foreground_color_legend = nothing))
        #savefig("Figs/sum_distr.png")

        sum_cum_distr = copy(sum_distr)
        for a_i = 2:number_asset_grid
            sum_cum_distr[a_i] += sum_cum_distr[a_i-1]
        end
        display(plot(#=asset_grid,=#sum_cum_distr,#=xticks = a_min:5.0:a_max,=#xrotation=-45,legend=:outertopleft,foreground_color_legend = nothing))
        =#
    end

    if text_output
        print_sameline("Calculation of aggregate quantities for capital and labour")
    end

    return distr, number_asset_grid, asset_grid, policy, a1_indices, lottery_prob_1, lottery_prob_2
end

function find_aggregate_capital_labour_demand_supply_fixed_occ(number_asset_grid,asset_grid,policy,density_distr,r,w, labour_excess, labour_d, labour_s, capital_excess, capital_d)

    # find capital demand and supply
    K_demand = sum( [sum(density_distr[occ].*capital_d[occ]) for occ in 1:3])
    K_supply = sum( [sum(density_distr[occ].*policy[occ]) for occ in 1:3])
    # find excess capital demand on capital market
    Capital_excess = sum( [sum(density_distr[occ].*capital_excess[occ]) for occ in 1:3])/sum( [sum(density_distr[occ].*(capital_d[occ].+policy[occ])) for occ in 1:3])

    # find labour demand and supply
    L_demand = sum( [sum(density_distr[occ].*labour_d[occ]) for occ in 1:3])
    L_supply = sum( [sum(density_distr[occ].*labour_s[occ]) for occ in 1:3])
    # find excess labour demand on labour market
    Labor_excess = sum( [sum(density_distr[occ].*labour_excess[occ]) for occ in 1:3])/sum( [sum(density_distr[occ].*(labour_d[occ].+labour_s[occ])) for occ in 1:3])

    if text_output
        println_sameline()
        println_sameline("\nCurrent capital demand excess: $(Capital_excess)")
        println_sameline("Current aggregate capital demand: $(K_demand)")
        println_sameline("Current aggregate capital supply: $(K_supply)")

        println_sameline("\nCurrent labour demand excess: $(Labor_excess)")
        println_sameline("Current aggregate labour demand: $(L_demand)")
        println_sameline("Current aggregate labour supply: $(L_supply)")
    end

    return K_demand, K_supply, L_demand, L_supply, Capital_excess, Labor_excess
end

function quick_calculate_results(fixed_occ_shares, number_asset_grid,asset_grid,policy, a1_indices, lottery_prob_1, lottery_prob_2,density_distr,r,w, occ_choice, income, earnings, capital_d, credit, output, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, p_alpha, P_zeta, P_u, stat_P_u, P_alpha, number_u_nodes, cost_of_employing, z_m_nodes, z_w_nodes, lambda,delta,gamma,eta,theta,c_e, labour_d)

    number_non_zero_asset_grid = 0
    for (occ,(u_i,(zeta_i,(alpha_m_i,alpha_w_i)))) in collect(Iterators.product(1:3,Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes)))))
        temp = findlast(round.(density_distr[occ][:,u_i,zeta_i,alpha_m_i,alpha_w_i]; digits=7) .!= 0.0)
        if isnothing(temp)
            temp = 0
        end
        if temp > number_non_zero_asset_grid
            number_non_zero_asset_grid = temp
        end
    end

    agg_capital = sum( [sum(density_distr[occ].*capital_d[occ]) for occ in 1:3 ])
    agg_output = sum( [sum(density_distr[occ].*output[occ]) for occ in 1:3 ])
    capital_to_output = agg_capital/agg_output

    agg_credit = sum( [sum(density_distr[occ].*credit[occ]) for occ in 1:3 ])
    credit_to_output = agg_credit/agg_output

    share_of_workers = sum(density_distr[1])
    share_of_self_employed = sum(density_distr[2])
    share_of_employers = sum(density_distr[3])

    consumption = Array{Any}(undef,3)
    for occ = 1:3
        consumption[occ] = income[occ].-policy[occ]
    end
    avg_log_c = sum( [sum(density_distr[occ].*log.(max.(1e-12,consumption[occ]))) for occ in 1:3 ])
    var_log_c = sum( [sum(density_distr[occ].*(log.(max.(1e-12,consumption[occ])).- avg_log_c).^2) for occ in 1:3 ])

    #=
    avg_log_income = sum(density_distr.*log.(max.(1e-12,income)))
    var_log_income = sum(density_distr.*(log.(max.(1e-12,income)).- avg_log_income).^2)
    =#
    avg_log_income = sum( [sum(density_distr[occ].*log.(max.(1e-12,earnings[occ]))) for occ in 1:3 ])
    var_log_income = sum( [sum(density_distr[occ].*(log.(max.(1e-12,earnings[occ])).- avg_log_income).^2) for occ in 1:3 ])

    density_distr_w = density_distr[1]
    density_distr_se = density_distr[2]
    density_distr_emp = density_distr[3]

    density_distr_workers = density_distr_w./sum(density_distr_w)
    density_distr_workers_vec = vec(density_distr_workers)
    income_workers = income[1]
    income_workers_vec = vec(income_workers)

    index_w_non_zero = findall(x-> x>1e-5,density_distr_workers_vec)
    density_distr_workers_vec_non_zero = density_distr_workers_vec[index_w_non_zero]
    income_workers_vec_non_zero = income_workers_vec[index_w_non_zero]

    index_w = sortperm(income_workers_vec_non_zero)
    ddwvs = density_distr_workers_vec_non_zero[index_w]
    iwvs = income_workers_vec_non_zero[index_w]

    S_w_income = cumsum(ddwvs.*iwvs)./sum(ddwvs.*iwvs)
    gini_y_w = 1.0 - ( S_w_income[1]*ddwvs[1] + sum((S_w_income[2:end] .+ S_w_income[1:end-1]).*ddwvs[2:end]) )

    density_distr_entrepreneurs = [density_distr_se; density_distr_emp]./sum([density_distr_se; density_distr_emp])
    density_distr_entrepreneurs_vec = vec(density_distr_entrepreneurs)#reshape(density_distr_entrepreneurs, 1, :)
    income_entrepreneurs = [income[2]; income[3]]
    income_entrepreneurs_vec = vec(income_entrepreneurs)#reshape(income_entrepreneurs,1,:)

    index_ent_non_zero = findall(x-> x>1e-5,density_distr_entrepreneurs_vec)
    density_distr_entrepreneurs_vec_non_zero = density_distr_entrepreneurs_vec[index_ent_non_zero]
    income_entrepreneurs_vec_non_zero = income_entrepreneurs_vec[index_ent_non_zero]

    index_ent = sortperm(income_entrepreneurs_vec_non_zero)
    ddevs = density_distr_entrepreneurs_vec_non_zero[index_ent]
    ievs = income_entrepreneurs_vec_non_zero[index_ent]

    S_ent_income = cumsum(ddevs.*ievs)./sum(ddevs.*ievs)
    gini_y_ent = 1.0 - ( S_ent_income[1]*ddevs[1] + sum((S_ent_income[2:end] .+ S_ent_income[1:end-1]).*ddevs[2:end]) )

    occ_trans_matrix = zeros(3,3)
    occ_trans_matrix[1,1] = 1.0
    occ_trans_matrix[2,2] = 1.0
    occ_trans_matrix[3,3] = 1.0

    # productivity of managerial input, capital input, labour input
    agg_labour = sum( [sum(density_distr[occ].*labour_d[occ]) for occ in 1:3] )

    TFP_ideal = agg_output/(agg_capital^(eta) * agg_labour^(theta))
    TFP_data = agg_output/(agg_capital^(eta/(theta+eta)))

    Labour_productivity = agg_output/agg_labour
    Capital_productivity = agg_output/agg_capital

    agg_ent_earnings = sum( [sum(density_distr[occ].*earnings[occ]) for occ in 2:3] )
    Investment_productivity = (r+delta)*agg_capital/agg_ent_earnings

    #                    n-4        n-3       n-2                  n-1                   n
    additional_results = TFP_ideal, TFP_data, Labour_productivity, Capital_productivity, Investment_productivity

    if calc_add_results
        agg_cost_of_employing = sum(density_distr.*cost_of_employing)
        agg_cost_of_employing_as_share_of_output = agg_cost_of_employing/agg_output

        agg_credit_to_capital = agg_credit/agg_capital

        #Distribution of assets
        asset_distr = zeros(number_non_zero_asset_grid)
        asset_distr_w = zeros(number_non_zero_asset_grid)
        asset_distr_se = zeros(number_non_zero_asset_grid)
        asset_distr_emp = zeros(number_non_zero_asset_grid)
        for a_i in 1:number_non_zero_asset_grid
            asset_distr[a_i] = sum(density_distr[a_i,:,:,:,:])
            asset_distr_w[a_i] = sum(density_distr_w[a_i,:,:,:,:])
            asset_distr_se[a_i] = sum(density_distr_se[a_i,:,:,:,:])
            asset_distr_emp[a_i] = sum(density_distr_emp[a_i,:,:,:,:])
        end
        if fig_output
            main_title = "Distribution of assets"
            display(plot(asset_grid[1:number_non_zero_asset_grid],asset_distr,title=main_title,xrotation=-45,legend=:topright,foreground_color_legend = nothing))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)asset_distr_w_grid_lambda_$(lambda).png")
            display(plot(asset_distr,title=main_title,xrotation=-45,legend=:topright,foreground_color_legend = nothing))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)asset_distr_lambda_$(lambda).png")

            labels = ["W" "SP" "EMP"]
            display(plot(asset_grid[1:number_non_zero_asset_grid],[asset_distr_w,asset_distr_se,asset_distr_emp],label=labels,title=main_title,xrotation=-45,legend=:topright,foreground_color_legend = nothing))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)asset_distr_W_SP_EMP_w_grid_lambda_$(lambda).png")
            display(plot([asset_distr_w,asset_distr_se,asset_distr_emp],label=labels,title=main_title,xrotation=-45,legend=:topright,foreground_color_legend = nothing))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)asset_distr_W_SP_EMP_lambda_$(lambda).png")
        end

        # draw occupational choice
        agg_assets_is = zeros(4)
        agg_assets_is[1] = sum(sum(asset_grid[1:number_non_zero_asset_grid].*asset_distr) .>= asset_grid)
        agg_assets_is[2] = sum(sum(asset_grid[1:number_non_zero_asset_grid].*asset_distr_w)/sum(asset_distr_w) .>= asset_grid)
        agg_assets_is[3] = sum(sum(asset_grid[1:number_non_zero_asset_grid].*asset_distr_se)/sum(asset_distr_se) .>= asset_grid)
        agg_assets_is[4] = sum(sum(asset_grid[1:number_non_zero_asset_grid].*asset_distr_emp)/sum(asset_distr_emp) .>= asset_grid)
        agg_assets_is = Int64.(agg_assets_is)

        if fig_output
            for a_i in agg_assets_is
                draw_occupational_choice(z_m_nodes,z_w_nodes,occ_choice, a_i, asset_grid,r,w,lambda,delta,gamma,eta,theta,c_e)
            end
        end

        #Distribution of earnings
        number_earnings_grid = Int64(round(number_asset_grid/10;digits=0))#number_non_zero_asset_grid
        earnings_grid = exp.(collect(range(log(0.0+1.0); stop=log(maximum(earnings[1:number_non_zero_asset_grid,:,:,:,:])+1.0), length=number_earnings_grid))).-1.0
        earnings_distr = zeros(number_earnings_grid)
        earnings_distr_w = zeros(number_earnings_grid)
        earnings_distr_se = zeros(number_earnings_grid)
        earnings_distr_emp = zeros(number_earnings_grid)
        for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            #non_zero_asset_grid_iters = findall(!iszero,density_distr[:,u_i,zeta_i,alpha_m_i,alpha_w_i])
            for a_i in 1:number_non_zero_asset_grid#non_zero_asset_grid_iters#1:number_asset_grid#
                e = earnings[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                e_i_1 = sum(e .>= earnings_grid)
                e_i = e_i_1 + 1

                if e_i <= number_earnings_grid && e_i_1 >= 1
                    lp1 = (earnings_grid[e_i] - e             )/(earnings_grid[e_i]-earnings_grid[e_i_1])
                    lp2 = (e            - earnings_grid[e_i_1])/(earnings_grid[e_i]-earnings_grid[e_i_1])
                elseif e_i_1 == number_earnings_grid
                    lp1 = 1.0
                    lp2 = 0.0
                elseif e_i_1 < 1
                    e_i_1 = 1
                    e_i = e_i_1 + 1
                    lp1 = 1.0
                    lp2 = 0.0
                end

                earnings_distr[e_i_1] += density_distr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp1
                earnings_distr_w[e_i_1] += density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp1
                earnings_distr_se[e_i_1] += density_distr_se[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp1
                earnings_distr_emp[e_i_1] += density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp1
                if e_i_1 != number_earnings_grid
                    earnings_distr[e_i] += density_distr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp2
                    earnings_distr_w[e_i] += density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp2
                    earnings_distr_se[e_i] += density_distr_se[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp2
                    earnings_distr_emp[e_i] += density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp2
                end
            end
        end
        number_non_zero_earnings_grid = 0
        for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            temp = findlast(round.(earnings_distr; digits=7) .!= 0.0)
            if isnothing(temp)
                temp = 0
            end
            if temp > number_non_zero_earnings_grid
                number_non_zero_earnings_grid = temp
            end
        end
        if fig_output
            main_title = "Distribution of earnings"
            display(plot(earnings_grid[1:number_non_zero_earnings_grid],earnings_distr[1:number_non_zero_earnings_grid],title=main_title,xrotation=-45,legend=:topright,foreground_color_legend = nothing))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)earnings_distr_w_grid_lambda_$(lambda).png")
            display(plot(earnings_distr[1:number_non_zero_earnings_grid],title=main_title,xrotation=-45,legend=:topright,foreground_color_legend = nothing))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)earnings_distr_lambda_$(lambda).png")

            labels = ["W" "SP" "EMP"]
            display(plot(earnings_grid[1:number_non_zero_earnings_grid],[earnings_distr_w[1:number_non_zero_earnings_grid],earnings_distr_se[1:number_non_zero_earnings_grid],earnings_distr_emp[1:number_non_zero_earnings_grid]],label=labels,title=main_title,xrotation=-45,legend=:topright,foreground_color_legend = nothing))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)earnings_distr_W_SP_EMP_w_grid_lambda_$(lambda).png")
            display(plot([earnings_distr_w[1:number_non_zero_earnings_grid],earnings_distr_se[1:number_non_zero_earnings_grid],earnings_distr_emp[1:number_non_zero_earnings_grid]],label=labels,title=main_title,xrotation=-45,legend=:topright,foreground_color_legend = nothing))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)earnings_distr_W_SP_EMP_lambda_$(lambda).png")
        end

        #Distribution of income
        number_income_grid = Int64(round(number_asset_grid/10;digits=0))#number_non_zero_asset_grid
        income_grid = exp.(collect(range(log(0.0+1.0); stop=log(maximum(income[1:number_non_zero_asset_grid,:,:,:,:])+1.0), length=number_income_grid))).-1.0
        income_distr = zeros(number_income_grid)
        income_distr_w = zeros(number_income_grid)
        income_distr_se = zeros(number_income_grid)
        income_distr_emp = zeros(number_income_grid)
        for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            #non_zero_asset_grid_iters = findall(!iszero,density_distr[:,u_i,zeta_i,alpha_m_i,alpha_w_i])
            for a_i in 1:number_non_zero_asset_grid#non_zero_asset_grid_iters#1:number_asset_grid#
                e = income[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                e_i_1 = sum(e .>= income_grid)
                e_i = e_i_1 + 1

                if e_i <= number_income_grid && e_i_1 >= 1
                    lp1 = (income_grid[e_i] - e             )/(income_grid[e_i]-income_grid[e_i_1])
                    lp2 = (e            - income_grid[e_i_1])/(income_grid[e_i]-income_grid[e_i_1])
                elseif e_i_1 == number_income_grid
                    lp1 = 1.0
                    lp2 = 0.0
                elseif e_i_1 < 1
                    e_i_1 = 1
                    e_i = e_i_1 + 1
                    lp1 = 1.0
                    lp2 = 0.0
                end

                income_distr[e_i_1] += density_distr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp1
                income_distr_w[e_i_1] += density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp1
                income_distr_se[e_i_1] += density_distr_se[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp1
                income_distr_emp[e_i_1] += density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp1
                if e_i_1 != number_income_grid
                    income_distr[e_i] += density_distr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp2
                    income_distr_w[e_i] += density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp2
                    income_distr_se[e_i] += density_distr_se[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp2
                    income_distr_emp[e_i] += density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp2
                end
            end
        end
        number_non_zero_income_grid = 0
        for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            temp = findlast(round.(income_distr; digits=7) .!= 0.0)
            if isnothing(temp)
                temp = 0
            end
            if temp > number_non_zero_income_grid
                number_non_zero_income_grid = temp
            end
        end
        if fig_output
            main_title = "Distribution of income"
            display(plot(income_grid[1:number_non_zero_income_grid],income_distr[1:number_non_zero_income_grid],title=main_title,xrotation=-45,legend=:topright,foreground_color_legend = nothing))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)income_distr_w_grid_lambda_$(lambda).png")
            display(plot(income_distr[1:number_non_zero_income_grid],title=main_title,xrotation=-45,legend=:topright,foreground_color_legend = nothing))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)income_distr_lambda_$(lambda).png")

            labels = ["W" "SP" "EMP"]
            display(plot(income_grid[1:number_non_zero_income_grid],[income_distr_w[1:number_non_zero_income_grid],income_distr_se[1:number_non_zero_income_grid],income_distr_emp[1:number_non_zero_income_grid]],label=labels,title=main_title,xrotation=-45,legend=:topright,foreground_color_legend = nothing))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)income_distr_W_SP_EMP_w_grid_lambda_$(lambda).png")
            display(plot([income_distr_w[1:number_non_zero_income_grid],income_distr_se[1:number_non_zero_income_grid],income_distr_emp[1:number_non_zero_income_grid]],label=labels,title=main_title,xrotation=-45,legend=:topright,foreground_color_legend = nothing))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)income_distr_W_SP_EMP_lambda_$(lambda).png")
        end

        #Distribution of consumption
        number_consumption_grid = Int64(round(number_asset_grid/10;digits=0))#number_non_zero_asset_grid
        consumption_grid = exp.(collect(range(log(0.0+1.0); stop=log(maximum(consumption[1:number_non_zero_asset_grid,:,:,:,:])+1.0), length=number_consumption_grid))).-1.0
        consumption_distr = zeros(number_consumption_grid)
        consumption_distr_w = zeros(number_consumption_grid)
        consumption_distr_se = zeros(number_consumption_grid)
        consumption_distr_emp = zeros(number_consumption_grid)
        for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            #non_zero_asset_grid_iters = findall(!iszero,density_distr[:,u_i,zeta_i,alpha_m_i,alpha_w_i])
            for a_i in 1:number_non_zero_asset_grid#non_zero_asset_grid_iters#1:number_asset_grid#
                e = consumption[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                e_i_1 = sum(e .>= consumption_grid)
                e_i = e_i_1 + 1

                if e_i <= number_consumption_grid && e_i_1 >= 1
                    lp1 = (consumption_grid[e_i] - e             )/(consumption_grid[e_i]-consumption_grid[e_i_1])
                    lp2 = (e            - consumption_grid[e_i_1])/(consumption_grid[e_i]-consumption_grid[e_i_1])
                elseif e_i_1 == number_consumption_grid
                    lp1 = 1.0
                    lp2 = 0.0
                elseif e_i_1 < 1
                    e_i_1 = 1
                    e_i = e_i_1 + 1
                    lp1 = 1.0
                    lp2 = 0.0
                end

                consumption_distr[e_i_1] += density_distr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp1
                consumption_distr_w[e_i_1] += density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp1
                consumption_distr_se[e_i_1] += density_distr_se[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp1
                consumption_distr_emp[e_i_1] += density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp1
                if e_i_1 != number_consumption_grid
                    consumption_distr[e_i] += density_distr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp2
                    consumption_distr_w[e_i] += density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp2
                    consumption_distr_se[e_i] += density_distr_se[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp2
                    consumption_distr_emp[e_i] += density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lp2
                end
            end
        end
        number_non_zero_consumption_grid = 0
        for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            temp = findlast(round.(consumption_distr; digits=7) .!= 0.0)
            if isnothing(temp)
                temp = 0
            end
            if temp > number_non_zero_consumption_grid
                number_non_zero_consumption_grid = temp
            end
        end
        if fig_output
            main_title = "Distribution of consumption"
            display(plot(consumption_grid[1:number_non_zero_consumption_grid],consumption_distr[1:number_non_zero_consumption_grid],title=main_title,xrotation=-45,legend=:topright,foreground_color_legend = nothing))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)consumption_distr_w_grid_lambda_$(lambda).png")
            display(plot(consumption_distr[1:number_non_zero_consumption_grid],title=main_title,xrotation=-45,legend=:topright,foreground_color_legend = nothing))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)consumption_distr_lambda_$(lambda).png")

            labels = ["W" "SP" "EMP"]
            display(plot(consumption_grid[1:number_non_zero_consumption_grid],[consumption_distr_w[1:number_non_zero_consumption_grid],consumption_distr_se[1:number_non_zero_consumption_grid],consumption_distr_emp[1:number_non_zero_consumption_grid]],label=labels,title=main_title,xrotation=-45,legend=:topright,foreground_color_legend = nothing))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)consumption_distr_W_SP_EMP_w_grid_lambda_$(lambda).png")
            display(plot([consumption_distr_w[1:number_non_zero_consumption_grid],consumption_distr_se[1:number_non_zero_consumption_grid],consumption_distr_emp[1:number_non_zero_consumption_grid]],label=labels,title=main_title,xrotation=-45,legend=:topright,foreground_color_legend = nothing))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)consumption_distr_W_SP_EMP_lambda_$(lambda).png")
        end

        # draw policy
        if fig_output
            labels = permutedims(["$(u),$(zeta),$(alpha_m),$(alpha_w)" for u in 1:number_u_nodes for zeta in 1:number_zeta_nodes for alpha_m in 1:number_alpha_m_nodes for alpha_w in 1:number_alpha_w_nodes])
            p2 = plot(asset_grid[1:number_non_zero_asset_grid],[policy[1:number_non_zero_asset_grid,u,zeta,alpha_m,alpha_w] for u in 1:number_u_nodes for zeta in 1:number_zeta_nodes for alpha_m in 1:number_alpha_m_nodes for alpha_w in 1:number_alpha_w_nodes], label=labels, legend=:outertopleft)
            display(plot(p2, title="Policy functions (all) at lambda=$(lambda)"))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)policy_lambda_$(lambda).png")

            # high vs low fixed effect (managerial and working)
            labels = permutedims(["$(u),$(zeta),$(alpha_m),$(alpha_w)" for u in [Int64(round((1+number_u_nodes)/2;digits=0))] for zeta in [Int64(round((1+number_zeta_nodes)/2;digits=0))] for alpha_m in 1:number_alpha_m_nodes for alpha_w in 1:number_alpha_w_nodes])
            p2 = plot(asset_grid[1:number_non_zero_asset_grid],[policy[1:number_non_zero_asset_grid,u,zeta,alpha_m,alpha_w] for u in [Int64(round((1+number_u_nodes)/2;digits=0))] for zeta in [Int64(round((1+number_zeta_nodes)/2;digits=0))] for alpha_m in 1:number_alpha_m_nodes for alpha_w in 1:number_alpha_w_nodes], label=labels, legend=:outertopleft)
            display(plot(p2, title="Policy functions (fixed effects) at lambda=$(lambda)"))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)policy_fixed_lambda_$(lambda).png")

            labels = permutedims(["$(u),$(zeta),$(alpha_m),$(alpha_w)" for u in [Int64(round((1+number_u_nodes)/2;digits=0))] for zeta in [Int64(round((1+number_zeta_nodes)/2;digits=0))] for alpha_m in [1,number_alpha_m_nodes] for alpha_w in [1,number_alpha_w_nodes]])
            p2 = plot(asset_grid[1:number_non_zero_asset_grid],[policy[1:number_non_zero_asset_grid,u,zeta,alpha_m,alpha_w] for u in [Int64(round((1+number_u_nodes)/2;digits=0))] for zeta in [Int64(round((1+number_zeta_nodes)/2;digits=0))] for alpha_m in [1,number_alpha_m_nodes] for alpha_w in [1,number_alpha_w_nodes]], label=labels, legend=:outertopleft)
            display(plot(p2, title="Policy functions (high vs low fixed effects) at lambda=$(lambda)"))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)policy_hvsl_fixed_lambda_$(lambda).png")

            # high vs low persistent (managerial and working)
            labels = permutedims(["$(u),$(zeta),$(alpha_m),$(alpha_w)" for u in [1, Int64(sqrt(number_u_nodes)), Int64(round((1+number_u_nodes)/2;digits=0)), number_u_nodes+1-Int64(sqrt(number_u_nodes)), number_u_nodes] for zeta in [Int64(round((1+number_zeta_nodes)/2;digits=0))] for alpha_m in [1,number_alpha_m_nodes] for alpha_w in [Int64(round((1+number_alpha_w_nodes)/2;digits=0))]])
            p2 = plot(asset_grid[1:number_non_zero_asset_grid],[policy[1:number_non_zero_asset_grid,u,zeta,alpha_m,alpha_w] for u in [1, Int64(sqrt(number_u_nodes)), Int64(round((1+number_u_nodes)/2;digits=0)), number_u_nodes+1-Int64(sqrt(number_u_nodes)), number_u_nodes] for zeta in [Int64(round((1+number_zeta_nodes)/2;digits=0))] for alpha_m in [1,number_alpha_m_nodes] for alpha_w in [Int64(round((1+number_alpha_w_nodes)/2;digits=0))]], label=labels, legend=:outertopleft)
            display(plot(p2, title="Policy functions (high vs low persistent) at lambda=$(lambda)"))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)policy_hvsl_persistent_lambda_$(lambda).png")

            # high vs low transititory
            labels = permutedims(["$(u),$(zeta),$(alpha_m),$(alpha_w)" for u in [Int64(round((1+number_u_nodes)/2;digits=0))] for zeta in 1:number_zeta_nodes for alpha_m in [1,number_alpha_m_nodes] for alpha_w in [Int64(round((1+number_alpha_w_nodes)/2;digits=0))]])
            p2 = plot(asset_grid[1:number_non_zero_asset_grid],[policy[1:number_non_zero_asset_grid,u,zeta,alpha_m,alpha_w] for u in [Int64(round((1+number_u_nodes)/2;digits=0))] for zeta in 1:number_zeta_nodes for alpha_m in [1,number_alpha_m_nodes] for alpha_w in [Int64(round((1+number_alpha_w_nodes)/2;digits=0))]], label=labels, legend=:outertopleft)
            display(plot(p2, title="Policy functions (high vs low transitory) at lambda=$(lambda)"))
            path = "/Users/maxsolodarenko/OneDrive - University of Glasgow/PhD/Year 2/AllubErosa Transitional dynamics/Code/Test_firm_r_w/Results/Figures/"
            if Sys.iswindows()
                path = "C:\\Users\\2288144s\\OneDrive - University of Glasgow\\PhD\\Year 2\\AllubErosa Transitional dynamics\\Code\\Test_firm_r_w\\Results\\Figures\\"
            end
            savefig("$(path)policy_hvsl_transitory_lambda_$(lambda).png")
        end

        # mobility of:
        #               occupations
        #               assets in quintiles between
        #                               and within occupations
        #               income in quintiles between
        #                               and within occupations
        #               earnings in quintiles between
        #                                 and within occupations
        #               consumption in quintiles between
        #                                    and within occupations
        density_distr_pr = copy(density_distr)
        density_distr_w_pr = copy(density_distr_w)
        density_distr_se_pr = copy(density_distr_se)
        density_distr_emp_pr = copy(density_distr_emp)

        number_of_quintiles = 5
        quintiles = collect(range(0.0; stop=1.0, length=number_of_quintiles+1))[2:end-1]

        # CALCULATE 0 PERIOD ASSETS, EARNINGS IN QUINTILES
        # calculate levels of assets for each quintile
        density_distr_marginal_assets = sum(sum(sum(sum(density_distr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
        cum_density_marginal_assets = cumsum(density_distr_marginal_assets)
        inverse_cum_density_marginal_assets = Schumaker(cum_density_marginal_assets,asset_grid; extrapolation = (Linear,Linear))
        asset_quintiles = zeros(number_of_quintiles-1)

        # calculate levels of earnings for each quintile
        cum_density_marginal_earnings = cumsum(earnings_distr)
        inverse_cum_density_marginal_earnings = Schumaker(cum_density_marginal_earnings,earnings_grid; extrapolation = (Linear,Linear))
        earnings_quintiles = zeros(number_of_quintiles-1)

        # calculate levels of income for each quintile
        cum_density_marginal_income = cumsum(income_distr)
        #display(plot(income_grid[1:number_non_zero_income_grid],income_distr[1:number_non_zero_income_grid]))
        #display(plot(income_grid[1:number_non_zero_income_grid],cum_density_marginal_income[1:number_non_zero_income_grid]))
        #throw(error)
        inverse_cum_density_marginal_income = Schumaker(cum_density_marginal_income,income_grid; extrapolation = (Linear,Linear))
        income_quintiles = zeros(number_of_quintiles-1)

        # calculate levels of consumption for each quintile
        cum_density_marginal_consumption = cumsum(consumption_distr)
        inverse_cum_density_marginal_consumption = Schumaker(cum_density_marginal_consumption,consumption_grid; extrapolation = (Linear,Linear))
        consumption_quintiles = zeros(number_of_quintiles-1)

        Threads.@threads for q in 1:number_of_quintiles-1
            asset_quintiles[q] = evaluate(inverse_cum_density_marginal_assets, quintiles[q])
            earnings_quintiles[q] = evaluate(inverse_cum_density_marginal_earnings, quintiles[q])
            income_quintiles[q] = evaluate(inverse_cum_density_marginal_income, quintiles[q])
            consumption_quintiles[q] = evaluate(inverse_cum_density_marginal_consumption, quintiles[q])
        end
        #display(income_quintiles)
        #throw(error)
        # create variable that indicates quintile each (a_i,u_i,zeta_i,alpha_m_i,alpha_w_i) have
        asset_quintile_choice = copy(density_distr).*0.0
        earnings_quintile_choice = copy(density_distr).*0.0
        income_quintile_choice = copy(density_distr).*0.0
        consumption_quintile_choice = copy(density_distr).*0.0
        Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            for a_i in 1:number_asset_grid
                for aq in 1:number_of_quintiles-1
                    if asset_grid[a_i] <= asset_quintiles[aq] && asset_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == 0.0
                        asset_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = aq
                    end
                    if earnings[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] <= earnings_quintiles[aq] && earnings_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == 0.0
                        earnings_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = aq
                    end
                    if income[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] <= income_quintiles[aq] && income_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == 0.0
                        income_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = aq
                    end
                    if consumption[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] <= consumption_quintiles[aq] && consumption_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == 0.0
                        consumption_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = aq
                    end
                end
                if asset_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == 0.0
                    asset_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = number_of_quintiles
                end
                if earnings_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == 0.0
                    earnings_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = number_of_quintiles
                end
                if income_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == 0.0
                    income_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = number_of_quintiles
                end
                if consumption_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == 0.0
                    consumption_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = number_of_quintiles
                end
            end
        end
        asset_quintile_choice_pr = copy(asset_quintile_choice)
        earnings_quintile_choice_pr = copy(earnings_quintile_choice)
        income_quintile_choice_pr = copy(income_quintile_choice)
        consumption_quintile_choice_pr = copy(consumption_quintile_choice)

        denisty_distr_assets_quintiles = zeros(number_of_quintiles, number_asset_grid, number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        denisty_distr_assets_quintiles_w = copy(denisty_distr_assets_quintiles)
        denisty_distr_assets_quintiles_se = copy(denisty_distr_assets_quintiles)
        denisty_distr_assets_quintiles_emp = copy(denisty_distr_assets_quintiles)

        denisty_distr_earnings_quintiles = copy(denisty_distr_assets_quintiles)
        denisty_distr_earnings_quintiles_w = copy(denisty_distr_assets_quintiles)
        denisty_distr_earnings_quintiles_se = copy(denisty_distr_assets_quintiles)
        denisty_distr_earnings_quintiles_emp = copy(denisty_distr_assets_quintiles)

        denisty_distr_income_quintiles = copy(denisty_distr_assets_quintiles)
        denisty_distr_income_quintiles_w = copy(denisty_distr_assets_quintiles)
        denisty_distr_income_quintiles_se = copy(denisty_distr_assets_quintiles)
        denisty_distr_income_quintiles_emp = copy(denisty_distr_assets_quintiles)

        denisty_distr_consumption_quintiles = copy(denisty_distr_assets_quintiles)
        denisty_distr_consumption_quintiles_w = copy(denisty_distr_assets_quintiles)
        denisty_distr_consumption_quintiles_se = copy(denisty_distr_assets_quintiles)
        denisty_distr_consumption_quintiles_emp = copy(denisty_distr_assets_quintiles)

        for aq = 1:number_of_quintiles
            Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
                for a_i in 1:number_asset_grid
                    if asset_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == aq
                        denisty_distr_assets_quintiles[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = density_distr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    end
                    if earnings_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == aq
                        denisty_distr_earnings_quintiles[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = density_distr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    end
                    if income_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == aq
                        denisty_distr_income_quintiles[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = density_distr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    end
                    if consumption_quintile_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == aq
                        denisty_distr_consumption_quintiles[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = density_distr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    end
                end
            end
            denisty_distr_assets_quintiles_w[  aq,:,:,:,:,:] = w_choice  .*denisty_distr_assets_quintiles[aq,:,:,:,:,:]
            denisty_distr_assets_quintiles_se[ aq,:,:,:,:,:] = se_choice .*denisty_distr_assets_quintiles[aq,:,:,:,:,:]
            denisty_distr_assets_quintiles_emp[aq,:,:,:,:,:] = emp_choice.*denisty_distr_assets_quintiles[aq,:,:,:,:,:]

            denisty_distr_earnings_quintiles_w[  aq,:,:,:,:,:] = w_choice  .*denisty_distr_earnings_quintiles[aq,:,:,:,:,:]
            denisty_distr_earnings_quintiles_se[ aq,:,:,:,:,:] = se_choice .*denisty_distr_earnings_quintiles[aq,:,:,:,:,:]
            denisty_distr_earnings_quintiles_emp[aq,:,:,:,:,:] = emp_choice.*denisty_distr_earnings_quintiles[aq,:,:,:,:,:]

            denisty_distr_income_quintiles_w[  aq,:,:,:,:,:] = w_choice  .*denisty_distr_income_quintiles[aq,:,:,:,:,:]
            denisty_distr_income_quintiles_se[ aq,:,:,:,:,:] = se_choice .*denisty_distr_income_quintiles[aq,:,:,:,:,:]
            denisty_distr_income_quintiles_emp[aq,:,:,:,:,:] = emp_choice.*denisty_distr_income_quintiles[aq,:,:,:,:,:]

            denisty_distr_consumption_quintiles_w[  aq,:,:,:,:,:] = w_choice  .*denisty_distr_consumption_quintiles[aq,:,:,:,:,:]
            denisty_distr_consumption_quintiles_se[ aq,:,:,:,:,:] = se_choice .*denisty_distr_consumption_quintiles[aq,:,:,:,:,:]
            denisty_distr_consumption_quintiles_emp[aq,:,:,:,:,:] = emp_choice.*denisty_distr_income_quintiles[aq,:,:,:,:,:]
        end

        denisty_distr_assets_quintiles_pr = copy(denisty_distr_assets_quintiles)
        denisty_distr_assets_quintiles_w_pr = copy(denisty_distr_assets_quintiles_w)
        denisty_distr_assets_quintiles_se_pr = copy(denisty_distr_assets_quintiles_se)
        denisty_distr_assets_quintiles_emp_pr = copy(denisty_distr_assets_quintiles_emp)

        denisty_distr_earnings_quintiles_pr = copy(denisty_distr_earnings_quintiles)
        denisty_distr_earnings_quintiles_w_pr = copy(denisty_distr_earnings_quintiles_w)
        denisty_distr_earnings_quintiles_se_pr = copy(denisty_distr_earnings_quintiles_se)
        denisty_distr_earnings_quintiles_emp_pr = copy(denisty_distr_earnings_quintiles_emp)

        denisty_distr_income_quintiles_pr = copy(denisty_distr_income_quintiles)
        denisty_distr_income_quintiles_w_pr = copy(denisty_distr_income_quintiles_w)
        denisty_distr_income_quintiles_se_pr = copy(denisty_distr_income_quintiles_se)
        denisty_distr_income_quintiles_emp_pr = copy(denisty_distr_income_quintiles_emp)

        denisty_distr_consumption_quintiles_pr = copy(denisty_distr_consumption_quintiles)
        denisty_distr_consumption_quintiles_w_pr = copy(denisty_distr_consumption_quintiles_w)
        denisty_distr_consumption_quintiles_se_pr = copy(denisty_distr_consumption_quintiles_se)
        denisty_distr_consumption_quintiles_emp_pr = copy(denisty_distr_consumption_quintiles_emp)

        TT = 100
        occ_mobility_matrices = zeros(TT, 3,3)

        asset_quintile_mobility_matrices = zeros(TT, number_of_quintiles,number_of_quintiles)
        asset_quintile_w_mobility_matrices = zeros(TT, number_of_quintiles,number_of_quintiles+1)
        asset_quintile_se_mobility_matrices = zeros(TT, number_of_quintiles,number_of_quintiles+1)
        asset_quintile_emp_mobility_matrices = zeros(TT, number_of_quintiles,number_of_quintiles+1)

        earnings_quintile_mobility_matrices = zeros(TT, number_of_quintiles,number_of_quintiles)
        earnings_quintile_w_mobility_matrices = zeros(TT, number_of_quintiles,number_of_quintiles+1)
        earnings_quintile_se_mobility_matrices = zeros(TT, number_of_quintiles,number_of_quintiles+1)
        earnings_quintile_emp_mobility_matrices = zeros(TT, number_of_quintiles,number_of_quintiles+1)

        income_quintile_mobility_matrices = zeros(TT, number_of_quintiles,number_of_quintiles)
        income_quintile_w_mobility_matrices = zeros(TT, number_of_quintiles,number_of_quintiles+1)
        income_quintile_se_mobility_matrices = zeros(TT, number_of_quintiles,number_of_quintiles+1)
        income_quintile_emp_mobility_matrices = zeros(TT, number_of_quintiles,number_of_quintiles+1)

        consumption_quintile_mobility_matrices = zeros(TT, number_of_quintiles,number_of_quintiles)
        consumption_quintile_w_mobility_matrices = zeros(TT, number_of_quintiles,number_of_quintiles+1)
        consumption_quintile_se_mobility_matrices = zeros(TT, number_of_quintiles,number_of_quintiles+1)
        consumption_quintile_emp_mobility_matrices = zeros(TT, number_of_quintiles,number_of_quintiles+1)

        for tt = 1:TT
            print_sameline("Mobility - $(tt)/$(TT)")

            temp_density_distr = copy(density_distr_pr)
            density_distr_pr .*= 0.0
            temp_density_distr_w = copy(density_distr_w_pr)
            density_distr_w_pr .*= 0.0
            temp_density_distr_se = copy(density_distr_se_pr)
            density_distr_se_pr .*= 0.0
            temp_density_distr_emp = copy(density_distr_emp_pr)
            density_distr_emp_pr .*= 0.0

            temp_denisty_distr_assets_quintiles = copy(denisty_distr_assets_quintiles_pr)
            denisty_distr_assets_quintiles_pr .*= 0.0
            temp_denisty_distr_assets_quintiles_w = copy(denisty_distr_assets_quintiles_w_pr)
            denisty_distr_assets_quintiles_w_pr .*= 0.0
            temp_denisty_distr_assets_quintiles_se = copy(denisty_distr_assets_quintiles_se_pr)
            denisty_distr_assets_quintiles_se_pr .*= 0.0
            temp_denisty_distr_assets_quintiles_emp = copy(denisty_distr_assets_quintiles_emp_pr)
            denisty_distr_assets_quintiles_emp_pr .*= 0.0

            temp_denisty_distr_earnings_quintiles = copy(denisty_distr_earnings_quintiles_pr)
            denisty_distr_earnings_quintiles_pr .*= 0.0
            temp_denisty_distr_earnings_quintiles_w = copy(denisty_distr_earnings_quintiles_w_pr)
            denisty_distr_earnings_quintiles_w_pr .*= 0.0
            temp_denisty_distr_earnings_quintiles_se = copy(denisty_distr_earnings_quintiles_se_pr)
            denisty_distr_earnings_quintiles_se_pr .*= 0.0
            temp_denisty_distr_earnings_quintiles_emp = copy(denisty_distr_earnings_quintiles_emp_pr)
            denisty_distr_earnings_quintiles_emp_pr .*= 0.0

            temp_denisty_distr_income_quintiles = copy(denisty_distr_income_quintiles_pr)
            denisty_distr_income_quintiles_pr .*= 0.0
            temp_denisty_distr_income_quintiles_w = copy(denisty_distr_income_quintiles_w_pr)
            denisty_distr_income_quintiles_w_pr .*= 0.0
            temp_denisty_distr_income_quintiles_se = copy(denisty_distr_income_quintiles_se_pr)
            denisty_distr_income_quintiles_se_pr .*= 0.0
            temp_denisty_distr_income_quintiles_emp = copy(denisty_distr_income_quintiles_emp_pr)
            denisty_distr_income_quintiles_emp_pr .*= 0.0

            temp_denisty_distr_consumption_quintiles = copy(denisty_distr_consumption_quintiles_pr)
            denisty_distr_consumption_quintiles_pr .*= 0.0
            temp_denisty_distr_consumption_quintiles_w = copy(denisty_distr_consumption_quintiles_w_pr)
            denisty_distr_consumption_quintiles_w_pr .*= 0.0
            temp_denisty_distr_consumption_quintiles_se = copy(denisty_distr_consumption_quintiles_se_pr)
            denisty_distr_consumption_quintiles_se_pr .*= 0.0
            temp_denisty_distr_consumption_quintiles_emp = copy(denisty_distr_consumption_quintiles_emp_pr)
            denisty_distr_consumption_quintiles_emp_pr .*= 0.0

            Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
                for a_i in 1:number_asset_grid
                    j_1 = a1_indices[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    #println_sameline([j_1,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i])
                    density_distr_w_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    density_distr_se_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_density_distr_se[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    density_distr_emp_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    for aq = 1:number_of_quintiles
                        denisty_distr_assets_quintiles_pr[aq, j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_assets_quintiles[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        denisty_distr_assets_quintiles_w_pr[aq, j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_assets_quintiles_w[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        denisty_distr_assets_quintiles_se_pr[aq, j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_assets_quintiles_se[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        denisty_distr_assets_quintiles_emp_pr[aq, j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_assets_quintiles_emp[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                        denisty_distr_earnings_quintiles_pr[aq, j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_earnings_quintiles[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        denisty_distr_earnings_quintiles_w_pr[aq, j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_earnings_quintiles_w[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        denisty_distr_earnings_quintiles_se_pr[aq, j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_earnings_quintiles_se[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        denisty_distr_earnings_quintiles_emp_pr[aq, j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_earnings_quintiles_emp[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                        denisty_distr_income_quintiles_pr[aq, j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_income_quintiles[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        denisty_distr_income_quintiles_w_pr[aq, j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_income_quintiles_w[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        denisty_distr_income_quintiles_se_pr[aq, j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_income_quintiles_se[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        denisty_distr_income_quintiles_emp_pr[aq, j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_income_quintiles_emp[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                        denisty_distr_consumption_quintiles_pr[aq, j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_consumption_quintiles[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        denisty_distr_consumption_quintiles_w_pr[aq, j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_consumption_quintiles_w[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        denisty_distr_consumption_quintiles_se_pr[aq, j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_consumption_quintiles_se[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        denisty_distr_consumption_quintiles_emp_pr[aq, j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_consumption_quintiles_emp[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    end
                    if j_1 != number_asset_grid
                        density_distr_w_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        density_distr_se_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_density_distr_se[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        density_distr_emp_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        for aq = 1:number_of_quintiles
                            denisty_distr_assets_quintiles_pr[aq, j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_assets_quintiles[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                            denisty_distr_assets_quintiles_w_pr[aq, j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_assets_quintiles_w[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                            denisty_distr_assets_quintiles_se_pr[aq, j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_assets_quintiles_se[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                            denisty_distr_assets_quintiles_emp_pr[aq, j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_assets_quintiles_emp[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                            denisty_distr_earnings_quintiles_pr[aq, j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_earnings_quintiles[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                            denisty_distr_earnings_quintiles_w_pr[aq, j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_earnings_quintiles_w[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                            denisty_distr_earnings_quintiles_se_pr[aq, j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_earnings_quintiles_se[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                            denisty_distr_earnings_quintiles_emp_pr[aq, j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_earnings_quintiles_emp[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                            denisty_distr_income_quintiles_pr[aq, j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_income_quintiles[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                            denisty_distr_income_quintiles_w_pr[aq, j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_income_quintiles_w[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                            denisty_distr_income_quintiles_se_pr[aq, j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_income_quintiles_se[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                            denisty_distr_income_quintiles_emp_pr[aq, j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_income_quintiles_emp[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                            denisty_distr_consumption_quintiles_pr[aq, j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_consumption_quintiles[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                            denisty_distr_consumption_quintiles_w_pr[aq, j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_consumption_quintiles_w[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                            denisty_distr_consumption_quintiles_se_pr[aq, j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_consumption_quintiles_se[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                            denisty_distr_consumption_quintiles_emp_pr[aq, j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += temp_denisty_distr_consumption_quintiles_emp[aq, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                        end
                    end
                end
            end

            density_distr_marginal_assets_w = sum(sum(sum(sum(density_distr_w_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
            density_distr_marginal_assets_se = sum(sum(sum(sum(density_distr_se_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
            density_distr_marginal_assets_emp = sum(sum(sum(sum(density_distr_emp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
            density_distr_marginal_assets = density_distr_marginal_assets_w .+ density_distr_marginal_assets_se .+ density_distr_marginal_assets_emp

            density_distr_marginal_assets_quintiles = zeros(number_of_quintiles, number_asset_grid)
            density_distr_marginal_assets_quintiles_w = copy(density_distr_marginal_assets_quintiles)
            density_distr_marginal_assets_quintiles_se = copy(density_distr_marginal_assets_quintiles)
            density_distr_marginal_assets_quintiles_emp = copy(density_distr_marginal_assets_quintiles)

            density_distr_marginal_earnings_quintiles = copy(density_distr_marginal_assets_quintiles)
            density_distr_marginal_earnings_quintiles_w = copy(density_distr_marginal_assets_quintiles)
            density_distr_marginal_earnings_quintiles_se = copy(density_distr_marginal_assets_quintiles)
            density_distr_marginal_earnings_quintiles_emp = copy(density_distr_marginal_assets_quintiles)

            density_distr_marginal_income_quintiles = copy(density_distr_marginal_assets_quintiles)
            density_distr_marginal_income_quintiles_w = copy(density_distr_marginal_assets_quintiles)
            density_distr_marginal_income_quintiles_se = copy(density_distr_marginal_assets_quintiles)
            density_distr_marginal_income_quintiles_emp = copy(density_distr_marginal_assets_quintiles)

            density_distr_marginal_consumption_quintiles = copy(density_distr_marginal_assets_quintiles)
            density_distr_marginal_consumption_quintiles_w = copy(density_distr_marginal_assets_quintiles)
            density_distr_marginal_consumption_quintiles_se = copy(density_distr_marginal_assets_quintiles)
            density_distr_marginal_consumption_quintiles_emp = copy(density_distr_marginal_assets_quintiles)
            Threads.@threads for aq = 1:number_of_quintiles
                density_distr_marginal_assets_quintiles[aq,:] = sum(sum(sum(sum(denisty_distr_assets_quintiles_pr[aq,:,:,:,:,:],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
                density_distr_marginal_assets_quintiles_w[aq,:] = sum(sum(sum(sum(denisty_distr_assets_quintiles_w_pr[aq,:,:,:,:,:],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
                density_distr_marginal_assets_quintiles_se[aq,:] = sum(sum(sum(sum(denisty_distr_assets_quintiles_se_pr[aq,:,:,:,:,:],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
                density_distr_marginal_assets_quintiles_emp[aq,:] = sum(sum(sum(sum(denisty_distr_assets_quintiles_emp_pr[aq,:,:,:,:,:],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]

                density_distr_marginal_earnings_quintiles[aq,:] = sum(sum(sum(sum(denisty_distr_earnings_quintiles_pr[aq,:,:,:,:,:],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
                density_distr_marginal_earnings_quintiles_w[aq,:] = sum(sum(sum(sum(denisty_distr_earnings_quintiles_w_pr[aq,:,:,:,:,:],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
                density_distr_marginal_earnings_quintiles_se[aq,:] = sum(sum(sum(sum(denisty_distr_earnings_quintiles_se_pr[aq,:,:,:,:,:],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
                density_distr_marginal_earnings_quintiles_emp[aq,:] = sum(sum(sum(sum(denisty_distr_earnings_quintiles_emp_pr[aq,:,:,:,:,:],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]

                density_distr_marginal_income_quintiles[aq,:] = sum(sum(sum(sum(denisty_distr_income_quintiles_pr[aq,:,:,:,:,:],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
                density_distr_marginal_income_quintiles_w[aq,:] = sum(sum(sum(sum(denisty_distr_income_quintiles_w_pr[aq,:,:,:,:,:],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
                density_distr_marginal_income_quintiles_se[aq,:] = sum(sum(sum(sum(denisty_distr_income_quintiles_se_pr[aq,:,:,:,:,:],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
                density_distr_marginal_income_quintiles_emp[aq,:] = sum(sum(sum(sum(denisty_distr_income_quintiles_emp_pr[aq,:,:,:,:,:],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]

                density_distr_marginal_consumption_quintiles[aq,:] = sum(sum(sum(sum(denisty_distr_consumption_quintiles_pr[aq,:,:,:,:,:],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
                density_distr_marginal_consumption_quintiles_w[aq,:] = sum(sum(sum(sum(denisty_distr_consumption_quintiles_w_pr[aq,:,:,:,:,:],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
                density_distr_marginal_consumption_quintiles_se[aq,:] = sum(sum(sum(sum(denisty_distr_consumption_quintiles_se_pr[aq,:,:,:,:,:],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
                density_distr_marginal_consumption_quintiles_emp[aq,:] = sum(sum(sum(sum(denisty_distr_consumption_quintiles_emp_pr[aq,:,:,:,:,:],dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
            end

            Threads.@threads for (alpha_m_i,alpha_w_i) in collect(Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))
                temp_w = sum(density_distr_w_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
                temp_w = temp_w*P_u
                temp_se = sum(density_distr_se_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
                temp_se = temp_se*P_u
                temp_emp = sum(density_distr_emp_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
                temp_emp = temp_emp*P_u

                temp_asset_quintiles = zeros(number_of_quintiles, number_asset_grid, number_u_nodes)
                temp_asset_quintiles_w = copy(temp_asset_quintiles)
                temp_asset_quintiles_se = copy(temp_asset_quintiles)
                temp_asset_quintiles_emp = copy(temp_asset_quintiles)

                temp_earnings_quintiles = copy(temp_asset_quintiles)
                temp_earnings_quintiles_w = copy(temp_earnings_quintiles)
                temp_earnings_quintiles_se = copy(temp_earnings_quintiles)
                temp_earnings_quintiles_emp = copy(temp_earnings_quintiles)

                temp_income_quintiles = copy(temp_asset_quintiles)
                temp_income_quintiles_w = copy(temp_income_quintiles)
                temp_income_quintiles_se = copy(temp_income_quintiles)
                temp_income_quintiles_emp = copy(temp_income_quintiles)

                temp_consumption_quintiles = copy(temp_asset_quintiles)
                temp_consumption_quintiles_w = copy(temp_consumption_quintiles)
                temp_consumption_quintiles_se = copy(temp_consumption_quintiles)
                temp_consumption_quintiles_emp = copy(temp_consumption_quintiles)
                for aq = 1:number_of_quintiles
                    temp_asset_quintiles[aq, :,:] = sum(denisty_distr_assets_quintiles_pr[aq,:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
                    temp_asset_quintiles_w[aq, :,:] = sum(denisty_distr_assets_quintiles_w_pr[aq,:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
                    temp_asset_quintiles_se[aq, :,:] = sum(denisty_distr_assets_quintiles_se_pr[aq,:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
                    temp_asset_quintiles_emp[aq, :,:] = sum(denisty_distr_assets_quintiles_emp_pr[aq,:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]

                    temp_earnings_quintiles[aq, :,:] = sum(denisty_distr_earnings_quintiles_pr[aq,:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
                    temp_earnings_quintiles_w[aq, :,:] = sum(denisty_distr_earnings_quintiles_w_pr[aq,:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
                    temp_earnings_quintiles_se[aq, :,:] = sum(denisty_distr_earnings_quintiles_se_pr[aq,:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
                    temp_earnings_quintiles_emp[aq, :,:] = sum(denisty_distr_earnings_quintiles_emp_pr[aq,:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]

                    temp_income_quintiles[aq, :,:] = sum(denisty_distr_income_quintiles_pr[aq,:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
                    temp_income_quintiles_w[aq, :,:] = sum(denisty_distr_income_quintiles_w_pr[aq,:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
                    temp_income_quintiles_se[aq, :,:] = sum(denisty_distr_income_quintiles_se_pr[aq,:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
                    temp_income_quintiles_emp[aq, :,:] = sum(denisty_distr_income_quintiles_emp_pr[aq,:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]

                    temp_consumption_quintiles[aq, :,:] = sum(denisty_distr_consumption_quintiles_pr[aq,:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
                    temp_consumption_quintiles_w[aq, :,:] = sum(denisty_distr_consumption_quintiles_w_pr[aq,:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
                    temp_consumption_quintiles_se[aq, :,:] = sum(denisty_distr_consumption_quintiles_se_pr[aq,:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
                    temp_consumption_quintiles_emp[aq, :,:] = sum(denisty_distr_consumption_quintiles_emp_pr[aq,:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
                end
                for zeta_prime_i in 1:number_zeta_nodes
                    density_distr_w_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_w
                    density_distr_se_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_se
                    density_distr_emp_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_emp
                    for aq = 1:number_of_quintiles
                        denisty_distr_assets_quintiles_pr[aq, :,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_asset_quintiles[aq, :,:]
                        denisty_distr_assets_quintiles_w_pr[  aq, :,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_asset_quintiles_w[aq, :,:]
                        denisty_distr_assets_quintiles_se_pr[ aq, :,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_asset_quintiles_se[aq, :,:]
                        denisty_distr_assets_quintiles_emp_pr[aq, :,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_asset_quintiles_emp[aq, :,:]

                        denisty_distr_earnings_quintiles_pr[aq, :,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_earnings_quintiles[aq, :,:]
                        denisty_distr_earnings_quintiles_w_pr[  aq, :,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_earnings_quintiles_w[aq, :,:]
                        denisty_distr_earnings_quintiles_se_pr[ aq, :,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_earnings_quintiles_se[aq, :,:]
                        denisty_distr_earnings_quintiles_emp_pr[aq, :,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_earnings_quintiles_emp[aq, :,:]

                        denisty_distr_income_quintiles_pr[aq, :,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_income_quintiles[aq, :,:]
                        denisty_distr_income_quintiles_w_pr[  aq, :,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_income_quintiles_w[aq, :,:]
                        denisty_distr_income_quintiles_se_pr[ aq, :,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_income_quintiles_se[aq, :,:]
                        denisty_distr_income_quintiles_emp_pr[aq, :,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_income_quintiles_emp[aq, :,:]

                        denisty_distr_consumption_quintiles_pr[aq, :,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_consumption_quintiles[aq, :,:]
                        denisty_distr_consumption_quintiles_w_pr[  aq, :,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_consumption_quintiles_w[aq, :,:]
                        denisty_distr_consumption_quintiles_se_pr[ aq, :,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_consumption_quintiles_se[aq, :,:]
                        denisty_distr_consumption_quintiles_emp_pr[aq, :,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_consumption_quintiles_emp[aq, :,:]
                    end
                end
            end
            Threads.@threads for (u_prime_i,(zeta_prime_i,(alpha_m_prime_i,alpha_w_prime_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
                density_distr_w_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_assets_w
                density_distr_se_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_assets_se
                density_distr_emp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_assets_emp
                for aq = 1:number_of_quintiles
                    denisty_distr_assets_quintiles_pr[aq, :,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_assets_quintiles[aq,:]
                    denisty_distr_assets_quintiles_w_pr[  aq, :,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_assets_quintiles_w[aq,:]
                    denisty_distr_assets_quintiles_se_pr[ aq, :,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_assets_quintiles_se[aq,:]
                    denisty_distr_assets_quintiles_emp_pr[aq, :,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_assets_quintiles_emp[aq,:]

                    denisty_distr_earnings_quintiles_pr[aq, :,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_earnings_quintiles[aq,:]
                    denisty_distr_earnings_quintiles_w_pr[  aq, :,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_earnings_quintiles_w[aq,:]
                    denisty_distr_earnings_quintiles_se_pr[ aq, :,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_earnings_quintiles_se[aq,:]
                    denisty_distr_earnings_quintiles_emp_pr[aq, :,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_earnings_quintiles_emp[aq,:]

                    denisty_distr_income_quintiles_pr[aq, :,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_income_quintiles[aq,:]
                    denisty_distr_income_quintiles_w_pr[  aq, :,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_income_quintiles_w[aq,:]
                    denisty_distr_income_quintiles_se_pr[ aq, :,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_income_quintiles_se[aq,:]
                    denisty_distr_income_quintiles_emp_pr[aq, :,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_income_quintiles_emp[aq,:]

                    denisty_distr_consumption_quintiles_pr[aq, :,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_consumption_quintiles[aq,:]
                    denisty_distr_consumption_quintiles_w_pr[  aq, :,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_consumption_quintiles_w[aq,:]
                    denisty_distr_consumption_quintiles_se_pr[ aq, :,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_consumption_quintiles_se[aq,:]
                    denisty_distr_consumption_quintiles_emp_pr[aq, :,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_consumption_quintiles_emp[aq,:]
                end
            end
            density_distr_pr = density_distr_w_pr .+ density_distr_se_pr .+ density_distr_emp_pr

            occ_mobility_matrices[tt,1,1] = sum(w_choice.*density_distr_w_pr)/max(1e-6,sum(density_distr_w))
            occ_mobility_matrices[tt,1,2] = sum(se_choice.*density_distr_w_pr)/max(1e-6,sum(density_distr_w))
            occ_mobility_matrices[tt,1,3] = sum(emp_choice.*density_distr_w_pr)/max(1e-6,sum(density_distr_w))
            occ_mobility_matrices[tt,2,1] = sum(w_choice.*density_distr_se_pr)/max(1e-6,sum(density_distr_se))
            occ_mobility_matrices[tt,2,2] = sum(se_choice.*density_distr_se_pr)/max(1e-6,sum(density_distr_se))
            occ_mobility_matrices[tt,2,3] = sum(emp_choice.*density_distr_se_pr)/max(1e-6,sum(density_distr_se))
            occ_mobility_matrices[tt,3,1] = sum(w_choice.*density_distr_emp_pr)/max(1e-6,sum(density_distr_emp))
            occ_mobility_matrices[tt,3,2] = sum(se_choice.*density_distr_emp_pr)/max(1e-6,sum(density_distr_emp))
            occ_mobility_matrices[tt,3,3] = sum(emp_choice.*density_distr_emp_pr)/max(1e-6,sum(density_distr_emp))

            # create variable that indicates quintile each (a_i,u_i,zeta_i,alpha_m_i,alpha_w_i) have
            asset_quintile_choice_pr = copy(density_distr).*0.0
            earnings_quintile_choice_pr = copy(density_distr).*0.0
            income_quintile_choice_pr = copy(density_distr).*0.0
            consumption_quintile_choice_pr = copy(density_distr).*0.0
            Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
                for a_i in 1:number_asset_grid
                    for aq in 1:number_of_quintiles-1
                        if asset_grid[a_i] <= asset_quintiles[aq] && asset_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == 0.0
                            asset_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = aq
                        end
                        if earnings[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] <= earnings_quintiles[aq] && earnings_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == 0.0
                            earnings_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = aq
                        end
                        if income[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] <= income_quintiles[aq] && income_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == 0.0
                            income_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = aq
                        end
                        if consumption[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] <= consumption_quintiles[aq] && consumption_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == 0.0
                            consumption_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = aq
                        end
                    end
                    if asset_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == 0.0
                        asset_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = number_of_quintiles
                    end
                    if earnings_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == 0.0
                        earnings_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = number_of_quintiles
                    end
                    if income_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == 0.0
                        income_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = number_of_quintiles
                    end
                    if consumption_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == 0.0
                        consumption_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = number_of_quintiles
                    end
                end
            end
            # calculate transition probability between quintiles
            for aq1 = 1:number_of_quintiles
                for aq2 = 1:number_of_quintiles
                    temp_distr2 = copy(density_distr).*0.0
                    temp_distr2_w = copy(density_distr).*0.0
                    temp_distr2_se = copy(density_distr).*0.0
                    temp_distr2_emp = copy(density_distr).*0.0

                    temp_distr2_earnings = copy(density_distr).*0.0
                    temp_distr2_earnings_w = copy(density_distr).*0.0
                    temp_distr2_earnings_se = copy(density_distr).*0.0
                    temp_distr2_earnings_emp = copy(density_distr).*0.0

                    temp_distr2_income = copy(density_distr).*0.0
                    temp_distr2_income_w = copy(density_distr).*0.0
                    temp_distr2_income_se = copy(density_distr).*0.0
                    temp_distr2_income_emp = copy(density_distr).*0.0

                    temp_distr2_consumption = copy(density_distr).*0.0
                    temp_distr2_consumption_w = copy(density_distr).*0.0
                    temp_distr2_consumption_se = copy(density_distr).*0.0
                    temp_distr2_consumption_emp = copy(density_distr).*0.0

                    Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
                        for a_i in 1:number_asset_grid
                            if asset_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == aq2
                                temp_distr2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = denisty_distr_assets_quintiles_pr[aq1,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                                temp_distr2_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = denisty_distr_assets_quintiles_w_pr[aq1,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                                temp_distr2_se[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = denisty_distr_assets_quintiles_se_pr[aq1,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                                temp_distr2_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = denisty_distr_assets_quintiles_emp_pr[aq1,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                            end
                            if earnings_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == aq2
                                temp_distr2_earnings[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = denisty_distr_earnings_quintiles_pr[aq1,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                                temp_distr2_earnings_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = denisty_distr_earnings_quintiles_w_pr[aq1,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                                temp_distr2_earnings_se[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = denisty_distr_earnings_quintiles_se_pr[aq1,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                                temp_distr2_earnings_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = denisty_distr_earnings_quintiles_emp_pr[aq1,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                            end
                            if income_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == aq2
                                temp_distr2_income[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = denisty_distr_income_quintiles_pr[aq1,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                                temp_distr2_income_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = denisty_distr_income_quintiles_w_pr[aq1,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                                temp_distr2_income_se[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = denisty_distr_income_quintiles_se_pr[aq1,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                                temp_distr2_income_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = denisty_distr_income_quintiles_emp_pr[aq1,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                            end
                            if consumption_quintile_choice_pr[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] == aq2
                                temp_distr2_consumption[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = denisty_distr_consumption_quintiles_pr[aq1,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                                temp_distr2_consumption_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = denisty_distr_consumption_quintiles_w_pr[aq1,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                                temp_distr2_consumption_se[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = denisty_distr_consumption_quintiles_se_pr[aq1,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                                temp_distr2_consumption_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = denisty_distr_consumption_quintiles_emp_pr[aq1,a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                            end
                        end
                    end
                    # probability of households with aq1 assets to move to aq2
                    asset_quintile_mobility_matrices[tt,aq1,aq2] = sum(temp_distr2)/max(1e-6,sum(denisty_distr_assets_quintiles[aq1,:,:,:,:,:]))
                    asset_quintile_w_mobility_matrices[tt,aq1,aq2] = sum(w_choice.*temp_distr2_w)/max(1e-6,sum(denisty_distr_assets_quintiles_w[aq1,:,:,:,:,:]))
                    asset_quintile_se_mobility_matrices[tt,aq1,aq2] = sum(se_choice.*temp_distr2_se)/max(1e-6,sum(denisty_distr_assets_quintiles_se[aq1,:,:,:,:,:]))
                    asset_quintile_emp_mobility_matrices[tt,aq1,aq2] = sum(emp_choice.*temp_distr2_emp)/max(1e-6,sum(denisty_distr_assets_quintiles_emp[aq1,:,:,:,:,:]))

                    earnings_quintile_mobility_matrices[tt,aq1,aq2] = sum(temp_distr2_earnings)/max(1e-6,sum(denisty_distr_earnings_quintiles[aq1,:,:,:,:,:]))
                    earnings_quintile_w_mobility_matrices[tt,aq1,aq2] = sum(w_choice.*temp_distr2_earnings_w)/max(1e-6,sum(denisty_distr_earnings_quintiles_w[aq1,:,:,:,:,:]))
                    earnings_quintile_se_mobility_matrices[tt,aq1,aq2] = sum(se_choice.*temp_distr2_earnings_se)/max(1e-6,sum(denisty_distr_earnings_quintiles_se[aq1,:,:,:,:,:]))
                    earnings_quintile_emp_mobility_matrices[tt,aq1,aq2] = sum(emp_choice.*temp_distr2_earnings_emp)/max(1e-6,sum(denisty_distr_earnings_quintiles_emp[aq1,:,:,:,:,:]))

                    income_quintile_mobility_matrices[tt,aq1,aq2] = sum(temp_distr2_income)/max(1e-6,sum(denisty_distr_income_quintiles[aq1,:,:,:,:,:]))
                    income_quintile_w_mobility_matrices[tt,aq1,aq2] = sum(w_choice.*temp_distr2_income_w)/max(1e-6,sum(denisty_distr_income_quintiles_w[aq1,:,:,:,:,:]))
                    income_quintile_se_mobility_matrices[tt,aq1,aq2] = sum(se_choice.*temp_distr2_income_se)/max(1e-6,sum(denisty_distr_income_quintiles_se[aq1,:,:,:,:,:]))
                    income_quintile_emp_mobility_matrices[tt,aq1,aq2] = sum(emp_choice.*temp_distr2_income_emp)/max(1e-6,sum(denisty_distr_income_quintiles_emp[aq1,:,:,:,:,:]))

                    consumption_quintile_mobility_matrices[tt,aq1,aq2] = sum(temp_distr2_consumption)/max(1e-6,sum(denisty_distr_consumption_quintiles[aq1,:,:,:,:,:]))
                    consumption_quintile_w_mobility_matrices[tt,aq1,aq2] = sum(w_choice.*temp_distr2_consumption_w)/max(1e-6,sum(denisty_distr_consumption_quintiles_w[aq1,:,:,:,:,:]))
                    consumption_quintile_se_mobility_matrices[tt,aq1,aq2] = sum(se_choice.*temp_distr2_consumption_se)/max(1e-6,sum(denisty_distr_consumption_quintiles_se[aq1,:,:,:,:,:]))
                    consumption_quintile_emp_mobility_matrices[tt,aq1,aq2] = sum(emp_choice.*temp_distr2_consumption_emp)/max(1e-6,sum(denisty_distr_consumption_quintiles_emp[aq1,:,:,:,:,:]))

                end
                # probability of households with aq1 assets and occupation to move to another occupation
                asset_quintile_w_mobility_matrices[tt,aq1,number_of_quintiles+1] = 1.0-sum(asset_quintile_w_mobility_matrices[tt,aq1,1:number_of_quintiles])
                asset_quintile_se_mobility_matrices[tt,aq1,number_of_quintiles+1] = 1.0-sum(asset_quintile_se_mobility_matrices[tt,aq1,1:number_of_quintiles])
                asset_quintile_emp_mobility_matrices[tt,aq1,number_of_quintiles+1] = 1.0-sum(asset_quintile_emp_mobility_matrices[tt,aq1,1:number_of_quintiles])

                earnings_quintile_w_mobility_matrices[tt,aq1,number_of_quintiles+1] = 1.0-sum(earnings_quintile_w_mobility_matrices[tt,aq1,1:number_of_quintiles])
                earnings_quintile_se_mobility_matrices[tt,aq1,number_of_quintiles+1] = 1.0-sum(earnings_quintile_se_mobility_matrices[tt,aq1,1:number_of_quintiles])
                earnings_quintile_emp_mobility_matrices[tt,aq1,number_of_quintiles+1] = 1.0-sum(earnings_quintile_emp_mobility_matrices[tt,aq1,1:number_of_quintiles])

                income_quintile_w_mobility_matrices[tt,aq1,number_of_quintiles+1] = 1.0-sum(income_quintile_w_mobility_matrices[tt,aq1,1:number_of_quintiles])
                income_quintile_se_mobility_matrices[tt,aq1,number_of_quintiles+1] = 1.0-sum(income_quintile_se_mobility_matrices[tt,aq1,1:number_of_quintiles])
                income_quintile_emp_mobility_matrices[tt,aq1,number_of_quintiles+1] = 1.0-sum(income_quintile_emp_mobility_matrices[tt,aq1,1:number_of_quintiles])

                consumption_quintile_w_mobility_matrices[tt,aq1,number_of_quintiles+1] = 1.0-sum(consumption_quintile_w_mobility_matrices[tt,aq1,1:number_of_quintiles])
                consumption_quintile_se_mobility_matrices[tt,aq1,number_of_quintiles+1] = 1.0-sum(consumption_quintile_se_mobility_matrices[tt,aq1,1:number_of_quintiles])
                consumption_quintile_emp_mobility_matrices[tt,aq1,number_of_quintiles+1] = 1.0-sum(consumption_quintile_emp_mobility_matrices[tt,aq1,1:number_of_quintiles])
            end
            #display(occ_mobility_matrices[tt,:,:])

            #display(asset_quintile_mobility_matrices[tt,:,:])

            #display(asset_quintile_w_mobility_matrices[tt,:,:])
            #display(sum(asset_quintile_w_mobility_matrices[tt,:,:]))
            #throw(error)
        end

        if text_output
            println_sameline("Occupation mobility matrix after $(TT) periods")
            display(occ_mobility_matrices[TT,:,:])

            println_sameline("Asset quintile mobility matrix after $(TT) periods")
            display(asset_quintile_mobility_matrices[TT,:,:])
            println_sameline("Asset quintile for workers mobility matrix after $(TT) periods")
            display(asset_quintile_w_mobility_matrices[TT,:,:])
            println_sameline("Asset quintile for sole proprietors mobility matrix after $(TT) periods")
            display(asset_quintile_se_mobility_matrices[TT,:,:])
            println_sameline("Asset quintile for employers mobility matrix after $(TT) periods")
            display(asset_quintile_emp_mobility_matrices[TT,:,:])

            println_sameline("Earnings quintile mobility matrix after $(TT) periods")
            display(earnings_quintile_mobility_matrices[TT,:,:])
            println_sameline("Earnings quintile for workers mobility matrix after $(TT) periods")
            display(earnings_quintile_w_mobility_matrices[TT,:,:])
            println_sameline("Earnings quintile for sole proprietors mobility matrix after $(TT) periods")
            display(earnings_quintile_se_mobility_matrices[TT,:,:])
            println_sameline("Earnings quintile for employers mobility matrix after $(TT) periods")
            display(earnings_quintile_emp_mobility_matrices[TT,:,:])

            println_sameline("Income quintile mobility matrix after $(TT) periods")
            display(income_quintile_mobility_matrices[TT,:,:])
            println_sameline("Income quintile for workers mobility matrix after $(TT) periods")
            display(income_quintile_w_mobility_matrices[TT,:,:])
            println_sameline("Income quintile for sole proprietors mobility matrix after $(TT) periods")
            display(income_quintile_se_mobility_matrices[TT,:,:])
            println_sameline("Income quintile for employers mobility matrix after $(TT) periods")
            display(income_quintile_emp_mobility_matrices[TT,:,:])

            println_sameline("Consumption quintile mobility matrix after $(TT) periods")
            display(consumption_quintile_mobility_matrices[TT,:,:])
            println_sameline("Consumption quintile for workers mobility matrix after $(TT) periods")
            display(consumption_quintile_w_mobility_matrices[TT,:,:])
            println_sameline("Consumption quintile for sole proprietors mobility matrix after $(TT) periods")
            display(consumption_quintile_se_mobility_matrices[TT,:,:])
            println_sameline("Consumption quintile for employers mobility matrix after $(TT) periods")
            display(consumption_quintile_emp_mobility_matrices[TT,:,:])
        end

        open("$(@__DIR__)\\Results\\Mobility\\occ_mobility_matrices.txt", "a") do f
            for tt=1:TT
                write(f, "0-$(tt):\n")
                write(f, "$(occ_mobility_matrices[tt,:,:])\n\n")
            end
        end

        open("$(@__DIR__)\\Results\\Mobility\\asset_quintile_mobility_matrices.txt", "a") do f
            for tt=1:TT
                write(f, "0-$(tt):\n")
                write(f, "$(asset_quintile_mobility_matrices[tt,:,:])\n")
                write(f, "0-$(tt) - w:\n")
                write(f, "$(asset_quintile_w_mobility_matrices[tt,:,:])\n")
                write(f, "0-$(tt) - sp:\n")
                write(f, "$(asset_quintile_se_mobility_matrices[tt,:,:])\n")
                write(f, "0-$(tt) - emp:\n")
                write(f, "$(asset_quintile_emp_mobility_matrices[tt,:,:])\n\n")
            end
        end

        open("$(@__DIR__)\\Results\\Mobility\\earnings_quintile_mobility_matrices.txt", "a") do f
            for tt=1:TT
                write(f, "0-$(tt):\n")
                write(f, "$(earnings_quintile_mobility_matrices[tt,:,:])\n")
                write(f, "0-$(tt) - w:\n")
                write(f, "$(earnings_quintile_w_mobility_matrices[tt,:,:])\n")
                write(f, "0-$(tt) - sp:\n")
                write(f, "$(earnings_quintile_se_mobility_matrices[tt,:,:])\n")
                write(f, "0-$(tt) - emp:\n")
                write(f, "$(earnings_quintile_emp_mobility_matrices[tt,:,:])\n\n")
            end
        end

        open("$(@__DIR__)\\Results\\Mobility\\income_quintile_mobility_matrices.txt", "a") do f
            for tt=1:TT
                write(f, "0-$(tt):\n")
                write(f, "$(income_quintile_mobility_matrices[tt,:,:])\n")
                write(f, "0-$(tt) - w:\n")
                write(f, "$(income_quintile_w_mobility_matrices[tt,:,:])\n")
                write(f, "0-$(tt) - sp:\n")
                write(f, "$(income_quintile_se_mobility_matrices[tt,:,:])\n")
                write(f, "0-$(tt) - emp:\n")
                write(f, "$(income_quintile_emp_mobility_matrices[tt,:,:])\n\n")
            end
        end

        open("$(@__DIR__)\\Results\\Mobility\\consumption_quintile_mobility_matrices.txt", "a") do f
            for tt=1:TT
                write(f, "0-$(tt):\n")
                write(f, "$(consumption_quintile_mobility_matrices[tt,:,:])\n")
                write(f, "0-$(tt) - w:\n")
                write(f, "$(consumption_quintile_w_mobility_matrices[tt,:,:])\n")
                write(f, "0-$(tt) - sp:\n")
                write(f, "$(consumption_quintile_se_mobility_matrices[tt,:,:])\n")
                write(f, "0-$(tt) - emp:\n")
                write(f, "$(consumption_quintile_emp_mobility_matrices[tt,:,:])\n\n")
            end
        end

        occ_trans_matrix = copy(occ_mobility_matrices[1,:,:])

        #                        1                      2                                         3                           4           5             6              7                8                    9             10             11               12                13                 14                             15                 16          17           18             19              20               21
        #22                     23               24                25                  26                   27                    28                                29                    30                               31                                 32                                  33                                    34                                  35                                    36                                     37
        #38                               39                                  40                                   41                                     42                                     43                                       44                                        45
        add_additional_results = agg_cost_of_employing, agg_cost_of_employing_as_share_of_output, number_non_zero_asset_grid, asset_distr,asset_distr_w,asset_distr_se,asset_distr_emp, number_earnings_grid,earnings_grid,earnings_distr,earnings_distr_w,earnings_distr_se,earnings_distr_emp,number_non_zero_earnings_grid, number_income_grid,income_grid,income_distr,income_distr_w,income_distr_se,income_distr_emp,number_non_zero_income_grid, number_consumption_grid,consumption_grid,consumption_distr,consumption_distr_w,consumption_distr_se,consumption_distr_emp,number_non_zero_consumption_grid, occ_mobility_matrices,asset_quintile_mobility_matrices,asset_quintile_w_mobility_matrices,asset_quintile_se_mobility_matrices,asset_quintile_emp_mobility_matrices, earnings_quintile_mobility_matrices,earnings_quintile_w_mobility_matrices,earnings_quintile_se_mobility_matrices,earnings_quintile_emp_mobility_matrices, income_quintile_mobility_matrices,income_quintile_w_mobility_matrices,income_quintile_se_mobility_matrices,income_quintile_emp_mobility_matrices, consumption_quintile_mobility_matrices,consumption_quintile_w_mobility_matrices,consumption_quintile_se_mobility_matrices,consumption_quintile_emp_mobility_matrices

        additional_results = add_additional_results, additional_results
    else

        #additional_results = []
    end

    if text_output
        println_sameline()
        println_sameline("\nCapital-to-output: $(capital_to_output)")
        println_sameline("Credit-to-output: $(credit_to_output)")

        println_sameline("Share of workers: $(share_of_workers)")
        println_sameline("Share of self-employed: $(share_of_self_employed)")
        println_sameline("Share of employers: $(share_of_employers)")

        println_sameline("Variance of log-consumption: $(var_log_c)")

        println_sameline("Occupational transition probabilities:")
        display(occ_trans_matrix)

        println_sameline("Variance of log-earnings: $(var_log_income)")

        println_sameline("Gini for worker's income: $(gini_y_w)")
        println_sameline("Gini for entrepreneurial's income: $(gini_y_ent)")
        if calc_add_results
            println_sameline("\n------Additional results---------")
            println_sameline("Aggregate cost of employing: $(agg_cost_of_employing) ($(agg_cost_of_employing_as_share_of_output*100)%)")
            println_sameline("Aggregate Loan-to-Value (Credit-to-Capital): $(agg_credit_to_capital)")
        end
        println_sameline("\n------Productivity results---------")
        println_sameline("TFP_ideal: $(TFP_ideal*100)%")
        println_sameline("TFP_data: $(TFP_data*100)%")
        println_sameline("Labour productivity: $(Labour_productivity)")
        println_sameline("Capital productivity (Output-to-capital): $(Capital_productivity)")
        println_sameline("Capital productivity (Expenses-to-revenue = (r+delta)*Capital-earnings): $(Investment_productivity)")

    end

    #return capital_to_output, credit_to_output, avg_y_w, var_y_w, avg_y_ent, var_y_ent, avg_y_emp, avg_y_se, share_of_workers, share_of_self_employed, share_of_employers, occ_trans_matrix, avg_log_c, var_log_c
    return capital_to_output, credit_to_output, share_of_workers, share_of_self_employed, share_of_employers, var_log_c, occ_trans_matrix, var_log_income, gini_y_w, gini_y_ent, agg_capital, agg_output, agg_credit, avg_log_income, avg_log_c, additional_results
end

function AllubErosa(r,w, gps, gaps, ps, ao, fixed_occ_shares)
    distr_tol = gps[3]
    val_tol = gps[4]

    number_a_nodes = gaps[1]
    number_zeta_nodes = gaps[4]
    number_alpha_m_nodes = gaps[5]
    number_alpha_w_nodes = gaps[6]

    crra = ps[20]
    beta = ps[2]
    p_alpha = ps[13]

    z_m_nodes = ao[1]
    z_w_nodes = ao[2]
    P_zeta = ao[3]
    P_u = ao[4]
    stat_P_u = ao[5]
    P_alpha = ao[6]
    number_u_nodes = ao[7]

    lambda = ps[1]
    delta = ps[3]
    gamma = ps[4]
    eta = ps[5]
    theta = ps[6]
    c_e = ps[7]

    z_m_nodes = ao[1]
    z_w_nodes = ao[2]
    number_u_nodes = ao[7]

    ############################
    a_min = 0.01#0.0
    a_max = 2*entrep_optimal(maximum(z_m_nodes),maximum(z_w_nodes),r,w, delta, gamma, eta, theta, c_e)[1]
    #a_max = 2500.0#37288869442.4668#400000.0#500.0#50.0

    a_nodes = exp.(collect(range(log(a_min+1); stop=log(a_max+1), length=number_a_nodes))).-1
    #a_nodes = collect(range(a_min; stop=a_max, length=number_a_nodes))
    ###########################
    # calculate income profiles for different households
    income, earnings = compute_income_and_earnings_fixed_occ(a_nodes,number_a_nodes,r, w, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, lambda, delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)
    # VFI that finds policy function for future capital
    a1_nodes, value = find_policy_fixed_occ(a_min,a_max,a_nodes,r,w, income,earnings, val_tol, number_a_nodes, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, beta, p_alpha, P_zeta, P_u, stat_P_u, P_alpha, number_u_nodes, crra)

    ##############Compute K_supply

    ########## Compute the stationary distribution of assets
    density_distr, number_asset_grid, asset_grid, policy, a1_indices, lottery_prob_1, lottery_prob_2 = find_stationary_distribution_pdf(fixed_occ_shares, a1_nodes,a_min,a_max,a_nodes, distr_tol, number_a_nodes, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, p_alpha, P_zeta, P_u, stat_P_u, P_alpha, number_u_nodes)
    # calculate income profiles for different households
    occ_choice, income, earnings, capital_excess, capital_d, credit, labour_excess, labour_d, labour_s, deposit, output, cost_of_employing, managerial_input = compute_income_profile_fixed_occ(asset_grid,number_asset_grid,r,w, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, lambda, delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)

    K_demand, K_supply, L_demand, L_supply, Capital_excess, Labor_excess = find_aggregate_capital_labour_demand_supply_fixed_occ(number_asset_grid,asset_grid,policy,density_distr,r,w, labour_excess, labour_d, labour_s, capital_excess, capital_d)

    capital_to_output, credit_to_output, share_of_workers, share_of_self_employed, share_of_employers, var_log_c, occ_trans_matrix, var_log_income, gini_y_w, gini_y_ent, agg_capital, agg_output, agg_credit, avg_log_income, avg_log_c, additional_results = quick_calculate_results(fixed_occ_shares, number_asset_grid,asset_grid,policy, a1_indices, lottery_prob_1, lottery_prob_2,density_distr,r,w, occ_choice, income, earnings, capital_d, credit, output, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, p_alpha, P_zeta, P_u, stat_P_u, P_alpha, number_u_nodes, cost_of_employing, z_m_nodes, z_w_nodes, lambda,delta,gamma,eta,theta,c_e, labour_d)
    #                       1                       2                                       3                           4
    # additional_results = agg_cost_of_employing, agg_cost_of_employing_as_share_of_output, number_non_zero_asset_grid, occ_mobility_matrices
    #           1           2               3       4           5           6           7       8           9           10              11              12              13              14                      15                  16                  17          18          19              20         21          22        23       24          25          26            27          28      29        30          31      32      33      34            35              36      37              38          39          40              41              42                  43                 44 45 46   47   48   49  50
    return [a_max, number_asset_grid, asset_grid, policy, density_distr, K_demand, K_supply, L_demand, L_supply, Capital_excess, Labor_excess, capital_to_output, credit_to_output, share_of_workers, share_of_self_employed, share_of_employers, var_log_c, occ_trans_matrix, var_log_income, gini_y_w, gini_y_ent, occ_choice, income, earnings, capital_excess, capital_d, credit, labour_excess, labour_d, labour_s, deposit, output, value, agg_capital, agg_output, agg_credit, avg_log_income, avg_log_c, a1_indices, lottery_prob_1, lottery_prob_2, cost_of_employing, additional_results, r, w, gps, gaps, ps, ao, a_nodes]

end
