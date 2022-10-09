##############
# Script to generate the probabilities of occupational mobility
##############

include("Functions/print_sameline.jl")
include("Functions/profit.jl")
############################
######### Libraries ########
print_sameline("Initialising libraries")
using Plots

#using SchumakerSpline
#using Statistics
using JLD2
using ProgressMeter
using SchumakerSpline

CODENAME = "SS_2092_69"
country = "Italy"
LOCAL_DIR = "$(@__DIR__)/Results/Transitional_dynamics_big_grid/$(country)_$(CODENAME)/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Transitional_dynamics_big_grid\\$(country)_$(CODENAME)\\"
end

# analytical part
function f()


    TIME_PERIODS = 50#100
    @load "$(LOCAL_DIR)trans_$(TIME_PERIODS)_results.jld2" trans_res ss_1 ss_2

    global_params = ss_1[1][46]
    global_approx_params = ss_1[1][47]
    model_params = ss_1[1][48]
    approx_object = ss_1[1][49]

    delta = model_params[3]
    gamma = model_params[4]
    eta = model_params[5]
    theta = model_params[6]
    c_e = model_params[7]

    number_u_nodes = approx_object[7]
    number_zeta_nodes = global_approx_params[4]
    number_alpha_m_nodes = global_approx_params[5]
    number_alpha_w_nodes = global_approx_params[6]
    number_asset_grid = trans_res[9]
    asset_grid = trans_res[10]
    T = TIME_PERIODS

    z_m_nodes = approx_object[1]
    z_w_nodes = approx_object[2]
    number_u_nodes = approx_object[7]

    P_u = approx_object[4]
    p_alpha = model_params[13]
    P_zeta = approx_object[3]
    stat_P_u = approx_object[5]
    P_alpha = approx_object[6]

    lambda_s = trans_res[2]
    r_s = trans_res[3]
    w_s = trans_res[4]
    capital_s_distr_s = trans_res[5][7]
    policy_s = trans_res[5][8]

    function calculate_from_policies(policies)
        a1_indices = Array{Int64}(undef,T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        lottery_prob_1 = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)
        lottery_prob_2 = zeros(T, number_asset_grid,number_u_nodes,number_zeta_nodes,number_alpha_m_nodes,number_alpha_w_nodes)

        p = Progress((T*number_asset_grid*number_u_nodes*number_zeta_nodes*number_alpha_m_nodes*number_alpha_w_nodes), dt=0.5,
                     barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                     barlen=25)
        Threads.@threads for (t,(a_i,(u_i,(zeta_i,(alpha_m_i,alpha_w_i))))) in collect(Iterators.product(1:T,Iterators.product(1:number_asset_grid,Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))))
            try
                a1 = policies[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
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
            next!(p)
        end

        return a1_indices, lottery_prob_1, lottery_prob_2
    end
    a1_indices, lottery_prob_1, lottery_prob_2 = calculate_from_policies(policy_s)
    println_sameline("calculate_from_policies - Done")

    init_density_distr = capital_s_distr_s[1,:,:,:,:,:]
    init_occ_choice,  = compute_income_profile(asset_grid,number_asset_grid,r_s[1], w_s[1], number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, lambda_s[1], delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)
    init_w_choice = Float64.(init_occ_choice.==1.0)
    init_sp_choice = Float64.(init_occ_choice.==2.0)
    init_emp_choice = Float64.(init_occ_choice.==3.0)

    init_density_distr_w = init_density_distr.*init_w_choice
    init_density_distr_sp = init_density_distr.*init_sp_choice
    init_density_distr_emp = init_density_distr.*init_emp_choice

    # compute cumulative conditional occupational choice probability after reform
    cum_density_distr_w = copy(init_density_distr_w)
    cum_density_distr_sp = copy(init_density_distr_sp)
    cum_density_distr_emp = copy(init_density_distr_emp)

    cum_occ_trans_matrix = zeros(T+1,3,3)

    # inequality measures for occupations fixed in time t=1 after reform
    function compute_household_measures(policy, r, w, lambda, a_grid)

        occ_choice, income, earnings, = compute_income_profile(a_grid,number_asset_grid,r, w, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes, lambda, delta, gamma, eta, theta, c_e, z_m_nodes, z_w_nodes, number_u_nodes)

        consumption = (income.-policy)
        wealth = ones(size(init_density_distr)).*asset_grid
        income = income .- wealth

        w_choice = Float64.(occ_choice.==1.0)
        sp_choice = Float64.(occ_choice.==2.0)
        emp_choice = Float64.(occ_choice.==3.0)

        return consumption, earnings, income, wealth, occ_choice,w_choice,sp_choice,emp_choice
    end
    # time periods, [consumption, earnings, income, wealth], [W,SP,EMP]
    cum_mean_household_measures = zeros(T+1,4,3)
    # time periods, [consumption, earnings, income, wealth], [W,SP,EMP]*[W,SP,EMP]
    cum_mean_trans_household_measures = zeros(T+1,4,3,3)

    # compute conditional occupational choice probability after reform
    density_distr_w = copy(init_density_distr_w)
    density_distr_sp = copy(init_density_distr_sp)
    density_distr_emp = copy(init_density_distr_emp)
    occ_trans_matrix = zeros(T+1,3,3)

    # Compute conditional occupational choice probability for 50 years in a steady state pre-reform
    ss_density_distr = ss_1[1][5]
    ss_occ_choice = ss_1[1][22]
    ss_w_choice = Float64.(ss_occ_choice.==1.0)
    ss_sp_choice = Float64.(ss_occ_choice.==2.0)
    ss_emp_choice = Float64.(ss_occ_choice.==3.0)
    ss_init_density_distr_w = ss_density_distr.*ss_w_choice
    ss_init_density_distr_sp = ss_density_distr.*ss_sp_choice
    ss_init_density_distr_emp = ss_density_distr.*ss_emp_choice
    ss_policy = ss_1[1][4]
    ss_a1_indices = ss_1[1][39]
    ss_lottery_prob_1 = ss_1[1][40]
    ss_lottery_prob_2 = ss_1[1][41]

    ss_cum_density_distr_w = copy(ss_init_density_distr_w)
    ss_cum_density_distr_sp = copy(ss_init_density_distr_sp)
    ss_cum_density_distr_emp = copy(ss_init_density_distr_emp)
    ss_cum_occ_trans_matrix = zeros(T+1,3,3)

    ss_consumption, ss_earnings, ss_income, ss_wealth,  = compute_household_measures(ss_policy, ss_1[1][44], ss_1[1][45], lambda_s[1], ss_1[1][3])
    # time periods, [consumption, earnings, income, wealth], [W,SP,EMP]
    ss_cum_mean_household_measures = zeros(T+1,4,3)
    # time periods, [consumption, earnings, income, wealth], [W,SP,EMP]*[W,SP,EMP]
    ss_cum_mean_trans_household_measures = zeros(T+1,4,3,3)

    # t = 0 & iterator = 1
    ss_cum_occ_trans_matrix[1,1,1] = 1.0
    ss_cum_occ_trans_matrix[1,2,2] = 1.0
    ss_cum_occ_trans_matrix[1,3,3] = 1.0

    for (m,o) = collect(Iterators.product(1:4,1:3))
        dd = ss_cum_density_distr_w
        if o == 2
            dd = ss_cum_density_distr_sp
        elseif o == 3
            dd = ss_cum_density_distr_emp
        end
        measure = ss_consumption
        if m==2
            measure = ss_earnings
        elseif m==3
            measure = ss_income
        elseif m==4
            measure = ss_wealth
        end
        ss_cum_mean_household_measures[1,m,o] = sum(measure.*dd)/sum(dd)
        for o2 = 1:3
            occ2 = ss_w_choice
            if o2 == 2
                occ2 = ss_sp_choice
            elseif o2 == 3
                occ2 = ss_emp_choice
            end
            if sum(dd.*occ2) == 0.0
                ss_cum_mean_trans_household_measures[1,m,o,o2] = 0.0
            else
                ss_cum_mean_trans_household_measures[1,m,o,o2] = sum(measure.*dd.*occ2)/sum(dd.*occ2)
            end
            #println_sameline(sum(dd.*occ2))
        end
    end

    cum_occ_trans_matrix[1,1,1] = 1.0
    cum_occ_trans_matrix[1,2,2] = 1.0
    cum_occ_trans_matrix[1,3,3] = 1.0

    cum_mean_household_measures[1,:,:] = copy(ss_cum_mean_household_measures[1,:,:])
    cum_mean_trans_household_measures[1,:,:,:] = copy(ss_cum_mean_trans_household_measures[1,:,:,:])

    occ_trans_matrix[1,:,:] = copy(ss_1[1][18])

    # t = 1 & iterator = 2

    ss_cum_density_distr_w_pr = copy(ss_cum_density_distr_w).*0.0
    ss_cum_density_distr_sp_pr = copy(ss_cum_density_distr_sp).*0.0
    ss_cum_density_distr_emp_pr = copy(ss_cum_density_distr_emp).*0.0

    Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
        for a_i in 1:number_asset_grid
            ss_j_1 = ss_a1_indices[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

            ss_cum_density_distr_w_pr[ss_j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
            ss_cum_density_distr_sp_pr[ss_j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
            ss_cum_density_distr_emp_pr[ss_j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

            if ss_j_1 != number_asset_grid
                ss_cum_density_distr_w_pr[ss_j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                ss_cum_density_distr_sp_pr[ss_j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                ss_cum_density_distr_emp_pr[ss_j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

            end
        end
    end

    ss_cum_density_distr_marginal_assets_w = sum(sum(sum(sum(ss_cum_density_distr_w_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
    ss_cum_density_distr_marginal_assets_sp = sum(sum(sum(sum(ss_cum_density_distr_sp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
    ss_cum_density_distr_marginal_assets_emp = sum(sum(sum(sum(ss_cum_density_distr_emp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]

    Threads.@threads for (alpha_m_i,alpha_w_i) in collect(Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))
        ss_cum_temp_w = sum(ss_cum_density_distr_w_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
        ss_cum_temp_w = ss_cum_temp_w*P_u
        ss_cum_temp_sp = sum(ss_cum_density_distr_sp_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
        ss_cum_temp_sp = ss_cum_temp_sp*P_u
        ss_cum_temp_emp = sum(ss_cum_density_distr_emp_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
        ss_cum_temp_emp = ss_cum_temp_emp*P_u

        for zeta_prime_i in 1:number_zeta_nodes
            ss_cum_density_distr_w_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*ss_cum_temp_w
            ss_cum_density_distr_sp_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*ss_cum_temp_sp
            ss_cum_density_distr_emp_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*ss_cum_temp_emp

        end
    end

    Threads.@threads for (u_prime_i,(zeta_prime_i,(alpha_m_prime_i,alpha_w_prime_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
        ss_cum_density_distr_w_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*ss_cum_density_distr_marginal_assets_w
        ss_cum_density_distr_sp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*ss_cum_density_distr_marginal_assets_sp
        ss_cum_density_distr_emp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*ss_cum_density_distr_marginal_assets_emp
    end

    ss_cum_occ_trans_matrix[2, 1,1] = sum(ss_w_choice.*ss_cum_density_distr_w_pr)/sum(ss_init_density_distr_w)
    ss_cum_occ_trans_matrix[2, 1,2] = sum(ss_sp_choice.*ss_cum_density_distr_w_pr)/sum(ss_init_density_distr_w)
    ss_cum_occ_trans_matrix[2, 1,3] = sum(ss_emp_choice.*ss_cum_density_distr_w_pr)/sum(ss_init_density_distr_w)

    ss_cum_occ_trans_matrix[2, 2,1] = sum(ss_w_choice.*ss_cum_density_distr_sp_pr)/sum(ss_init_density_distr_sp)
    ss_cum_occ_trans_matrix[2, 2,2] = sum(ss_sp_choice.*ss_cum_density_distr_sp_pr)/sum(ss_init_density_distr_sp)
    ss_cum_occ_trans_matrix[2, 2,3] = sum(ss_emp_choice.*ss_cum_density_distr_sp_pr)/sum(ss_init_density_distr_sp)

    ss_cum_occ_trans_matrix[2, 3,1] = sum(ss_w_choice.*ss_cum_density_distr_emp_pr)/sum(ss_init_density_distr_emp)
    ss_cum_occ_trans_matrix[2, 3,2] = sum(ss_sp_choice.*ss_cum_density_distr_emp_pr)/sum(ss_init_density_distr_emp)
    ss_cum_occ_trans_matrix[2, 3,3] = sum(ss_emp_choice.*ss_cum_density_distr_emp_pr)/sum(ss_init_density_distr_emp)

    ss_cum_density_distr_w = copy(ss_cum_density_distr_w_pr)
    ss_cum_density_distr_sp = copy(ss_cum_density_distr_sp_pr)
    ss_cum_density_distr_emp = copy(ss_cum_density_distr_emp_pr)

    for (m,o) = collect(Iterators.product(1:4,1:3))
        dd = ss_cum_density_distr_w
        if o == 2
            dd = ss_cum_density_distr_sp
        elseif o == 3
            dd = ss_cum_density_distr_emp
        end
        measure = ss_consumption
        if m==2
            measure = ss_earnings
        elseif m==3
            measure = ss_income
        elseif m==4
            measure = ss_wealth
        end
        ss_cum_mean_household_measures[2,m,o] = sum(measure.*dd)/sum(dd)
        for o2 = 1:3
            occ2 = ss_w_choice
            if o2 == 2
                occ2 = ss_sp_choice
            elseif o2 == 3
                occ2 = ss_emp_choice
            end
            ss_cum_mean_trans_household_measures[2,m,o,o2] = sum(measure.*dd.*occ2)/sum(dd.*occ2)
        end
    end

    cum_occ_trans_matrix[2,:,:] = copy(ss_cum_occ_trans_matrix[2,:,:])

    Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
        temp_d_w = Schumaker(ss_1[1][3], cumsum(ss_cum_density_distr_w[:,u_i,zeta_i,alpha_m_i,alpha_w_i]); extrapolation=(Linear,Linear))
        cum_density_distr_w[:,u_i,zeta_i,alpha_m_i,alpha_w_i] = diff([0.0;evaluate.(temp_d_w, asset_grid)])
        cum_density_distr_w[:,u_i,zeta_i,alpha_m_i,alpha_w_i] .*= sum(ss_cum_density_distr_w[:,u_i,zeta_i,alpha_m_i,alpha_w_i])/sum(cum_density_distr_w[:,u_i,zeta_i,alpha_m_i,alpha_w_i])

        temp_d_sp = Schumaker(ss_1[1][3], cumsum(ss_cum_density_distr_sp[:,u_i,zeta_i,alpha_m_i,alpha_w_i]); extrapolation=(Linear,Linear))
        cum_density_distr_sp[:,u_i,zeta_i,alpha_m_i,alpha_w_i] = diff([0.0;evaluate.(temp_d_sp, asset_grid)])
        cum_density_distr_sp[:,u_i,zeta_i,alpha_m_i,alpha_w_i] .*= sum(ss_cum_density_distr_sp[:,u_i,zeta_i,alpha_m_i,alpha_w_i])/sum(cum_density_distr_sp[:,u_i,zeta_i,alpha_m_i,alpha_w_i])

        temp_d_emp = Schumaker(ss_1[1][3], cumsum(ss_cum_density_distr_emp[:,u_i,zeta_i,alpha_m_i,alpha_w_i]); extrapolation=(Linear,Linear))
        cum_density_distr_emp[:,u_i,zeta_i,alpha_m_i,alpha_w_i] = diff([0.0;evaluate.(temp_d_emp, asset_grid)])
        cum_density_distr_emp[:,u_i,zeta_i,alpha_m_i,alpha_w_i] .*= sum(ss_cum_density_distr_emp[:,u_i,zeta_i,alpha_m_i,alpha_w_i])/sum(cum_density_distr_emp[:,u_i,zeta_i,alpha_m_i,alpha_w_i])

    end

    cum_consumption, cum_earnings, cum_income, cum_wealth, cum_occ_choice,cum_w_choice,cum_sp_choice,cum_emp_choice = compute_household_measures(policy_s[1,:,:,:,:,:], r_s[1], w_s[1], lambda_s[1], asset_grid)
    for (m,o) = collect(Iterators.product(1:4,1:3))
        dd = cum_density_distr_w
        if o == 2
            dd = cum_density_distr_sp
        elseif o == 3
            dd = cum_density_distr_emp
        end
        measure = cum_consumption
        if m==2
            measure = cum_earnings
        elseif m==3
            measure = cum_income
        elseif m==4
            measure = cum_wealth
        end
        cum_mean_household_measures[2,m,o] = sum(measure.*dd)/sum(dd)
        for o2 = 1:3
            occ2 = cum_w_choice
            if o2 == 2
                occ2 = cum_sp_choice
            elseif o2 == 3
                occ2 = cum_emp_choice
            end
            cum_mean_trans_household_measures[2,m,o,o2] = sum(measure.*dd.*occ2)/sum(dd.*occ2)
        end
    end

    occ_trans_matrix[2,:,:] = copy(ss_1[1][18])

    # t = 2:T
    p = Progress((T-1)*((number_asset_grid*number_u_nodes+2)*(number_zeta_nodes*number_alpha_m_nodes*number_alpha_w_nodes)), dt=0.5,
                 barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                 barlen=100)
    for t=1:T-1 # & i = 3:T+1
        i = t+2

        cum_density_distr_w_pr = copy(cum_density_distr_w).*0.0
        cum_density_distr_sp_pr = copy(cum_density_distr_sp).*0.0
        cum_density_distr_emp_pr = copy(cum_density_distr_emp).*0.0

        ss_cum_density_distr_w_pr = copy(ss_cum_density_distr_w).*0.0
        ss_cum_density_distr_sp_pr = copy(ss_cum_density_distr_sp).*0.0
        ss_cum_density_distr_emp_pr = copy(ss_cum_density_distr_emp).*0.0

        density_distr_w_pr = copy(density_distr_w).*0.0
        density_distr_sp_pr = copy(density_distr_sp).*0.0
        density_distr_emp_pr = copy(density_distr_emp).*0.0

        Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            for a_i in 1:number_asset_grid
                j_1 = a1_indices[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                cum_density_distr_w_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += cum_density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                cum_density_distr_sp_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += cum_density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                cum_density_distr_emp_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += cum_density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                density_distr_w_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                density_distr_sp_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                density_distr_emp_pr[j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_1[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                if j_1 != number_asset_grid
                    cum_density_distr_w_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += cum_density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    cum_density_distr_sp_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += cum_density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    cum_density_distr_emp_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += cum_density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                    density_distr_w_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    density_distr_sp_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    density_distr_emp_pr[j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*lottery_prob_2[t, a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                end

                ss_j_1 = ss_a1_indices[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                ss_cum_density_distr_w_pr[ss_j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                ss_cum_density_distr_sp_pr[ss_j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                ss_cum_density_distr_emp_pr[ss_j_1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_1[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                if ss_j_1 != number_asset_grid
                    ss_cum_density_distr_w_pr[ss_j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_w[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    ss_cum_density_distr_sp_pr[ss_j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_sp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
                    ss_cum_density_distr_emp_pr[ss_j_1+1,u_i,zeta_i,alpha_m_i,alpha_w_i] += ss_cum_density_distr_emp[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]*ss_lottery_prob_2[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]

                end
                next!(p)
            end
        end

        cum_density_distr_marginal_assets_w = sum(sum(sum(sum(cum_density_distr_w_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
        cum_density_distr_marginal_assets_sp = sum(sum(sum(sum(cum_density_distr_sp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
        cum_density_distr_marginal_assets_emp = sum(sum(sum(sum(cum_density_distr_emp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]

        ss_cum_density_distr_marginal_assets_w = sum(sum(sum(sum(ss_cum_density_distr_w_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
        ss_cum_density_distr_marginal_assets_sp = sum(sum(sum(sum(ss_cum_density_distr_sp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
        ss_cum_density_distr_marginal_assets_emp = sum(sum(sum(sum(ss_cum_density_distr_emp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]

        density_distr_marginal_assets_w = sum(sum(sum(sum(density_distr_w_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
        density_distr_marginal_assets_sp = sum(sum(sum(sum(density_distr_sp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]
        density_distr_marginal_assets_emp = sum(sum(sum(sum(density_distr_emp_pr,dims=5)[:,:,:,:,1],dims=4)[:,:,:,1],dims=3)[:,:,1],dims=2)[:,1]

        Threads.@threads for (alpha_m_i,alpha_w_i) in collect(Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))
            cum_temp_w = sum(cum_density_distr_w_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            cum_temp_w = cum_temp_w*P_u
            cum_temp_sp = sum(cum_density_distr_sp_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            cum_temp_sp = cum_temp_sp*P_u
            cum_temp_emp = sum(cum_density_distr_emp_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            cum_temp_emp = cum_temp_emp*P_u

            ss_cum_temp_w = sum(ss_cum_density_distr_w_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            ss_cum_temp_w = ss_cum_temp_w*P_u
            ss_cum_temp_sp = sum(ss_cum_density_distr_sp_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            ss_cum_temp_sp = ss_cum_temp_sp*P_u
            ss_cum_temp_emp = sum(ss_cum_density_distr_emp_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            ss_cum_temp_emp = ss_cum_temp_emp*P_u

            temp_w = sum(density_distr_w_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            temp_w = temp_w*P_u
            temp_sp = sum(density_distr_sp_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            temp_sp = temp_sp*P_u
            temp_emp = sum(density_distr_emp_pr[:,:,:,alpha_m_i,alpha_w_i],dims=3)[:,:,1]
            temp_emp = temp_emp*P_u

            for zeta_prime_i in 1:number_zeta_nodes
                cum_density_distr_w_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*cum_temp_w
                cum_density_distr_sp_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*cum_temp_sp
                cum_density_distr_emp_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*cum_temp_emp

                ss_cum_density_distr_w_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*ss_cum_temp_w
                ss_cum_density_distr_sp_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*ss_cum_temp_sp
                ss_cum_density_distr_emp_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*ss_cum_temp_emp

                density_distr_w_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_w
                density_distr_sp_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_sp
                density_distr_emp_pr[:,:,zeta_prime_i,alpha_m_i,alpha_w_i] = ((1-p_alpha)*P_zeta[zeta_prime_i]).*temp_emp

                next!(p)
            end
        end

        Threads.@threads for (u_prime_i,(zeta_prime_i,(alpha_m_prime_i,alpha_w_prime_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
            cum_density_distr_w_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*cum_density_distr_marginal_assets_w
            cum_density_distr_sp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*cum_density_distr_marginal_assets_sp
            cum_density_distr_emp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*cum_density_distr_marginal_assets_emp

            ss_cum_density_distr_w_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*ss_cum_density_distr_marginal_assets_w
            ss_cum_density_distr_sp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*ss_cum_density_distr_marginal_assets_sp
            ss_cum_density_distr_emp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*ss_cum_density_distr_marginal_assets_emp

            density_distr_w_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_assets_w
            density_distr_sp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_assets_sp
            density_distr_emp_pr[:,u_prime_i,zeta_prime_i,alpha_m_prime_i,alpha_w_prime_i] .+= (p_alpha*stat_P_u[u_prime_i]*P_zeta[zeta_prime_i]*P_alpha[alpha_m_prime_i,alpha_w_prime_i]).*density_distr_marginal_assets_emp

            next!(p)
        end

        # compute inequality measures
        cum_consumption, cum_earnings, cum_income, cum_wealth, occ_choice,cum_w_choice,cum_sp_choice,cum_emp_choice = compute_household_measures(policy_s[t+1,:,:,:,:,:], r_s[t+1], w_s[t+1], lambda_s[t+1], asset_grid)
        for (m,o) = collect(Iterators.product(1:4,1:3))
            dd = cum_density_distr_w_pr
            if o == 2
                dd = cum_density_distr_sp_pr
            elseif o == 3
                dd = cum_density_distr_emp_pr
            end
            measure = cum_consumption
            if m==2
                measure = cum_earnings
            elseif m==3
                measure = cum_income
            elseif m==4
                measure = cum_wealth
            end
            cum_mean_household_measures[i,m,o] = sum(measure.*dd)/sum(dd)
            for o2 = 1:3
                occ2 = cum_w_choice
                if o2 == 2
                    occ2 = cum_sp_choice
                elseif o2 == 3
                    occ2 = cum_emp_choice
                end
                cum_mean_trans_household_measures[i,m,o,o2] = sum(measure.*dd.*occ2)/sum(dd.*occ2)
            end
        end

        w_choice = Float64.(occ_choice.==1.0)
        sp_choice = Float64.(occ_choice.==2.0)
        emp_choice = Float64.(occ_choice.==3.0)

        # conditional cumulative probability of households - P( ot = OCCj | o1 = OCCi )
        # W->W
        cum_occ_trans_matrix[i, 1,1] = sum(w_choice.*cum_density_distr_w_pr)/sum(init_density_distr_w)
        # W->SP
        cum_occ_trans_matrix[i, 1,2] = sum(sp_choice.*cum_density_distr_w_pr)/sum(init_density_distr_w)
        # W->EMP
        cum_occ_trans_matrix[i, 1,3] = sum(emp_choice.*cum_density_distr_w_pr)/sum(init_density_distr_w)

        # SP->W
        cum_occ_trans_matrix[i, 2,1] = sum(w_choice.*cum_density_distr_sp_pr)/sum(init_density_distr_sp)
        # W->SP
        cum_occ_trans_matrix[i, 2,2] = sum(sp_choice.*cum_density_distr_sp_pr)/sum(init_density_distr_sp)
        # W->EMP
        cum_occ_trans_matrix[i, 2,3] = sum(emp_choice.*cum_density_distr_sp_pr)/sum(init_density_distr_sp)

        # EMP->W
        cum_occ_trans_matrix[i, 3,1] = sum(w_choice.*cum_density_distr_emp_pr)/sum(init_density_distr_emp)
        # EMP->SP
        cum_occ_trans_matrix[i, 3,2] = sum(sp_choice.*cum_density_distr_emp_pr)/sum(init_density_distr_emp)
        # EMP->EMP
        cum_occ_trans_matrix[i, 3,3] = sum(emp_choice.*cum_density_distr_emp_pr)/sum(init_density_distr_emp)

        # conditional cumulative occupational choice probability for steady state
        ss_cum_occ_trans_matrix[i, 1,1] = sum(ss_w_choice.*ss_cum_density_distr_w_pr)/sum(ss_init_density_distr_w)
        ss_cum_occ_trans_matrix[i, 1,2] = sum(ss_sp_choice.*ss_cum_density_distr_w_pr)/sum(ss_init_density_distr_w)
        ss_cum_occ_trans_matrix[i, 1,3] = sum(ss_emp_choice.*ss_cum_density_distr_w_pr)/sum(ss_init_density_distr_w)

        ss_cum_occ_trans_matrix[i, 2,1] = sum(ss_w_choice.*ss_cum_density_distr_sp_pr)/sum(ss_init_density_distr_sp)
        ss_cum_occ_trans_matrix[i, 2,2] = sum(ss_sp_choice.*ss_cum_density_distr_sp_pr)/sum(ss_init_density_distr_sp)
        ss_cum_occ_trans_matrix[i, 2,3] = sum(ss_emp_choice.*ss_cum_density_distr_sp_pr)/sum(ss_init_density_distr_sp)

        ss_cum_occ_trans_matrix[i, 3,1] = sum(ss_w_choice.*ss_cum_density_distr_emp_pr)/sum(ss_init_density_distr_emp)
        ss_cum_occ_trans_matrix[i, 3,2] = sum(ss_sp_choice.*ss_cum_density_distr_emp_pr)/sum(ss_init_density_distr_emp)
        ss_cum_occ_trans_matrix[i, 3,3] = sum(ss_emp_choice.*ss_cum_density_distr_emp_pr)/sum(ss_init_density_distr_emp)

        # inequality measures for steady state
        for (m,o) = collect(Iterators.product(1:4,1:3))
            dd = ss_cum_density_distr_w_pr
            if o == 2
                dd = ss_cum_density_distr_sp_pr
            elseif o == 3
                dd = ss_cum_density_distr_emp_pr
            end
            measure = ss_consumption
            if m==2
                measure = ss_earnings
            elseif m==3
                measure = ss_income
            elseif m==4
                measure = ss_wealth
            end
            ss_cum_mean_household_measures[i,m,o] = sum(measure.*dd)/sum(dd)
            for o2 = 1:3
                occ2 = ss_w_choice
                if o2 == 2
                    occ2 = ss_sp_choice
                elseif o2 == 3
                    occ2 = ss_emp_choice
                end
                ss_cum_mean_trans_household_measures[i,m,o,o2] = sum(measure.*dd.*occ2)/sum(dd.*occ2)
            end
        end

        # conditional occupational choice probability for steady state
        occ_trans_matrix[i, 1,1] = sum(w_choice.*density_distr_w_pr)/sum(density_distr_w)
        occ_trans_matrix[i, 1,2] = sum(sp_choice.*density_distr_w_pr)/sum(density_distr_w)
        occ_trans_matrix[i, 1,3] = sum(emp_choice.*density_distr_w_pr)/sum(density_distr_w)

        occ_trans_matrix[i, 2,1] = sum(w_choice.*density_distr_sp_pr)/sum(density_distr_sp)
        occ_trans_matrix[i, 2,2] = sum(sp_choice.*density_distr_sp_pr)/sum(density_distr_sp)
        occ_trans_matrix[i, 2,3] = sum(emp_choice.*density_distr_sp_pr)/sum(density_distr_sp)

        occ_trans_matrix[i, 3,1] = sum(w_choice.*density_distr_emp_pr)/sum(density_distr_emp)
        occ_trans_matrix[i, 3,2] = sum(sp_choice.*density_distr_emp_pr)/sum(density_distr_emp)
        occ_trans_matrix[i, 3,3] = sum(emp_choice.*density_distr_emp_pr)/sum(density_distr_emp)


        # final in-loop step
        cum_density_distr_w = copy(cum_density_distr_w_pr)
        cum_density_distr_sp = copy(cum_density_distr_sp_pr)
        cum_density_distr_emp = copy(cum_density_distr_emp_pr)

        ss_cum_density_distr_w = copy(ss_cum_density_distr_w_pr)
        ss_cum_density_distr_sp = copy(ss_cum_density_distr_sp_pr)
        ss_cum_density_distr_emp = copy(ss_cum_density_distr_emp_pr)

        density_distr_w = w_choice.*(density_distr_w_pr.+density_distr_sp_pr.+density_distr_emp_pr)
        density_distr_sp = sp_choice.*(density_distr_w_pr.+density_distr_sp_pr.+density_distr_emp_pr)
        density_distr_emp = emp_choice.*(density_distr_w_pr.+density_distr_sp_pr.+density_distr_emp_pr)

        #temp_res = [sum(density_distr_w), sum(density_distr_sp), sum(density_distr_emp)]
        #println_sameline(temp_res)
    end

    println_sameline("Final cycle - Done")

    return cum_occ_trans_matrix,ss_cum_occ_trans_matrix, occ_trans_matrix, cum_mean_household_measures,ss_cum_mean_household_measures, cum_mean_trans_household_measures,ss_cum_mean_trans_household_measures
end
cum_occ_trans_matrix,ss_cum_occ_trans_matrix, occ_trans_matrix, cum_mean_household_measures,ss_cum_mean_household_measures, cum_mean_trans_household_measures,ss_cum_mean_trans_household_measures = f()

# graphing part
LOCAL_DIR_OCCUPATION = "$(LOCAL_DIR)/Occupation/"
if Sys.iswindows()
    LOCAL_DIR_OCCUPATION = "$(LOCAL_DIR)\\Occupation\\"
end

occ_text = ["W","SP","EMP"]

plts_cum_otm = Array{Any}(undef,3,3)
plts_cum_ssotm = Array{Any}(undef,3,3)
plts_cum_diffotm = Array{Any}(undef,3,3)
plts_cum_wssotm = Array{Any}(undef,3,3)

plts_otm = Array{Any}(undef,3,3)

for (i,j) = collect(Iterators.product(1:3,1:3))
    T = length(cum_occ_trans_matrix[:,i,j])-1
    plts_cum_otm[i,j] = plot(0:T, cum_occ_trans_matrix[1:T+1,i,j], legend=false,ylabel="$(occ_text[i])->$(occ_text[j])")
    plts_cum_ssotm[i,j] = plot(0:T, ss_cum_occ_trans_matrix[1:T+1,i,j], legend=false,ylabel="$(occ_text[i])->$(occ_text[j])")
    plts_cum_diffotm[i,j] = plot(0:T, cum_occ_trans_matrix[1:T+1,i,j].-ss_cum_occ_trans_matrix[1:T+1,i,j], legend=false,ylabel="$(occ_text[i])->$(occ_text[j])")
    plts_cum_wssotm[i,j] = plot(0:T, [cum_occ_trans_matrix[1:T+1,i,j] ss_cum_occ_trans_matrix[1:T+1,i,j]], legend=false,ylabel="$(occ_text[i])->$(occ_text[j])")

    plts_otm[i,j] = plot(0:T, occ_trans_matrix[1:T+1,i,j], legend=false,ylabel="$(occ_text[i])->$(occ_text[j])")

end
#=
# it shows the cumulative share of occupation1 at time 1 who become occupation2 at time t
plt_cum_otm = plot(plts_cum_otm[1,1],plts_cum_otm[1,2],plts_cum_otm[1,3], plts_cum_otm[2,1],plts_cum_otm[2,2],plts_cum_otm[2,3], plts_cum_otm[3,1],plts_cum_otm[3,2],plts_cum_otm[3,3], layout=(3,3))
display(plt_cum_otm)
savefig(plt_cum_otm,"$(LOCAL_DIR_OCCUPATION)time_cumulative_occupational_mobility.png")

# it shows the cumulative share of occupation1 at time 1 who become occupation2 at time t for pre-reform steady state
plt_cum_ssotm = plot(plts_cum_ssotm[1,1],plts_cum_ssotm[1,2],plts_cum_ssotm[1,3], plts_cum_ssotm[2,1],plts_cum_ssotm[2,2],plts_cum_ssotm[2,3], plts_cum_ssotm[3,1],plts_cum_ssotm[3,2],plts_cum_ssotm[3,3], layout=(3,3))
display(plt_cum_ssotm)
savefig(plt_cum_ssotm,"$(LOCAL_DIR_OCCUPATION)time_cumulative_occupational_mobility_ss.png")

# it shows the difference between cumulative occupational mobility after reform and pre-reform
plt_cum_diffotm = plot(plts_cum_diffotm[1,1],plts_cum_diffotm[1,2],plts_cum_diffotm[1,3], plts_cum_diffotm[2,1],plts_cum_diffotm[2,2],plts_cum_diffotm[2,3], plts_cum_diffotm[3,1],plts_cum_diffotm[3,2],plts_cum_diffotm[3,3], layout=(3,3))
display(plt_cum_diffotm)
savefig(plt_cum_diffotm,"$(LOCAL_DIR_OCCUPATION)time_cumulative_occupational_mobility_diff_with_ss.png")
=#
# it shows cumulative occupational mobilities after reform and pre-reform
plt_cum_wssotm = plot(plts_cum_wssotm[1,1],plts_cum_wssotm[1,2],plts_cum_wssotm[1,3], plts_cum_wssotm[2,1],plts_cum_wssotm[2,2],plts_cum_wssotm[2,3], plts_cum_wssotm[3,1],plts_cum_wssotm[3,2],plts_cum_wssotm[3,3], layout=(3,3))
display(plt_cum_wssotm)
savefig(plt_cum_wssotm,"$(LOCAL_DIR_OCCUPATION)time_cumulative_occupational_mobility_with_ss.png")

# it shows the transitional share of occupation1 at time 1 who become occupation2 at time t
plt_otm = plot(plts_otm[1,1],plts_otm[1,2],plts_otm[1,3], plts_otm[2,1],plts_otm[2,2],plts_otm[2,3], plts_otm[3,1],plts_otm[3,2],plts_otm[3,3], layout=(3,3))
display(plt_otm)
savefig(plt_otm,"$(LOCAL_DIR_OCCUPATION)time_occupational_mobility.png")

# it shows the inequality measures for fixed occupation in time t=1

LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)/Inequality/Transition/"
if Sys.iswindows()
    LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)\\Inequality\\Transition\\"
end
mkpath(LOCAL_DIR_INEQUALITY)
#cum_mean_household_measures
measures = ["consumption","earnings","income","wealth"]
for m = 1:4
    plts_cum_hm = Array{Any}(undef,3)
    plts_cum_sshm = Array{Any}(undef,3)
    plts_cum_diffhm = Array{Any}(undef,3)
    plts_cum_wsshm = Array{Any}(undef,3)
    for o = 1:3
        T = length(cum_mean_household_measures[:,m,o])-1
        plts_cum_hm[o] = plot(0:T, cum_mean_household_measures[1:T+1,m,o], legend=false,ylabel="$(occ_text[o])'s mean $(measures[m])")
        plts_cum_sshm[o] = plot(0:T, ss_cum_mean_household_measures[1:T+1,m,o], legend=false,ylabel="$(occ_text[o])'s mean $(measures[m]) (pre-reform)")

        plts_cum_diffhm[o] = plot(0:T, cum_mean_household_measures[1:T+1,m,o].-ss_cum_mean_household_measures[1:T+1,m,o], legend=false,ylabel="$(occ_text[o])'s mean $(measures[m]) (diff)")
        plts_cum_wsshm[o] = plot(0:T, [cum_mean_household_measures[1:T+1,m,o] ss_cum_mean_household_measures[1:T+1,m,o]], legend=false,ylabel="$(occ_text[o])'s mean $(measures[m])")

    end
    #=
    plt_cum_hm = plot(plts_cum_hm[1],plts_cum_hm[2],plts_cum_hm[3], layout=(1,3))
    display(plt_cum_hm)
    savefig(plt_cum_hm,"$(LOCAL_DIR_INEQUALITY)time_mean_$(measures[m]).png")

    plt_cum_sshm = plot(plts_cum_sshm[1],plts_cum_sshm[2],plts_cum_sshm[3], layout=(1,3))
    display(plt_cum_sshm)
    savefig(plt_cum_sshm,"$(LOCAL_DIR_INEQUALITY)time_mean_$(measures[m])_ss.png")

    plt_cum_diffhm = plot(plts_cum_diffhm[1],plts_cum_diffhm[2],plts_cum_diffhm[3], layout=(1,3))
    display(plt_cum_diffhm)
    savefig(plt_cum_diffhm,"$(LOCAL_DIR_INEQUALITY)time_mean_$(measures[m])_diff.png")
    =#
    plt_cum_wsshm = plot(plts_cum_wsshm[1],plts_cum_wsshm[2],plts_cum_wsshm[3], layout=(1,3))
    display(plt_cum_wsshm)
    savefig(plt_cum_wsshm,"$(LOCAL_DIR_INEQUALITY)time_mean_$(measures[m])_with_ss.png")

    plts_cum_mthm = Array{Any}(undef,3,3)
    plts_cum_ssmthm = Array{Any}(undef,3,3)
    plts_cum_diffmthm = Array{Any}(undef,3,3)
    plts_cum_wssmthm = Array{Any}(undef,3,3)
    for (i,j) = collect(Iterators.product(1:3,1:3))
        T = length(cum_mean_trans_household_measures[:,m,i,j])-1
        plts_cum_mthm[i,j] = plot(0:T, cum_mean_trans_household_measures[1:T+1,m,i,j], legend=false,ylabel="$(occ_text[i])->$(occ_text[j])")
        plts_cum_ssmthm[i,j] = plot(0:T, ss_cum_mean_trans_household_measures[1:T+1,m,i,j], legend=false,ylabel="$(occ_text[i])->$(occ_text[j])")
        plts_cum_diffmthm[i,j] = plot(0:T, cum_mean_trans_household_measures[1:T+1,m,i,j].-ss_cum_mean_trans_household_measures[1:T+1,m,i,j], legend=false,ylabel="$(occ_text[i])->$(occ_text[j])")
        plts_cum_wssmthm[i,j] = plot(0:T, [cum_mean_trans_household_measures[1:T+1,m,i,j] ss_cum_mean_trans_household_measures[1:T+1,m,i,j]], legend=false,ylabel="$(occ_text[i])->$(occ_text[j])")

    end

    plt_cum_wssmthm = plot(plts_cum_wssmthm[1,1],plts_cum_wssmthm[1,2],plts_cum_wssmthm[1,3], plts_cum_wssmthm[2,1],plts_cum_wssmthm[2,2],plts_cum_wssmthm[2,3], plts_cum_wssmthm[3,1],plts_cum_wssmthm[3,2],plts_cum_wssmthm[3,3], layout=(3,3))
    display(plt_cum_wssmthm)
    savefig(plt_cum_wssmthm,"$(LOCAL_DIR_OCCUPATION)time_cumulative_mean_$(measures[m])_with_ss.png")
end
