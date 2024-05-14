include("Functions/print_sameline.jl")
############################
######### Libraries ########
print_sameline("Initialising libraries")
using Plots
#using PyPlot

using Interpolations
using Statistics
using JLD2
using ProgressMeter

country = "Italy"

# global parameters of the approximation objects
#                               1               2               3               4                       5                   6
#                       number_a_nodes, number_u_m_nodes, number_u_w_nodes, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes
GLOBAL_APPROX_PARAMS = [69,             3,                3,                3,                 6,                    3]

# parameters of the model's economy (Italy)
LAMBDA =            1.665907
LOCAL_DIR = "$(@__DIR__)/Results/Stationary/GE_lambda/"
if Sys.iswindows()
    LOCAL_DIR = "$(@__DIR__)\\Results\\Stationary\\GE_lambda\\"
end
global_approx_params = copy(GLOBAL_APPROX_PARAMS)

print_sameline("Loading data from Lambda_grid_big_grid_test")
#@load "$(LOCAL_DIR)SSS.jld2" SSS
SSS_names = readdir(LOCAL_DIR)
SSS_names = SSS_names[findall(x->occursin("SS_",x), SSS_names)]
SSS = Array{Any}(undef,length(SSS_names))
lambdas = zeros(length(SSS))
@showprogress for i = 1:length(SSS)
    @load "$(LOCAL_DIR)$(SSS_names[i])" SS
    SSS[i] = copy(SS)
    lambdas[i] = SS[5][1]
end
temp_is = sortperm(lambdas)
lambdas = lambdas[temp_is]
SSS = SSS[temp_is]

num_lambdas = length(lambdas)#20# #instead of length(lambda)
lambdas = lambdas[1:num_lambdas]
calibrated_lambda = findfirst(x -> x==LAMBDA, lambdas)

Rs = zeros(num_lambdas)
Ws = zeros(num_lambdas)
K_Ys = zeros(num_lambdas)
C_Ys = zeros(num_lambdas)
occWs = zeros(num_lambdas)
occSEs = zeros(num_lambdas)
occEMPs = zeros(num_lambdas)
occENTs = zeros(num_lambdas)
logcs = zeros(num_lambdas)
loges = zeros(num_lambdas)
giniWs = zeros(num_lambdas)
giniEnts = zeros(num_lambdas)

Outputs = zeros(num_lambdas)
Incomes = zeros(num_lambdas)
Consumptions = zeros(num_lambdas)
Income_to_outputs = zeros(num_lambdas)

# income, earnings, wealth, consumption
# All, W, SP, EMP, ENT
# lambdas
means = zeros(4,5,num_lambdas)
ginis = zeros(4,5,num_lambdas)
avglogs = zeros(4,5,num_lambdas)
varlogs = zeros(4,5,num_lambdas)
avgs = zeros(4,5,num_lambdas)
vars = zeros(4,5,num_lambdas)

function quantile_mean(stat, distr)
    distr ./= sum(distr)

    stat_grid = vec(stat)
    stat_distr = vec(distr)
    ids = findall(x->x>0.0, stat_distr)
    ud_stat_grid = stat_grid[ids]
    ud_stat_distr = stat_distr[ids]
    try
        stat_grid = exp.(collect(range(log(1); stop=log(maximum(ud_stat_grid)+1-minimum(ud_stat_grid)), length=length(stat[:,1,1,1,1,1])))).+minimum(ud_stat_grid).-1
    catch e
        temp111 = findmin(ud_stat_grid)
        display([ud_stat_grid[temp111[2]], ud_stat_distr[temp111[2]]])
        id111 = findall(x->x<=0.0, stat_grid)
        display(sum(stat_distr[id111]))
        throw(e)
    end
    stat_distr = zeros(length(stat_grid))
    Threads.@threads for i in 1:length(stat_grid)
        iters1 = []
        iters2 = []
        if i == 1
            iters2 = findall(x->stat_grid[i]<=x<=stat_grid[i+1],ud_stat_grid)
        elseif i == length(stat_grid)
            iters1 = findall(x->stat_grid[i-1]<=x<=stat_grid[i],ud_stat_grid)
        else
            iters1 = findall(x->stat_grid[i-1]<x<=stat_grid[i],ud_stat_grid)
            iters2 = findall(x->stat_grid[i]<x<stat_grid[i+1],ud_stat_grid)
        end
        if !isempty(iters1)
            cand1 = ud_stat_grid[iters1]
            lottery_prob1 = (cand1 .- stat_grid[i-1])./(stat_grid[i]-stat_grid[i-1])
            cand_distr1 = ud_stat_distr[iters1]
            stat_distr[i] += sum(cand_distr1.*lottery_prob1)
        end
        if !isempty(iters2)
            cand2 = ud_stat_grid[iters2]
            lottery_prob2 = (stat_grid[i+1] .- cand2)./(stat_grid[i+1]-stat_grid[i])
            cand_distr2 = ud_stat_distr[iters2]
            stat_distr[i] += sum(cand_distr2.*lottery_prob2)
        end

    end
    cdf = cumsum(stat_distr)

    idx = unique(z -> cdf[z], 1:length(cdf))
    cdf = cdf[idx]
    stat_grid = stat_grid[idx]
    stat_distr = stat_distr[idx]

    p0 = plot(stat_distr)
    max_a_i = findlast(x->x<1.0-1e-4, cdf)
    try
        cdf = cdf[1:max_a_i]./cdf[max_a_i]
    catch e1
        try
            max_a_i = length(cdf)-1
            cdf = cdf[1:max_a_i]./cdf[max_a_i]
        catch e
            display(cdf)
            throw(e)
        end
    end
    stat_grid = stat_grid[1:max_a_i]
    stat_distr = stat_distr[1:max_a_i]

    idx = unique(z -> cdf[z], 1:length(cdf))
    cdf = cdf[idx]
    stat_grid = stat_grid[idx]
    stat_distr = stat_distr[idx]

    p0 = plot(stat_distr)
    max_a_i = length(cdf)#findlast(x->x<1.0-1e-4, cdf)

    p1 = plot(cdf[1:max_a_i]./cdf[max_a_i])

    cdf_1 = []
    try
        #cdf_1 = Schumaker(cdf[1:max_a_i]./cdf[max_a_i],1:max_a_i; extrapolation=(Constant,Constant))
        cdf_1 = linear_interpolation(collect(cdf[1:max_a_i]./cdf[max_a_i]), collect(1:max_a_i), extrapolation_bc = Flat());
        p2 = plot(collect(range(0.0;stop=1.0,length=120)),cdf_1.(collect(range(0.0;stop=1.0,length=120))))
    catch e1
        try
            #cdf_1 = Schumaker(cdf[1:max_a_i]./cdf[max_a_i],1:max_a_i; extrapolation=(Constant,Constant))
            cdf_1 = linear_interpolation([0.0; collect(cdf[1:max_a_i]./cdf[max_a_i])], [1; collect(1:max_a_i)], extrapolation_bc = Flat());
            p2 = plot(collect(range(0.0;stop=1.0,length=120)),cdf_1.(collect(range(0.0;stop=1.0,length=120))))
        catch e2
            display(plot(p0,p1,legend=false))
            throw(e2)
        end
    end

    quantiles = collect(range(0.0;stop=1.0,length=6))
    #qwb = evaluate.(cdf_1,quantiles)    #quantile_wealth_bounds
    qwb = cdf_1.(quantiles)

    quantile_mean = zeros(5)
    for q = 1:5
        if qwb[q] > qwb[q+1]
            display(plot(p0,p1,p2,layout=(1,3),legend=false))
            display(qwb)
            throw(error)
        end
        _q = Int32(max(1,min(floor(qwb[q]),max_a_i)))
        q_ = Int32(max(1,min(ceil(qwb[q]),max_a_i)))
        _q1 = Int32(max(1,min(floor(qwb[q+1]),max_a_i)))
        q1_ = Int32(max(1,min(ceil(qwb[q+1]),max_a_i)))

        m_w_down_down = sum(stat_grid[_q:_q1] .* stat_distr[_q:_q1] / sum(stat_distr[_q:_q1]) )
        m_w_down_up = sum(stat_grid[_q:q1_] .* stat_distr[_q:q1_] / sum(stat_distr[_q:q1_]) )
        m_w_up_down = sum(stat_grid[q_:_q1] .* stat_distr[q_:_q1] / sum(stat_distr[q_:_q1]) )
        m_w_up_up = sum(stat_grid[q_:q1_] .* stat_distr[q_:q1_] / sum(stat_distr[q_:q1_]) )

        quantile_mean[q] = (m_w_down_down+m_w_down_up+m_w_up_down+m_w_up_up)/4
    end
    return quantile_mean
end
# 5 quantiles
# All, W, SP, EMP
# income, earnings, wealth, consumption
quantile_means = zeros(5,4,4,num_lambdas)

# All, W, SP, EMP
Capital = zeros(4,num_lambdas)

# ENT, SP, EMP
TFPis = zeros(3,num_lambdas)
TFPds = zeros(3,num_lambdas)

mean_MPL = zeros(3,num_lambdas)
var_MPL = zeros(3,num_lambdas)

mean_MPK = zeros(3,num_lambdas)
var_MPK = zeros(3,num_lambdas)

share_unbound = zeros(3,num_lambdas)

# W, ENT, SP,EMP
avg_w_skill = zeros(4,num_lambdas)
avg_m_skill = zeros(4,num_lambdas)

var_w_skill = zeros(4,num_lambdas)
var_m_skill = zeros(4,num_lambdas)

sum_w_skill = zeros(4,num_lambdas)
sum_m_skill = zeros(4,num_lambdas)

ratio_of_capitals = zeros(num_lambdas)
ratio_of_outputs  = zeros(num_lambdas)

old_TFPis = zeros(num_lambdas)
old_TFPds = zeros(num_lambdas)
Labour_productivity_s = zeros(num_lambdas)
Capital_productivity_s = zeros(num_lambdas)
Investment_productivity_s = zeros(num_lambdas)

# Distribution of output to different channels of earnigs/income to occupations
share_W_earnings_in_output = zeros(num_lambdas)
share_SP_earnings_in_output = zeros(num_lambdas)
share_EMP_earnings_in_output = zeros(num_lambdas)
share_W_capital_income_in_output = zeros(num_lambdas)
share_SP_capital_income_in_output = zeros(num_lambdas)
share_EMP_capital_income_in_output = zeros(num_lambdas)
#throw(error)

@showprogress for i = 1:num_lambdas
    #println("\n$(country) - $(i)/$(num_lambdas)")
    # ss_star = [res, r, w, approx_object, params]
    ss_star = copy(SSS[i])

    Rs[i] = ss_star[2]
    Ws[i] = ss_star[3]
    K_Ys[i] = ss_star[1][12]
    C_Ys[i] = ss_star[1][13]
    occWs[i] = ss_star[1][14]
    occSEs[i] = ss_star[1][15]
    occEMPs[i] = ss_star[1][16]
    occENTs[i] = occSEs[i] + occEMPs[i]
    logcs[i] = ss_star[1][17]
    loges[i] = ss_star[1][19]
    giniWs[i] = ss_star[1][20]
    giniEnts[i] = ss_star[1][21]

    eta = ss_star[5][5]
    theta = ss_star[5][6]

    z_m_nodes = ss_star[4][1]
    z_w_nodes = ss_star[4][2]

    asset_grid = ss_star[1][3]
    density_distr = ss_star[1][5]
    wealth = ones(size(density_distr)).*asset_grid
    policy = ss_star[1][4]
    income = ss_star[1][23] .- wealth
    earnings = ss_star[1][24]
    consumption = income .+ wealth .- policy
    occ_choice = ss_star[1][22]
    capital_d = ss_star[1][26]
    labour_d = ss_star[1][29]

    output = ss_star[1][32]

    Capital[1,i] = sum(asset_grid.*density_distr)

    Outputs[i] = sum(output .* density_distr)
    Incomes[i] = sum(income .* density_distr)
    Consumptions[i] = sum(consumption .* density_distr)#Incomes[i] - ss_star[1][7]
    Income_to_outputs[i] = Incomes[i]/Outputs[i]

    number_u_nodes = SSS[i][4][7]
    number_asset_grid = ss_star[1][2]
    number_zeta_nodes = global_approx_params[4]
    number_alpha_m_nodes = global_approx_params[5]
    number_alpha_w_nodes = global_approx_params[6]

    w_choice   = copy(occ_choice)
    sp_choice  = copy(occ_choice)
    emp_choice = copy(occ_choice)
    Threads.@threads for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
        for a_i in 1:number_asset_grid
            if w_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] != 1
                w_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 0
            end
            if sp_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] != 2
                sp_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 0
            else
                sp_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 1
            end
            if emp_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] != 3
                emp_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 0
            else
                emp_choice[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i] = 1
            end
        end
    end
    ent_choice = sp_choice.+emp_choice

    Capital[2,i] = sum(asset_grid.* (density_distr.*w_choice)./sum(density_distr.*w_choice) )
    Capital[3,i] = sum(asset_grid.* (density_distr.*sp_choice)./sum(density_distr.*sp_choice) )
    Capital[4,i] = sum(asset_grid.* (density_distr.*emp_choice)./sum(density_distr.*emp_choice) )

    # [1] = income, earnings, wealth, consumption
    # [2] = All, W, SP, EMP
    # [3] = lambdas
    # = means, ginis, varlogs
    for s = 1:4 # [1] = income, earnings, wealth, consumption
        # if s == 1
        stat_distr = copy(income)
        if s == 2
            stat_distr = copy(earnings)
        elseif s == 3
            stat_distr = copy(wealth)
        elseif s == 4
            stat_distr = copy(consumption)
        end

        for h = 1:5 # [2] = All, W, SP, EMP, ENT
            choice = copy(occ_choice)
            if h == 1 #All
                choice = choice.*0.0 .+ 1.0
            elseif h == 2 #W
                choice = copy(w_choice)
            elseif h == 3 #SP
                choice = copy(sp_choice)
            elseif h == 4 #EMP
                choice = copy(emp_choice)
            elseif h == 5 #ENT
                choice = copy(ent_choice)
            end

            # calculate mean
            means[s,h,i] = sum(stat_distr .* (density_distr.*choice)/sum(density_distr.*choice) )

            # calculate gini coefficent
            density_distr_choice = (density_distr.*choice)./sum(density_distr.*choice)
            stat_choice = stat_distr.*choice
            density_distr_choice_vec = vec(density_distr_choice)
            stat_choice_vec = vec(stat_choice)

            index_non_zero = findall(x-> x>1e-5,density_distr_choice_vec)
            density_distr_choice_vec_non_zero = density_distr_choice_vec[index_non_zero]
            stat_choice_vec_non_zero = stat_choice_vec[index_non_zero]

            ginis[s,h,i] = sum([ density_distr_choice_vec_non_zero[y_i]*density_distr_choice_vec_non_zero[y_j]*abs(stat_choice_vec_non_zero[y_i]-stat_choice_vec_non_zero[y_j]) for y_i=1:length(density_distr_choice_vec_non_zero), y_j=1:length(density_distr_choice_vec_non_zero) ]) / (2*sum(stat_choice_vec_non_zero.*density_distr_choice_vec_non_zero))

            # calculate variance of log-s
            # avglogs[s,h,i] = sum(density_distr.*choice.*log.(max.(1e-12,stat_distr))   )
            # varlogs[s,h,i] = sum(density_distr.*choice.*(log.(max.(1e-12,stat_distr)).-avglogs[s,h,i]).^2)/sum(density_distr.*choice)
            # avglogs[s,h,i] /= sum(density_distr.*choice)
            avglogs[s,h,i] = sum( (density_distr.*choice./sum(density_distr.*choice)) .* log.(max.(1e-12,stat_distr)) )
            varlogs[s,h,i] = sum( (density_distr.*choice./sum(density_distr.*choice)) .* (log.(max.(1e-12,stat_distr)) .- avglogs[s,h,i]).^2 )
            if s==3
                # avglogs[s,h,i] = sum(density_distr.*choice.*max.(1e-12,stat_distr)   )
                # varlogs[s,h,i] = sum(density_distr.*choice.*(max.(1e-12,stat_distr).- avglogs[s,h,i]).^2)/sum(density_distr.*choice)
                # avglogs[s,h,i] /= sum(density_distr.*choice)
                avglogs[s,h,i] = sum( (density_distr.*choice./sum(density_distr.*choice)) .* max.(1e-12,stat_distr) )
                varlogs[s,h,i] = sum( (density_distr.*choice./sum(density_distr.*choice)) .* (max.(1e-12,stat_distr) .- avglogs[s,h,i]).^2 )
            end

            # avgs[s,h,i] = sum(density_distr.*choice.*max.(1e-12,stat_distr)   )
            # vars[s,h,i] = sum(density_distr.*choice.*(max.(1e-12,stat_distr).- avgs[s,h,i]).^2)/sum(density_distr.*choice)
            # avgs[s,h,i] /= sum(density_distr.*choice)
            avgs[s,h,i] = sum( (density_distr.*choice./sum(density_distr.*choice)) .* max.(1e-12,stat_distr) )
            vars[s,h,i] = sum( (density_distr.*choice./sum(density_distr.*choice)) .* (max.(1e-12,stat_distr) .- avgs[s,h,i]).^2 )

            if h!=5
                try
                    quantile_means[:,h,s,i] .= quantile_mean(stat_distr, density_distr.*choice)
                catch e
                    quantile_means[:,h,s,i] .= NaN
                end
            end
        end
    end

    ratio_of_capitals[i] = ss_star[1][6]

    unlimited_capital_choice = capital_d.<(lambdas[i].*(asset_grid.*ones(size(capital_d))))
    for c = 1:4
        if c==1
            choice = w_choice
        elseif c==2
            choice = ent_choice
        elseif c==3
            choice = sp_choice
        elseif c==4
            choice = emp_choice
        else
            throw("c is out of bounds")
        end

        denumerator = sum(choice.*density_distr)
        if c>1
            #if c==2
            #    choice = ones(size(choice))
            #end
            #if c==2
            #    denumerator = sum(choice.*density_distr)
            #end
            Y = sum( (choice.*density_distr./denumerator) .* output)
            K = sum( (choice.*density_distr./denumerator) .* capital_d)
            L = sum( (choice.*density_distr./denumerator) .* labour_d)
            TFPis[c-1,i] = Y/(K^eta*L^theta)
            TFPds[c-1,i] = Y/K^(eta/(theta+eta))

            c_1 = capital_d.^(-1.0)
            replace!(c_1,Inf=>0.0)
            # mean_MPK[c-1,i] = sum(choice.*density_distr.*(eta.*output.*c_1))/denumerator
            # var_MPK[c-1,i] = sum(choice.*density_distr.*(eta.*output.*c_1 .- mean_MPK[c-1,i]*denumerator).^2)/denumerator
            mean_MPK[c-1,i] = sum( (choice.*density_distr./denumerator) .* (eta.*output.*c_1) )
            var_MPK[c-1,i] = sum( (choice.*density_distr./denumerator) .* (eta.*output.*c_1 .- mean_MPK[c-1,i]).^2 )

            l_1 = labour_d.^(-1.0)
            replace!(l_1,Inf=>0.0)
            # mean_MPL[c-1,i] = sum(choice.*density_distr.*(theta.*output.*l_1))/denumerator
            # var_MPL[c-1,i] = sum(choice.*density_distr.*(theta.*output.*l_1 .- mean_MPL[c-1,i]*denumerator).^2)/denumerator
            mean_MPL[c-1,i] = sum( (choice.*density_distr./denumerator) .* (theta.*output.*l_1) )
            var_MPL[c-1,i] = sum( (choice.*density_distr./denumerator) .* (theta.*output.*l_1 .- mean_MPL[c-1,i]).^2 )

            #if c==2
            #    choice = ent_choice
            #end
            share_unbound[c-1,i] = sum( (choice.*density_distr./sum(choice.*density_distr)) .* unlimited_capital_choice)
        end

        #[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
        #m = u_i,alpha_m_i,zeta_i
        m_skill_distr = permutedims(sum(density_distr.*choice, dims=[1,5])[1,:,:,:,1], [1,3,2])
        #w = u_i,alpha_m_i,alpha_w_i
        w_skill_distr = sum(density_distr.*choice, dims=[1,3])[1,:,1,:,:]

        # avg_m_skill[c,i] = sum(m_skill_distr.*z_m_nodes)/denumerator
        # var_m_skill[c,i] = sum(m_skill_distr.*(z_m_nodes .- avg_m_skill[c,i]*denumerator).^2)/denumerator
        avg_m_skill[c,i] = sum( (m_skill_distr./denumerator) .* z_m_nodes )
        var_m_skill[c,i] = sum( (m_skill_distr./denumerator) .* (z_m_nodes .- avg_m_skill[c,i]).^2 )
        sum_m_skill[c,i] = sum(m_skill_distr.*z_m_nodes)

        # avg_w_skill[c,i] = sum(w_skill_distr.*z_w_nodes)/denumerator
        # var_w_skill[c,i] = sum(w_skill_distr.*(z_w_nodes .- avg_w_skill[c,i]*denumerator).^2)/denumerator
        avg_w_skill[c,i] = sum( (w_skill_distr./denumerator) .* z_w_nodes )
        var_w_skill[c,i] = sum( (w_skill_distr./denumerator) .* (z_w_nodes .- avg_w_skill[c,i]).^2 )
        sum_w_skill[c,i] = sum(w_skill_distr.*z_w_nodes)

        #display([avg_w_skill[c,i],avg_m_skill[c,i], var_w_skill[c,i],var_m_skill[c,i]])
    end

    old_TFPis[i] = ss_star[1][43][end-4]
    old_TFPds[i] = ss_star[1][43][end-3]
    Labour_productivity_s[i] = ss_star[1][43][end-2]
    Capital_productivity_s[i] = ss_star[1][43][end-1]
    Investment_productivity_s[i] = ss_star[1][43][end]

    agg_c_e = ss_star[5][7]*sum(density_distr.*emp_choice)
    share_W_earnings_in_output[i] = sum(ss_star[1][24] .* density_distr.*w_choice)/(Outputs[i]-agg_c_e)
    share_SP_earnings_in_output[i] = sum(ss_star[1][24] .* density_distr.*sp_choice)/(Outputs[i]-agg_c_e)
    share_EMP_earnings_in_output[i] = sum(ss_star[1][24] .* density_distr.*emp_choice)/(Outputs[i]-agg_c_e)
    share_W_capital_income_in_output[i] = (ss_star[5][3]+ss_star[1][44])*sum(ones(size(ss_star[1][5])).*ss_star[1][3] .* density_distr.*w_choice)/(Outputs[i]-agg_c_e)
    share_SP_capital_income_in_output[i] = (ss_star[5][3]+ss_star[1][44])*sum(ones(size(ss_star[1][5])).*ss_star[1][3] .* density_distr.*sp_choice)/(Outputs[i]-agg_c_e)
    share_EMP_capital_income_in_output[i] = (ss_star[5][3]+ss_star[1][44])*sum(ones(size(ss_star[1][5])).*ss_star[1][3] .* density_distr.*emp_choice)/(Outputs[i]-agg_c_e)
    #throw(error)
end

asset_grid = SSS[1][1][3]
number_non_zero_asset_grid = 0
number_u_nodes = SSS[1][4][7]
number_zeta_nodes = global_approx_params[4]
number_alpha_m_nodes = global_approx_params[5]
number_alpha_w_nodes = global_approx_params[6]
for (u_i,(zeta_i,(alpha_m_i,alpha_w_i))) in collect(Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes))))
    global number_non_zero_asset_grid
    temp = findlast(round.(SSS[1][1][5][:,u_i,zeta_i,alpha_m_i,alpha_w_i]; digits=7) .!= 0.0)
    if isnothing(temp)
        temp = 0
    end
    if temp > number_non_zero_asset_grid
        number_non_zero_asset_grid = temp
    end
end
#policy = SSS[l][1][4]
num_lambdas_sp = findfirst(x -> x>0.75, C_Ys)-1

throw(error)

function create_plot(X::Vector{Float64},XLABEL::String,Y::Vector{Float64},YLABEL::String, IS_Y_PERCENTAGE::Bool=true, TICKSFONTSIZE::Int64=9, NUM_YTICKS::Int64=10, OCCUPATION::String="H", LOW_LIMIT=-Inf)
    COLOR="blue"
    if OCCUPATION=="W"
        COLOR="purple"
    elseif OCCUPATION=="SP"
        COLOR="red"
    elseif OCCUPATION=="EMP"
        COLOR="green"
    end

    YLIMS1 = max(minimum(Y), LOW_LIMIT)
    YLIMS2 = maximum(Y)
    #YLIMS=(YLIMS1-0.01, YLIMS2+0.01)
    YLIMMARGIN = abs(YLIMS2-YLIMS1)*0.015
    YLIMS=(YLIMS1-YLIMMARGIN, YLIMS2+YLIMMARGIN)
    YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=NUM_YTICKS))
    # YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=Int(round(num_lambdas/2;digits=0))))


    DIGITS = Int(max(2, round(log10(1/(YTICKS[2]-YTICKS[1]));digits=0) ))
    YPOS = mean(YTICKS[end-1:end])
    if maximum(Y[calibrated_lambda:calibrated_lambda+4]) > YTICKS[end-1]
        YPOS = mean(YTICKS[1:2])
    end
    if IS_Y_PERCENTAGE
        YTICKS = (YTICKS, ["$(round(100*y;digits=DIGITS))%" for y in YTICKS])
    else
        YTICKS = (YTICKS, [round(y;digits=DIGITS) for y in YTICKS])
    end
    plt = scatter(X,Y,
                    #=regression=true,=#
                    color=COLOR,
                    legend=false,
                    xlabel=XLABEL,
                    # ylabel=YLABEL,
                    xtickfontsize=TICKSFONTSIZE,
                    ytickfontsize=TICKSFONTSIZE,
                    xguidefontsize=TICKSFONTSIZE,
                    yguidefontsize=TICKSFONTSIZE,
                    yticks = YTICKS,
                    ylims = YLIMS )
    vline!([X[calibrated_lambda]], color="grey")
    if IS_Y_PERCENTAGE
        TEXT = " Calibrated economy ($(round(X[calibrated_lambda]; digits=DIGITS+1)),$(round(Y[calibrated_lambda]*100; digits=DIGITS+1))%) "
    else
        TEXT = " Calibrated economy ($(round(X[calibrated_lambda]; digits=DIGITS+1)),$(round(Y[calibrated_lambda]; digits=DIGITS+1))) "
    end

    annotate!([X[calibrated_lambda]], YPOS, text(TEXT, :grey, :left, 7))
    return plt
end
function create_plots(X::Vector{Float64},XLABEL::String,Ys,YLABELs, IS_Y_PERCENTAGE::Bool=true, TICKSFONTSIZE::Int64=9, NUM_YTICKS::Int64=10, OCCUPATION=["H", "W", "SP", "EMP"], LOW_LIMIT=-Inf)
    plts = []
    half_range = maximum([ maximum(filter(!isnan,copy(Ys[i])))-minimum(filter(!isnan,copy(Ys[i]))) for i in 1:length(Ys)])/2.0
    for y_i = 1:length(Ys)
        Y = Ys[y_i]
        Ynanless = copy(Y)
        Ynanless = filter(!isnan,Ynanless)
        ns = num_lambdas - (length(Y) - length(Ynanless))
        YLABEL = YLABELs[y_i]

        YLIMS1 = max(middle(Ynanless[1:ns])-half_range, LOW_LIMIT)
        if IS_Y_PERCENTAGE
            YLIMS1 = max(YLIMS1, min(0.0, minimum(filter(!isnan,copy(Y)))))
        end
        YLIMS2 = YLIMS1+2*half_range
        YLIMMARGIN = abs(YLIMS2-YLIMS1)*0.015
        YLIMS=(YLIMS1-YLIMMARGIN, YLIMS2+YLIMMARGIN)
        YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=NUM_YTICKS))
        # YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=Int(round(num_lambdas/2;digits=0))))

        DIGITS = Int(max(2, round(log10(1/(YTICKS[2]-YTICKS[1]));digits=0) ))+1
        YPOS = mean(YTICKS[end-1:end])
        if maximum(filter(!isnan,copy(Y[calibrated_lambda:calibrated_lambda+4]))) > YTICKS[end-1]
            YPOS = mean(YTICKS[1:2])
        end
        if IS_Y_PERCENTAGE
            YTICKS = (YTICKS, ["$(round(100*y;digits=DIGITS))%" for y in YTICKS])
        else
            YTICKS = (YTICKS, [round(y;digits=DIGITS) for y in YTICKS])
        end

        COLOR="blue"
        if OCCUPATION[y_i]=="W"
            COLOR="purple"
        elseif OCCUPATION[y_i]=="SP"
            COLOR="red"
        elseif OCCUPATION[y_i]=="EMP"
            COLOR="green"
        end
        plt = scatter(X,Y,
                        #=regression=true,=#
                        color = COLOR,
                        legend = false,
                        legendfontsize=ceil(Int64,TICKSFONTSIZE*0.72),
                        xlabel = XLABEL,
                        # ylabel = YLABEL,
                        xtickfontsize=TICKSFONTSIZE,
                        ytickfontsize=TICKSFONTSIZE,
                        xguidefontsize=TICKSFONTSIZE,
                        yguidefontsize=TICKSFONTSIZE,
                        yticks = YTICKS,
                        ylims = YLIMS )
        vline!([X[calibrated_lambda]], color="grey")
        if IS_Y_PERCENTAGE
            TEXT = " Calibrated economy ($(round(X[calibrated_lambda]; digits=DIGITS+1)),$(round(Y[calibrated_lambda]*100; digits=DIGITS+1))%) "
        else
            TEXT = " Calibrated economy ($(round(X[calibrated_lambda]; digits=DIGITS+1)),$(round(Y[calibrated_lambda]; digits=DIGITS+1))) "
        end
        annotate!([X[calibrated_lambda]], YPOS, text(TEXT, :grey, :left, 7))
        push!(plts, plt)
    end
    return plts
end
function create_combined_plot(X::Vector{Float64},XLABEL::String,Ys,YLABELs,YLABEL, IS_Y_PERCENTAGE::Bool=true, TICKSFONTSIZE::Int64=9, NUM_YTICKS::Int64=10, OCCUPATION=["H", "W", "SP", "EMP"], LOW_LIMIT=-Inf, LEGEND=:outertopright)


    YLIMS1 = max(minimum(minimum.(Ys)), LOW_LIMIT)
    YLIMS2 = maximum(maximum.(Ys))
    #YLIMS=(YLIMS1-0.01, YLIMS2+0.01)
    YLIMMARGIN = abs(YLIMS2-YLIMS1)*0.015
    YLIMS=(YLIMS1-YLIMMARGIN, YLIMS2+YLIMMARGIN)
    YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=NUM_YTICKS))
    # YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=Int(round(num_lambdas/2;digits=0))))

    DIGITS = Int(max(2, round(log10(1/(YTICKS[2]-YTICKS[1]));digits=0) ))
    YPOS = mean(YTICKS[end-1:end])

    COLORS=[]
    for y_i = 1:length(Ys)
        if OCCUPATION[y_i]=="W"
            push!(COLORS,"purple")
        elseif OCCUPATION[y_i]=="SP"
            push!(COLORS,"red")
        elseif OCCUPATION[y_i]=="EMP"
            push!(COLORS,"green")
        elseif OCCUPATION[y_i]=="ENT"
            push!(COLORS,"yellow")
        else
            push!(COLORS,"blue")
        end

        if maximum(Ys[y_i][calibrated_lambda:calibrated_lambda+4]) > YTICKS[end-1]
            YPOS = mean(YTICKS[1:2])
        end
    end
    if IS_Y_PERCENTAGE
        YTICKS = (YTICKS, ["$(round(100*y;digits=DIGITS))%" for y in YTICKS])
    else
        YTICKS = (YTICKS, [round(y;digits=DIGITS) for y in YTICKS])
    end
    plt = vline([X[calibrated_lambda]], color="grey",label="")
    if IS_Y_PERCENTAGE
        TEXT = " Calibrated economy ($(round(X[calibrated_lambda]; digits=DIGITS+1))) "
    else
        TEXT = " Calibrated economy ($(round(X[calibrated_lambda]; digits=DIGITS+1))) "
    end
    annotate!([X[calibrated_lambda]], YPOS, text(TEXT, :grey, :left, 7))

    for y_i = 1:length(Ys)
        plt = scatter!(X,Ys[y_i],
                        #=regression=true,=#
                        color=COLORS[y_i],
                        legend=LEGEND,
                        legendfontsize=ceil(Int64,TICKSFONTSIZE*0.72),
                        xlabel = XLABEL,
                        # ylabel = YLABEL,
                        xtickfontsize=TICKSFONTSIZE,
                        ytickfontsize=TICKSFONTSIZE,
                        xguidefontsize=TICKSFONTSIZE,
                        yguidefontsize=TICKSFONTSIZE,
                        label=YLABELs[y_i],
                        yticks = YTICKS,
                        ylims = YLIMS )
    end

    return plt
end

LOCAL_DIR_GENERAL = "$(LOCAL_DIR)/General/"
if Sys.iswindows()
    LOCAL_DIR_GENERAL = "$(LOCAL_DIR)\\General\\"
end
mkpath(LOCAL_DIR_GENERAL)
#Interest rate and wage
plt = create_plot(C_Ys,"Credit/Output", Rs,"Interest Rate", true, 12,7)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_credit_to_gdp_interest_rate.png")
plt = create_plot(C_Ys,"Credit/Output", Ws,"Wage", false, 12,7)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_credit_to_gdp_interest_wage.png")

#lambda to: Output, Consumption, Income/Output, Income
plt = create_plot(lambdas,"λ", Outputs,"Output", false, 12,7)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Outputs.png")
plt = create_plot(lambdas,"λ", Consumptions,"Consumption", false, 12,7)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Consumptions.png")
plt = create_plot(lambdas,"λ", Income_to_outputs,"Income/Output", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Income_to_outputs.png")
plt = create_plot(lambdas,"λ", Incomes,"Income", false, 12,7)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Incomes.png")

plt = create_plot(lambdas,"λ", Capital[1,:],"Capital", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Capitals.png")

plt = create_plot(lambdas,"λ", Capital[2,:],"Workers' Capital", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Capitals_W.png")
plt = create_plot(lambdas,"λ", Capital[3,:],"Sole Proprietors' Capital", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Capitals_SP.png")
plt = create_plot(lambdas,"λ", Capital[4,:],"Employers' Capital", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Capitals_EMP.png")

#Credit/Output to Income/Output
plt = create_plot(C_Ys,"Credit/Output", Income_to_outputs,"Income/Output", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_credit_to_gdp_Income_to_outputs.png")

#lambda to Credit/Output, Credit/Output to Capital/Output
plt = create_plot(lambdas,"λ", C_Ys,"Credit/Output", false, 12,7)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_credit_to_gdp.png")
plt = create_plot(C_Ys,"Credit/Output", K_Ys,"Capital/Output", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_credit_to_gdp_capital_to_gdp.png")

LOCAL_DIR_OCCUPATION = "$(LOCAL_DIR)/Occupation/"
if Sys.iswindows()
    LOCAL_DIR_OCCUPATION = "$(LOCAL_DIR)\\Occupation\\"
end
mkpath(LOCAL_DIR_OCCUPATION)
#Crdit/Output to share of W, SP, EMP
plts = create_plots(C_Ys,"Credit/Output", [occENTs, occSEs, occEMPs, occWs],["Share of Entrepreneurs","Share of Sole Proprietors","Share of Employers","Share of Workers"],true, 12,7,["ENT","SP","EMP","W"])
display(plts[1])
savefig(plts[1],"$(LOCAL_DIR_OCCUPATION)$(country)_credit_to_gdp_entrepreneurs.png")
display(plts[2])
savefig(plts[2],"$(LOCAL_DIR_OCCUPATION)$(country)_credit_to_gdp_sole_proprietors.png")
display(plts[3])
savefig(plts[3],"$(LOCAL_DIR_OCCUPATION)$(country)_credit_to_gdp_employers.png")
display(plts[4])
savefig(plts[4],"$(LOCAL_DIR_OCCUPATION)$(country)_credit_to_gdp_workers.png")

LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)/Inequality/"
if Sys.iswindows()
    LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)\\Inequality\\"
end
mkpath(LOCAL_DIR_INEQUALITY)
#Credit/Output to: log-consumption, log-earnings
plt = create_plot(C_Ys,"Credit/Output", logcs,"Variance of log-consumption", false, 12,7)
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_log_consumption.png")

plt = create_plot(C_Ys,"Credit/Output", loges,"Variance of log-earnings", false, 12,7)
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_log_earnings.png")

#Crdit/Output to Gini Workers and Entrepreneurs
plts = create_plots(C_Ys,"Credit/Output", [giniWs, giniEnts],["Gini for workers' income","Gini for entrepreneurs' income"],false, 12,7,["W","ENT"])
for i in 1:length(plts)
    display(plts[i])
    if i==1
        savefig(plts[i],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_gini_workers.png")
    elseif i==2
        savefig(plts[i],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_gini_entrepreneurs.png")
    end
end

# households' policy a_t+1
middle_point = Int64(round((1+num_lambdas)/2;digits=0))
middle_point1 = Int64(round((1+middle_point)/2;digits=0))
middle_point2 = Int64(round((middle_point+num_lambdas)/2;digits=0))
suitable_ls = [1,middle_point1,middle_point,middle_point2,num_lambdas]
#,u$(u),z$(zeta),aw$(alpha_w)
labels = permutedims(["λ=$(round(lambdas[l];digits=2)),am$(alpha_m)" for l in suitable_ls for u in [Int64(round((1+number_u_nodes)/2;digits=0))] for zeta in [Int64(round((1+number_zeta_nodes)/2;digits=0))] for alpha_m in [#=Int64(round((1+number_alpha_m_nodes)/2;digits=0))=#1,number_alpha_m_nodes] for alpha_w in [Int64(round((1+number_alpha_w_nodes)/2;digits=0))]])
p2 = plot(asset_grid[1:number_non_zero_asset_grid],[SSS[l][1][4][1:number_non_zero_asset_grid,u,zeta,alpha_m,alpha_w] for l in suitable_ls for u in [Int64(round((1+number_u_nodes)/2;digits=0))] for zeta in [Int64(round((1+number_zeta_nodes)/2;digits=0))] for alpha_m in [#=Int64(round((1+number_alpha_m_nodes)/2;digits=0))=#1,number_alpha_m_nodes] for alpha_w in [Int64(round((1+number_alpha_w_nodes)/2;digits=0))]], label=labels, legend=:outertopleft)
plt = plot(p2, title="Policy functions (high vs low fixed effects)")
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_policy_hvsl_fixed_lambda.png")

for s = 1:4 # [1] = income, earnings, wealth, consumption
    # if s == 1
    stat_name = "Income"
    if s == 2
        stat_name = "Earnings"
    elseif s == 3
        stat_name = "Wealth"
    elseif s == 4
        stat_name = "Consumption"
    end

    CHOICE_NAMES = ["Households","Workers","Sole Proprietors","Employers","Entrepreneurs"]
    OCCS = ["H","W","SP","EMP","ENT"]

    LABELS=["1st","2nd","3rd","4th","5th"]
    for h=1:5
        LABEL = "Mean of $(CHOICE_NAMES[h])' $(stat_name) (quantiles)"
        LEGENDPOS = false
        if h==4
            LEGENDPOS = :topright
        end
        try
            plt = create_combined_plot(C_Ys,"Credit/Output", [quantile_means[1,h,s,:],quantile_means[2,h,s,:],quantile_means[3,h,s,:],quantile_means[4,h,s,:],quantile_means[5,h,s,:]],LABELS,LABEL, false, 12,7, ["H","W","SP","EMP","ENT"], -Inf, LEGENDPOS)
            display(plt)
            savefig(plt,"$(LOCAL_DIR_INEQUALITY)time_combined_mean_$(stat_name)_$(OCCS[h])_quantiles.png")
        catch e

        end
    end
    # calculate mean
    #means[s,h,i]
    LABELS = ["Mean of $cn' $stat_name" for cn in CHOICE_NAMES]
    plts = create_plots(C_Ys,"Credit/Output", [means[s,1,:], means[s,2,:], means[s,3,:], means[s,4,:], means[s,5,:]],LABELS,false, 12,7,["H","W","SP","EMP","ENT"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_mean_$(choice_name)_$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", means[s,h,:],LABELS[h],false, 12,7,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_mean_$(choice_name)_$(stat_name).png")
    end

    # calculate gini coefficent
    #ginis[s,h,i]
    LABELS = ["Gini of $cn' $stat_name" for cn in CHOICE_NAMES]
    plts = create_plots(C_Ys,"Credit/Output", [ginis[s,1,:], ginis[s,2,:], ginis[s,3,:], ginis[s,4,:], ginis[s,5,:]],LABELS,false, 12,7,["H","W","SP","EMP","ENT"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_gini_$(choice_name)_$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", ginis[s,h,:],LABELS[h],false, 12,7,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_gini_$(choice_name)_$(stat_name).png")
    end
    # calculate variance of log-s
    #avgs[s,h,i]
    LABELS = ["Average of $cn' $stat_name" for cn in CHOICE_NAMES]
    plts = create_plots(C_Ys,"Credit/Output", [avgs[s,1,:], avgs[s,2,:], avgs[s,3,:], avgs[s,4,:], avgs[s,5,:]],LABELS,false, 12,7,["H","W","SP","EMP","ENT"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", avgs[s,h,:],LABELS[h],false, 12,7,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_$(stat_name).png")
    end

    #vars[s,h,i]
    LABELS = ["Variance of $cn' $stat_name" for cn in CHOICE_NAMES]
    plts = create_plots(C_Ys,"Credit/Output", [vars[s,1,:], vars[s,2,:], vars[s,3,:], vars[s,4,:], vars[s,5,:]],LABELS,false, 12,7,["H","W","SP","EMP","ENT"], 0.0)
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", vars[s,h,:],LABELS[h],false, 12,7,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_$(stat_name).png")
    end

    # calculate variance of log-s
    #avglogs[s,h,i]
    LABELS = ["Average of $cn' Log-$stat_name" for cn in CHOICE_NAMES]
    plts = create_plots(C_Ys,"Credit/Output", [avglogs[s,1,:], avglogs[s,2,:], avglogs[s,3,:], avglogs[s,4,:], avglogs[s,5,:]],LABELS,false, 12,7,["H","W","SP","EMP","ENT"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_log$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", avglogs[s,h,:],LABELS[h],false, 12,7,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_log$(stat_name).png")
    end

    #varlogs[s,h,i]
    if s!=3
        LABELS = ["Variance of $cn' Log-$stat_name" for cn in CHOICE_NAMES]
    else
        LABELS = ["Variance of $cn' $stat_name" for cn in CHOICE_NAMES]
    end
    plts = create_plots(C_Ys,"Credit/Output", [varlogs[s,1,:], varlogs[s,2,:], varlogs[s,3,:], varlogs[s,4,:], varlogs[s,5,:]],LABELS,false, 12,7,["H","W","SP","EMP","ENT"], 0.0)
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_log$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", varlogs[s,h,:],LABELS[h],false, 12,7,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_log$(stat_name).png")
    end
end
#=
# output by SP
plt = create_plot(C_Ys,"Credit/Output", output_SPs,"Aggregate output by Sole Proprietors", false, "SP")
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_output_sp.png")
# output by EMP
plt = create_plot(C_Ys,"Credit/Output", output_EMPs,"Aggregate output by Employers", false, "EMP")
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_output_emp.png")

# share of output by SP
plt = create_plot(C_Ys,"Credit/Output", share_of_output_SPs,"Share of output by Sole Proprietors", true, "SP")
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_share_of_output_sp.png")
# share of output by EMP
plt = create_plot(C_Ys,"Credit/Output", share_of_output_EMPs,"Share of output by Employers", true, "EMP")
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_share_of_output_emp.png")
=#

ratio_of_capitals = ratio_of_capitals./ratio_of_capitals[num_lambdas]
ratio_of_outputs = Outputs./Outputs[num_lambdas]
# ratio of capital restricted to unrestricted
plt = create_plot(C_Ys,"Credit/Output", ratio_of_capitals,"Ratio of Restricted Capital to Unrestricted", true)
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_ratio_of_capital.png")
# ratio of output restricted to unrestricted
plt = create_plot(C_Ys,"Credit/Output", ratio_of_outputs,"Ratio of Restricted Output to Unrestricted", true)
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_ratio_of_output.png")

# aggregate output ~~~ DOESN'T WORK
#=
using StatsPlots
xt = Plots.xticks(plt[1])
plt2 = groupedbar(C_Ys[1:num_lambdas],[output_SPs[1:num_lambdas] output_EMPs[1:num_lambdas]],
        bar_position = :stack,
        #bar_width=[abs(C_Ys[2]-C_Ys[1]); [abs(C_Ys[i+1]-C_Ys[i-1]) for i=2:num_lambdas-1]; abs(C_Ys[num_lambdas]-C_Ys[num_lambdas-1])]./2.0,
        bar_width=median([abs(C_Ys[i+1]-C_Ys[i-1]) for i=2:num_lambdas])./2.0,
        xticks=xt,
        label=["SP" "EMP"],
        color=["red" "green"])
scatter!(C_Ys[1:num_lambdas],[Outputs[1:num_lambdas] output_EMPs[1:num_lambdas]],#=regression=true,=#color=["blue" "green"],xlabel="Credit/Output",ylabel="Output",label=["Total" "EMP"],legend=:bottomright)
vline!([C_Ys[calibrated_lambda]], color="grey")
annotate!([C_Ys[calibrated_lambda]], Outputs[num_lambdas], text(" Calibrated economy ($(round(C_Ys[calibrated_lambda]; digits=2)),$(round(Outputs[num_lambdas]; digits=2))) ", :grey, :left, 7))
display(plt2)
savefig(plt2,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_Outputs_modified.png")
=#

LOCAL_DIR_PRODUCTIVITY = "$(LOCAL_DIR)/Productivity/"
if Sys.iswindows()
    LOCAL_DIR_PRODUCTIVITY = "$(LOCAL_DIR)\\Productivity\\"
end
mkpath(LOCAL_DIR_PRODUCTIVITY)
# TFP_ideal for ENT
plt = create_plot(C_Ys,"Credit/Output",TFPis[1,:],"TFP for Entrepreneurs",false, 12,7,"ENT")
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_tfp_ideal_ent.png")
# TFP_ideal for SP,EMP
plts = create_plots(C_Ys,"Credit/Output", [TFPis[2,:], TFPis[3,:]],["TFP for Sole Proprietors","TFP for Employers"],false, 12,7,["SP","EMP"])
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_tfp_ideal_sp.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_tfp_ideal_emp.png")
# TFP_data for SP,EMP
plts = create_plots(C_Ys,"Credit/Output", [TFPds[2,:], TFPds[3,:]],["TFP_data for Sole Proprietors","TFP_data for Employers"],false, 12,7,["SP","EMP"])
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_tfp_data_sp.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_tfp_data_emp.png")


# mean of MPL for ENT, SP,EMP
plt = create_plot(C_Ys,"Credit/Output",mean_MPL[1,:],"Mean of MPL for Entrepreneurs", false, 12,7, "ENT")
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_mean_mpl_ent.png")
plts = create_plots(C_Ys,"Credit/Output", [mean_MPL[2,:], mean_MPL[3,:]],["Mean of MPL for Sole Proprietors","Mean of MPL for Employers"],false, 12,7,["SP","EMP"])
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_mean_mpl_sp.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_mean_mpl_emp.png")
LABELS = ["SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [mean_MPL[2,:], mean_MPL[3,:]],LABELS,"Mean of MPL",false, 12,7,["SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_credit_to_gdp_mean_mpl.png")

# variance of MPL for ENT, SP,EMP
plt = create_plot(C_Ys,"Credit/Output",var_MPL[1,:],"Variance of MPL for Entrepreneurs", false, 12,7, "ENT", 0.0)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_var_mpl_ent.png")
plts = create_plots(C_Ys,"Credit/Output", [var_MPL[2,:], var_MPL[3,:]],["Variance of MPL for Sole Proprietors","Variance of MPL for Employers"],false, 12,7,["SP","EMP"],0.0)
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_var_mpl_sp.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_var_mpl_emp.png")
LABELS = ["SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [var_MPL[2,:], var_MPL[3,:]],LABELS,"Variance of MPL",false, 12,7,["SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_credit_to_gdp_var_mpl.png")

# mean of MPK for ENT, SP,EMP
plt = create_plot(C_Ys,"Credit/Output",mean_MPK[1,:],"Mean of MPK for Entrepreneurs", false, 12,7, "ENT")
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_mean_mpk_ent.png")
plts = create_plots(C_Ys,"Credit/Output", [mean_MPK[2,:], mean_MPK[3,:]],["Mean of MPK for Sole Proprietors","Mean of MPK for Employers"],false, 12,7,["SP","EMP"])
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_mean_mpk_sp.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_mean_mpk_emp.png")
LABELS = ["SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [mean_MPK[2,:], mean_MPK[3,:]],LABELS,"Mean of MPK",false, 12,7,["SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_credit_to_gdp_mean_mpk.png")

# variance of MPK for ENT, SP,EMP
plt = create_plot(C_Ys,"Credit/Output",var_MPK[1,:],"Variance of MPK for Entrepreneurs", false, 12,7, "ENT", 0.0)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_var_mpk_ent.png")
plts = create_plots(C_Ys,"Credit/Output", [var_MPK[2,:], var_MPK[3,:]],["Variance of MPK for Sole Proprietors","Variance of MPK for Employers"],false, 12,7,["SP","EMP"],0.0)
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_var_mpk_sp.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_var_mpk_emp.png")
LABELS = ["SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [var_MPK[2,:], var_MPK[3,:]],LABELS,"Variance of MPK",false, 12,7,["SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_credit_to_gdp_var_mpk.png")

# share of unbound ENT, SP,EMP
plts = create_plots(C_Ys,"Credit/Output", [share_unbound[1,:], share_unbound[2,:], share_unbound[3,:]],["Share of Unconstrained Entrepreneurs","Share of Unconstrained Sole Proprietors","Share of Unconstrained Employers"],true, 12,7,["ENT","SP","EMP"])
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_ent_unbound.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_se_unbound.png")
savefig(plts[3],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_emp_unbound.png")

# average and variance of working and managerial skills across occupations
CHOICE_NAMES = ["Workers","Entrepreneurs","Sole Proprietors","Employers"]
LABELS = ["Average Managerial Skill across $(occup)" for occup in CHOICE_NAMES]
plts = create_plots(C_Ys,"Credit/Output", [avg_m_skill[1,:], avg_m_skill[2,:], avg_m_skill[3,:], avg_m_skill[4,:]],LABELS,false, 12,7,["W","ENT","SP","EMP"])
for h in 1:length(plts)
    display(plts[h])
    occup = CHOICE_NAMES[h]
    if h==3
        occup = "SoleProprietors"
    end
    savefig(plts[h],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_avg_m_skill_$(occup).png")
end
LABELS = ["W","ENT","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [avg_m_skill[1,:], avg_m_skill[2,:], avg_m_skill[3,:], avg_m_skill[4,:]],LABELS,"Average Managerial Skill",false, 12,7,["W","ENT","SP","EMP"],-Inf, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_avg_m_skill_with_ENT.png")
LABELS = ["W","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [avg_m_skill[1,:], avg_m_skill[3,:], avg_m_skill[4,:]],LABELS,"Average Managerial Skill",false, 12,7,["W","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_avg_m_skill.png")


LABELS = ["Average Working Skill across $(occup)" for occup in CHOICE_NAMES]
plts = create_plots(C_Ys,"Credit/Output", [avg_w_skill[1,:], avg_w_skill[2,:], avg_w_skill[3,:], avg_w_skill[4,:]],LABELS,false, 12,7,["W","ENT","SP","EMP"])
for h in 1:length(plts)
    display(plts[h])
    occup = CHOICE_NAMES[h]
    if h==3
        occup = "SoleProprietors"
    end
    savefig(plts[h],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_avg_w_skill_$(occup).png")
end
LABELS = ["W","ENT","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [avg_w_skill[1,:], avg_w_skill[2,:], avg_w_skill[3,:], avg_w_skill[4,:]],LABELS,"Average Working Skill",false, 12,7,["W","ENT","SP","EMP"],-Inf, false)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_avg_w_skill_with_ENT.png")
LABELS = ["W","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [avg_w_skill[1,:], avg_w_skill[3,:], avg_w_skill[4,:]],LABELS,"Average Working Skill",false, 12,7,["W","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_avg_w_skill.png")


LABELS = ["Variance of Managerial Skill across $(occup)" for occup in CHOICE_NAMES]
plts = create_plots(C_Ys,"Credit/Output", [var_m_skill[1,:], var_m_skill[2,:], var_m_skill[3,:], var_m_skill[4,:]],LABELS,false, 12,7,["W","ENT","SP","EMP"], 0.0)
for h in 1:length(plts)
    display(plts[h])
    occup = CHOICE_NAMES[h]
    if h==3
        occup = "SoleProprietors"
    end
    savefig(plts[h],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_var_m_skill_$(occup).png")
end
LABELS = ["W","ENT","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [var_m_skill[1,:], var_m_skill[2,:], var_m_skill[3,:], var_m_skill[4,:]],LABELS,"Variance of Managerial Skill",false, 12,7,["W","ENT","SP","EMP"],0.0, :right)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_var_m_skill_with_ENT.png")
LABELS = ["W","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [var_m_skill[1,:], var_m_skill[3,:], var_m_skill[4,:]],LABELS,"Variance of Managerial Skill",false, 12,7,["W","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_var_m_skill.png")


LABELS = ["Variance of Working Skill across $(occup)" for occup in CHOICE_NAMES]
plts = create_plots(C_Ys,"Credit/Output", [var_w_skill[1,:], var_w_skill[2,:], var_w_skill[3,:], var_w_skill[4,:]],LABELS,false, 12,7,["W","ENT","SP","EMP"])
for h in 1:length(plts)
    display(plts[h])
    occup = CHOICE_NAMES[h]
    if h==3
        occup = "SoleProprietors"
    end
    savefig(plts[h],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_var_w_skill_$(occup).png")
end
LABELS = ["W","ENT","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [var_w_skill[1,:], var_w_skill[2,:], var_w_skill[3,:], var_w_skill[4,:]],LABELS,"Variance of Working Skill",false, 12,7,["W","ENT","SP","EMP"], 0.0, :topright)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_var_w_skill_with_ENT.png")
LABELS = ["W","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [var_w_skill[1,:], var_w_skill[3,:], var_w_skill[4,:]],LABELS,"Variance of Working Skill",false, 12,7,["W","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_var_w_skill.png")

LABELS = ["Sum of Managerial Skill across $(occup)" for occup in CHOICE_NAMES]
plts = create_plots(C_Ys,"Credit/Output", [sum_m_skill[1,:], sum_m_skill[2,:], sum_m_skill[3,:], sum_m_skill[4,:]],LABELS,false, 12,7,["W","ENT","SP","EMP"])
for h in 1:length(plts)
    display(plts[h])
    occup = CHOICE_NAMES[h]
    if h==3
        occup = "SoleProprietors"
    end
    savefig(plts[h],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_sum_m_skill_$(occup).png")
end
LABELS = ["W","ENT","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [sum_m_skill[1,:], sum_m_skill[2,:], sum_m_skill[3,:], sum_m_skill[4,:]],LABELS,"Sum of Managerial Skill",false, 12,7,["W","ENT","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_sum_m_skill_with_ENT.png")
LABELS = ["W","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [sum_m_skill[1,:], sum_m_skill[3,:], sum_m_skill[4,:]],LABELS,"Sum of Managerial Skill",false, 12,7,["W","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_sum_m_skill.png")


LABELS = ["Sum of Working Skill across $(occup)" for occup in CHOICE_NAMES]
plts = create_plots(C_Ys,"Credit/Output", [sum_w_skill[1,:], sum_w_skill[2,:], sum_w_skill[3,:], sum_w_skill[4,:]],LABELS,false, 12,7,["W","ENT","SP","EMP"])
for h in 1:length(plts)
    display(plts[h])
    occup = CHOICE_NAMES[h]
    if h==3
        occup = "SoleProprietors"
    end
    savefig(plts[h],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_sum_w_skill_$(occup).png")
end
LABELS = ["W","ENT","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [sum_w_skill[1,:], sum_w_skill[2,:], sum_w_skill[3,:], sum_w_skill[4,:]],LABELS,"Sum of Working Skill",false, 12,7,["W","ENT","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_sum_w_skill_with_ENT.png")
LABELS = ["W","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [sum_w_skill[1,:], sum_w_skill[3,:], sum_w_skill[4,:]],LABELS,"Sum of Working Skill",false, 12,7,["W","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_sum_w_skill.png")


plts = create_plots(C_Ys,"Credit/Output", [share_W_earnings_in_output, share_SP_earnings_in_output, share_EMP_earnings_in_output, share_W_capital_income_in_output, share_SP_capital_income_in_output, share_EMP_capital_income_in_output],["Share of output as Workers' Earnings","Share of output as Sole Proprietors' Earnings","Share of output as Employers' Earnings","Share of output as Workers' Capital Income","Share of output as Sole Proprietors' Capital Income","Share of output as Employers' Capital Income"],true, 12,7,["W","SP","EMP","W","SP","EMP"])
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_of_output_W_earnings.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_of_output_SP_earnings.png")
savefig(plts[3],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_of_output_EMP_earnings.png")
savefig(plts[4],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_of_output_W_capital_income.png")
savefig(plts[5],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_of_output_SP_capital_income.png")
savefig(plts[6],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_of_output_EMP_capital_income.png")


# ratio of capital restricted to unrestricted
plt = create_plot(C_Ys,"Credit/Output", ratio_of_capitals,"Ratio of Restricted Capital to Unrestricted", true)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_ratio_of_capital.png")
# ratio of output restricted to unrestricted
plt = create_plot(C_Ys,"Credit/Output", ratio_of_outputs,"Ratio of Restricted Output to Unrestricted", true)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_ratio_of_output.png")


#TFP_ideal, TFP_data
plt = create_plot(C_Ys,"Credit/Output", old_TFPis,"TFP", false, 12,7)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_tfp_ideal.png")
plt = create_plot(C_Ys,"Credit/Output", old_TFPds,"TFP_data", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_tfp_data.png")

#Labour_productivity, Capital_productivity, Investment_productivity
plt = create_plot(C_Ys,"Credit/Output", Labour_productivity_s,"Output-to-labour", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_labour_productivity.png")
plt = create_plot(C_Ys,"Credit/Output", Capital_productivity_s,"Output-to-capital", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_capital_productivity.png")
plt = create_plot(C_Ys,"Credit/Output", Investment_productivity_s,"Expenses-to-revenue", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_investment_productivity.png")

# Print statistics for calibrated economy

# income, earnings, wealth, consumption
# All, W, SP, EMP
# lambdas
#means = zeros(4,4,num_lambdas)
#ginis = zeros(4,4,num_lambdas)
#vars = zeros(4,4,num_lambdas)
for i=1:4

    if i==1
        stat_name="income"
    elseif i==2
        stat_name="earnings"
    elseif i==3
        stat_name="wealth"
    else  #i==4
        stat_name="consumption"
    end

    for j=[4]#1:4
        if j==1
            stat_type_name="Mean"
            stat_type = means
        elseif j==2
            stat_type_name="Variance"
            stat_type = vars
        elseif j==3
            stat_type_name="Varlog"
            stat_type = varlogs
        else  #i==4
            stat_type_name="Gini"
            stat_type = ginis
        end
        println("$(stat_type_name) of $(stat_name) for H,W,SP,EMP: $(round.(stat_type[i,:,calibrated_lambda];digits=2))")
    end
end

ss_star = SSS[10]
asset_grid = ss_star[1][3]
density_distr = ss_star[1][5]
wealth = ones(size(density_distr)).*asset_grid
policy = ss_star[1][4]
income = ss_star[1][23] .- wealth
earnings = ss_star[1][24]
consumption = income .+ wealth .- policy
occ_choice = ss_star[1][22]

h_choice = Float64.(ss_star[1][22].!=0.0)
w_choice = Float64.(ss_star[1][22].==1.0)
sp_choice = Float64.(ss_star[1][22].==2.0)
emp_choice = Float64.(ss_star[1][22].==3.0)

num_of_occs = 4

mean_consumption_to_income_wealth = Array{Any}(undef,num_of_occs)
var_consumption_to_income_wealth = Array{Any}(undef,num_of_occs)
std_consumption_to_income_wealth = Array{Any}(undef,num_of_occs)

mean_income_to_income_wealth = Array{Any}(undef,num_of_occs)
var_income_to_income_wealth = Array{Any}(undef,num_of_occs)
std_income_to_income_wealth = Array{Any}(undef,num_of_occs)

mean_consumption_to_income = Array{Any}(undef,num_of_occs)
var_consumption_to_income = Array{Any}(undef,num_of_occs)
std_consumption_to_income = Array{Any}(undef,num_of_occs)

mean_savings_to_income = Array{Any}(undef,num_of_occs)
var_savings_to_income = Array{Any}(undef,num_of_occs)
std_savings_to_income = Array{Any}(undef,num_of_occs)

mean_consumption_to_earnings = Array{Any}(undef,num_of_occs)
var_consumption_to_earnings = Array{Any}(undef,num_of_occs)
std_consumption_to_earnings = Array{Any}(undef,num_of_occs)

mean_savings_to_earnings = Array{Any}(undef,num_of_occs)
var_savings_to_earnings = Array{Any}(undef,num_of_occs)
std_savings_to_earnings = Array{Any}(undef,num_of_occs)

for i = 1:num_of_occs
    choice = h_choice
    if i==2
        choice = w_choice
    elseif i==3
        choice = sp_choice
    elseif i==4
        choice = emp_choice
    end

    mean_consumption_to_income_wealth[i] = sum((consumption./ss_star[1][23]).*choice.*density_distr)
    var_consumption_to_income_wealth[i] = sum((consumption./ss_star[1][23] .- mean_consumption_to_income_wealth[i]).^2 .*choice.*density_distr./sum(choice.*density_distr))
    mean_consumption_to_income_wealth[i] /= sum(choice.*density_distr)
    std_consumption_to_income_wealth[i] = sqrt(var_consumption_to_income_wealth[i])
    #display(mean_consumption_to_income_wealth[i])
    #display(var_consumption_to_income_wealth[i])
    #display(std_consumption_to_income_wealth[i])


    mean_income_to_income_wealth[i] = sum((income./ss_star[1][23]).*choice.*density_distr)
    var_income_to_income_wealth[i] = sum((income./ss_star[1][23] .- mean_income_to_income_wealth[i]).^2 .*choice.*density_distr./sum(choice.*density_distr))
    mean_income_to_income_wealth[i] /= sum(choice.*density_distr)
    std_income_to_income_wealth[i] = sqrt(var_income_to_income_wealth[i])
    #display(mean_income_to_income_wealth[i])
    #display(var_income_to_income_wealth[i])
    #display(std_income_to_income_wealth[i])


    mean_consumption_to_income[i] = sum((consumption./income).*choice.*density_distr)
    var_consumption_to_income[i] = sum((consumption./income .- mean_consumption_to_income[i]).^2 .*choice.*density_distr./sum(choice.*density_distr))
    mean_consumption_to_income[i] /= sum(choice.*density_distr)
    std_consumption_to_income[i] = sqrt(var_consumption_to_income[i])
    #display(mean_consumption_to_income[i])
    #display(var_consumption_to_income[i])
    #display(std_consumption_to_income[i])

    mean_savings_to_income[i] = sum(((policy.-wealth)./income).*choice.*density_distr)
    var_savings_to_income[i] = sum(((policy.-wealth)./income .- mean_savings_to_income[i]).^2 .*choice.*density_distr./sum(choice.*density_distr))
    mean_savings_to_income[i] /= sum(choice.*density_distr)
    std_savings_to_income[i] = sqrt(var_savings_to_income[i])
    #display(mean_savings_to_income[i])
    #display(var_savings_to_income[i])
    #display(std_savings_to_income[i])

    mean_consumption_to_earnings[i] = sum((consumption./earnings).*choice.*density_distr)
    var_consumption_to_earnings[i] = sum((consumption./earnings .- mean_consumption_to_earnings[i]).^2 .*choice.*density_distr./sum(choice.*density_distr))
    mean_consumption_to_earnings[i] /= sum(choice.*density_distr)
    std_consumption_to_earnings[i] = sqrt(var_consumption_to_earnings[i])
    #display(mean_consumption_to_income[i])
    #display(var_consumption_to_income[i])
    #display(std_consumption_to_income[i])

    mean_savings_to_earnings[i] = sum(((policy.-wealth)./earnings).*choice.*density_distr)
    var_savings_to_earnings[i] = sum(((policy.-wealth)./earnings .- mean_savings_to_earnings[i]).^2 .*choice.*density_distr./sum(choice.*density_distr))
    mean_savings_to_earnings[i] /= sum(choice.*density_distr)
    std_savings_to_earnings[i] = sqrt(var_savings_to_earnings[i])
    #display(mean_savings_to_income[i])
    #display(var_savings_to_income[i])
    #display(std_savings_to_income[i])

end

println("Mean of the ratio of C/(Income+Wealth) for H,W,SP,EMP: $(round.(mean_consumption_to_income_wealth;digits=2))")
println("Std of the ratio of C/(Income+Wealth) for H,W,SP,EMP: $(round.(std_consumption_to_income_wealth;digits=2))")

println("Mean of the ratio of Income/(Income+Wealth) for H,W,SP,EMP: $(round.(mean_income_to_income_wealth;digits=2))")
println("Std of the ratio of Income/(Income+Wealth) for H,W,SP,EMP: $(round.(std_income_to_income_wealth;digits=2))")

println("Mean of the ratio of C/Income for H,W,SP,EMP: $(round.(mean_consumption_to_income;digits=2))")
println("Std of the ratio of C/Income for H,W,SP,EMP: $(round.(std_consumption_to_income;digits=2))")

println("Mean of the ratio of Savings/Income for H,W,SP,EMP: $(round.(mean_savings_to_income;digits=2))")
println("Std of the ratio of Savings/Income for H,W,SP,EMP: $(round.(std_savings_to_income;digits=2))")

println("Mean of the ratio of C/Earnings for H,W,SP,EMP: $(round.(mean_consumption_to_earnings;digits=2))")
println("Std of the ratio of C/Earnings for H,W,SP,EMP: $(round.(std_consumption_to_earnings;digits=2))")

println("Mean of the ratio of Savings/Earnings for H,W,SP,EMP: $(round.(mean_savings_to_earnings;digits=2))")
println("Std of the ratio of Savings/Earnings for H,W,SP,EMP: $(round.(std_savings_to_earnings;digits=2))")

@save "$(LOCAL_DIR_GENERAL)SSS.jld2" SSS C_Ys Outputs Incomes Consumptions Rs Ws logcs loges giniWs giniEnts share_unbound means ginis avglogs varlogs avgs vars quantile_means TFPis TFPds mean_MPL var_MPL mean_MPK var_MPK Capital share_W_earnings_in_output share_SP_earnings_in_output share_EMP_earnings_in_output share_W_capital_income_in_output share_SP_capital_income_in_output share_EMP_capital_income_in_output
