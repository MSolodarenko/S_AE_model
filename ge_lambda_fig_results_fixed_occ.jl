include("Functions/print_sameline.jl")
############################
######### Libraries ########
print_sameline("Initialising libraries")
using Plots
#using PyPlot

using SchumakerSpline
using Statistics
using JLD2
using ProgressMeter

# global parameters of the approximation objects
#                               1               2               3               4                       5                   6
#                       number_a_nodes, number_u_m_nodes, number_u_w_nodes, number_zeta_nodes, number_alpha_m_nodes, number_alpha_w_nodes
GLOBAL_APPROX_PARAMS = [69,3,3,3,6,3]#[35,3,3,3,6,3]

CODENAME = #="SS_2642"#"SS_2065"=#"SS_2092"
CODENAME = "$(CODENAME)_$(GLOBAL_APPROX_PARAMS[1])"
# parameters of the model's economy (Italy)
                    #SS_2642 #SS_2065 #SS_2092 #prev.calibration
LAMBDA =            #=1.633951#1.405096=#1.665907 #1.513028
country = "Italy"
LOCAL_DIR = "$(@__DIR__)/Results/Fixed_occ_shares/Lambda_grid/$(country)_$(CODENAME)/"
if Sys.iswindows()
    LOCAL_DIR = "\\Results\\Fixed_occ_shares\\Lambda_grid\\$(country)_$(CODENAME)\\"
end
global_approx_params = copy(GLOBAL_APPROX_PARAMS)

print_sameline("Loading data from Fixed_occ_shares/Lambda_grid")
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

num_lambdas = length(lambda)#20# #instead of length(lambda)
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

quantile_mean_wealth = zeros(5,num_lambdas)
quantile_mean_wealth_w = zeros(5,num_lambdas)
quantile_mean_wealth_sp = zeros(5,num_lambdas)
quantile_mean_wealth_emp = zeros(5,num_lambdas)

interquartile_wealth = zeros(num_lambdas)
interquintile_wealth = zeros(num_lambdas)
interdecile_wealth = zeros(num_lambdas)
interquartile_wealth_w = zeros(num_lambdas)
interquintile_wealth_w = zeros(num_lambdas)
interdecile_wealth_w = zeros(num_lambdas)
interquartile_wealth_se = zeros(num_lambdas)
interquintile_wealth_se = zeros(num_lambdas)
interdecile_wealth_se = zeros(num_lambdas)
interquartile_wealth_emp = zeros(num_lambdas)
interquintile_wealth_emp = zeros(num_lambdas)
interdecile_wealth_emp = zeros(num_lambdas)

palma_wealth = zeros(num_lambdas)
palma_wealth_w = zeros(num_lambdas)
palma_wealth_se = zeros(num_lambdas)
palma_wealth_emp = zeros(num_lambdas)

ratio2020 = zeros(num_lambdas)
ratio2020_w = zeros(num_lambdas)
ratio2020_se = zeros(num_lambdas)
ratio2020_emp = zeros(num_lambdas)

gini_wealth = zeros(num_lambdas)
gini_wealth_w = zeros(num_lambdas)
gini_wealth_se = zeros(num_lambdas)
gini_wealth_emp = zeros(num_lambdas)

Capital = zeros(num_lambdas)

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

output_SPs  = zeros(num_lambdas)
output_EMPs = zeros(num_lambdas)

share_of_output_SPs  = zeros(num_lambdas)
share_of_output_EMPs = zeros(num_lambdas)

ratio_of_capitals = zeros(num_lambdas)
ratio_of_outputs  = zeros(num_lambdas)

old_TFPis = zeros(num_lambdas)
old_TFPds = zeros(num_lambdas)
Labour_productivity_s = zeros(num_lambdas)
Capital_productivity_s = zeros(num_lambdas)
Investment_productivity_s = zeros(num_lambdas)

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
    policy = ss_star[1][4]
    income = ss_star[1][23] .- wealth
    earnings = ss_star[1][24]
    capital_d = ss_star[1][26]
    labour_d = ss_star[1][29]
    output = ss_star[1][32]

    wealth = Array{Any}(undef,3)
    income = Array{Any}(undef,3)
    consumption = Array{Any}(undef,3)
    for occ in 1:3
        wealth[occ] = ones(size(density_distr[occ])).*asset_grid
        income[occ] = ss_star[1][23][occ] .- wealth[occ]
        consumption[occ] = income[occ] .+ wealth[occ] .- policy[occ]

        next!(p)
    end

    Capital[i] = sum( [sum(asset_grid.*density_distr[occ]) for occ in 1:3] )

    Outputs[i] = sum( [sum(output[occ] .* density_distr[occ]) for occ in 1:3] )
    Incomes[i] = sum( [sum(income[occ] .* density_distr[occ]) for occ in 1:3] )
    Consumptions[i] = sum( [sum(consumption[occ] .* density_distr[occ]) for occ in 1:3] )
    Income_to_outputs[i] = Incomes[i]/Outputs[i]
    next!(p)

    number_u_nodes = SSS[i][4][7]
    number_asset_grid = ss_star[1][2]
    number_zeta_nodes = global_approx_params[4]
    number_alpha_m_nodes = global_approx_params[5]
    number_alpha_w_nodes = global_approx_params[6]

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

        for h = 1:4 # [2] = All, W, SP, EMP
            temp_occ = 4
            if h == 1 #All
                # calculate mean
                means[s,h,i] = sum( [sum(stat_distr[occ] .* density_distr[occ]) for occ in 1:3])/sum( [sum(density_distr[occ]) for occ in 1:3])

                # calculate gini coefficent
                density_distr_choice = [density_distr[2]; density_distr[3]]./sum([density_distr[2]; density_distr[3]])
                stat_choice = [stat_distr[2]; stat_distr[3]]
                density_distr_choice_vec = vec(density_distr_choice)
                stat_choice_vec = vec(stat_choice)

                index_non_zero = findall(x-> x>1e-5,density_distr_choice_vec)
                density_distr_choice_vec_non_zero = density_distr_choice_vec[index_non_zero]
                stat_choice_vec_non_zero = stat_choice_vec[index_non_zero]

                ginis[s,h,i] = sum([ density_distr_choice_vec_non_zero[y_i]*density_distr_choice_vec_non_zero[y_j]*abs(stat_choice_vec_non_zero[y_i]-stat_choice_vec_non_zero[y_j]) for y_i=1:length(density_distr_choice_vec_non_zero), y_j=1:length(density_distr_choice_vec_non_zero) ]) / (2*sum(stat_choice_vec_non_zero.*density_distr_choice_vec_non_zero))

                # calculate variance of log-s
                avglogs[s,h,i] = sum([ sum(density_distr[occ].*log.(max.(1e-12,stat_distr[occ]))   ) for occ in 1:3])
                varlogs[s,h,i] = sum([ sum(density_distr[occ].*(log.(max.(1e-12,stat_distr[occ])).-avglogs[s,h,i]).^2) for occ in 1:3])/sum([ sum(density_distr[occ]) for occ in 1:3])
                avglogs[s,h,i] /= sum([ sum(density_distr[occ]) for occ in 1:3])
                if s==3
                    avglogs[s,h,i] = sum([ sum(density_distr[occ].*max.(1e-12,stat_distr[occ])   ) for occ in 1:3])
                    varlogs[s,h,i] = sum([ sum(density_distr[occ].*(max.(1e-12,stat_distr[occ]).-avglogs[s,h,i]).^2) for occ in 1:3])/sum([ sum(density_distr[occ]) for occ in 1:3])
                    avglogs[s,h,i] /= sum([ sum(density_distr[occ]) for occ in 1:3])
                end

                avgs[s,h,i] = sum([ sum(density_distr[occ].*max.(1e-12,stat_distr[occ])   ) for occ in 1:3])
                vars[s,h,i] = sum([ sum(density_distr[occ].*(max.(1e-12,stat_distr[occ]).-avgs[s,h,i]).^2) for occ in 1:3])/sum([ sum(density_distr[occ]) for occ in 1:3])
                avgs[s,h,i] /= sum([ sum(density_distr[occ]) for occ in 1:3])
            elseif h == 2 #W
                temp_occ = 1
            elseif h == 3 #SP
                temp_occ = 2
            elseif h == 4 #EMP
                temp_occ = 3
            end
            occ = temp_occ

            if h!=1
                # calculate mean
                means[s,h,i] = sum(stat_distr[occ] .* density_distr[occ])/sum(density_distr[occ])

                # calculate gini coefficent
                density_distr_choice = (density_distr[occ])./sum(density_distr[occ])
                stat_choice = stat_distr[occ]
                density_distr_choice_vec = vec(density_distr_choice)
                stat_choice_vec = vec(stat_choice)

                index_non_zero = findall(x-> x>1e-5,density_distr_choice_vec)
                density_distr_choice_vec_non_zero = density_distr_choice_vec[index_non_zero]
                stat_choice_vec_non_zero = stat_choice_vec[index_non_zero]

                ginis[s,h,i] = sum([ density_distr_choice_vec_non_zero[y_i]*density_distr_choice_vec_non_zero[y_j]*abs(stat_choice_vec_non_zero[y_i]-stat_choice_vec_non_zero[y_j]) for y_i=1:length(density_distr_choice_vec_non_zero), y_j=1:length(density_distr_choice_vec_non_zero) ]) / (2*sum(stat_choice_vec_non_zero.*density_distr_choice_vec_non_zero))

                # calculate variance of log-s
                avglogs[s,h,i] = sum(density_distr[occ].*log.(max.(1e-12,stat_distr[occ]))   )
                varlogs[s,h,i] = sum(density_distr[occ].*(log.(max.(1e-12,stat_distr[occ])).-avglogs[s,h,i]).^2)/sum(density_distr[occ])
                avglogs[s,h,i] /= sum(density_distr[occ])
                if s==3
                    avglogs[s,h,i] = sum(density_distr[occ].*max.(1e-12,stat_distr[occ])   )
                    varlogs[s,h,i] = sum(density_distr[occ].*(max.(1e-12,stat_distr[occ]).- avglogs[s,h,i]).^2)/sum(density_distr[occ])
                    avglogs[s,h,i] /= sum(density_distr[occ])
                end

                avgs[s,h,i] = sum(density_distr[occ].*max.(1e-12,stat_distr[occ])   )
                vars[s,h,i] = sum(density_distr[occ].*(max.(1e-12,stat_distr[occ]).- avgs[s,h,i]).^2)/sum(density_distr[occ])
                avgs[s,h,i] /= sum(density_distr[occ])
            end
            next!(p)
        end
    end

    output_SPs[i]  = sum(ss_star[1][32][2].* density_distr[2])
    output_EMPs[i] = sum(ss_star[1][32][3].* density_distr[3])

    share_of_output_SPs[i]  = output_SPs[i]/Outputs[i]
    share_of_output_EMPs[i] = output_EMPs[i]/Outputs[i]

    ratio_of_capitals[i] = ss_star[1][6]


    asset_distr = sum(sum(density_distr), dims=2:5)[:,1,1,1,1]
    asset_distr ./= sum(asset_distr)
    asset_distr_w = sum(density_distr[1], dims=2:5)[:,1,1,1,1]
    asset_distr_w ./= sum(asset_distr_w)
    asset_distr_se = sum(density_distr[2], dims=2:5)[:,1,1,1,1]
    asset_distr_se ./= sum(asset_distr_se)
    asset_distr_emp = sum(density_distr[3], dims=2:5)[:,1,1,1,1]
    asset_distr_emp ./= sum(asset_distr_emp)

    cdf = cumsum(asset_distr)
    cdf_w = cumsum(asset_distr_w)
    cdf_se = cumsum(asset_distr_se)
    cdf_emp = cumsum(asset_distr_emp)

    max_a_i = findlast(x->x<1.0-1e-5, cdf)

    cdf_1 = Schumaker(cdf[1:max_a_i]./cdf[max_a_i], asset_grid[1:max_a_i]; extrapolation = (Constant,Constant))
    cdf_w_1 = Schumaker(cdf_w[1:max_a_i]./cdf_w[max_a_i], asset_grid[1:max_a_i]; extrapolation = (Constant,Constant))
    cdf_se_1 = Schumaker(cdf_se[1:max_a_i]./cdf_se[max_a_i], asset_grid[1:max_a_i]; extrapolation = (Constant,Constant))
    cdf_emp_1 = Schumaker(cdf_emp[1:max_a_i]./cdf_emp[max_a_i], asset_grid[1:max_a_i]; extrapolation = (Constant,Constant))

    quantiles = collect(range(0.0;stop=1.0,length=6))
    quantile_wealth_bounds = evaluate.(cdf_1,quantiles)
    quantile_wealth_bounds_w = evaluate.(cdf_w_1,quantiles)
    quantile_wealth_bounds_sp = evaluate.(cdf_se_1,quantiles)
    quantile_wealth_bounds_emp = evaluate.(cdf_emp_1,quantiles)

    for occ = 1:4
        if occ==1
            qwb = quantile_wealth_bounds
            ad = asset_distr
        elseif occ==2
            qwb = quantile_wealth_bounds_w
            ad = asset_distr_w
        elseif occ==3
            qwb = quantile_wealth_bounds_sp
            ad = asset_distr_se
        else#if occ==4
            qwb = quantile_wealth_bounds_emp
            ad = asset_distr_emp
        end
        for q in 1:5
            m_w_down_down = sum(asset_grid[Int32(floor(qwb[q])):Int32(floor(qwb[q+1]))] .* ad[Int32(floor(qwb[q])):Int32(floor(qwb[q+1]))]) / sum(ad[Int32(floor(qwb[q])):Int32(floor(qwb[q+1]))])
            m_w_down_up = sum(asset_grid[Int32(floor(qwb[q])):Int32(ceil(qwb[q+1]))] .* ad[Int32(floor(qwb[q])):Int32(ceil(qwb[q+1]))]) / sum(ad[Int32(floor(qwb[q])):Int32(ceil(qwb[q+1]))])
            m_w_up_down = sum(asset_grid[Int32(ceil(qwb[q])):Int32(floor(qwb[q+1]))] .* ad[Int32(ceil(qwb[q])):Int32(floor(qwb[q+1]))]) / sum(ad[Int32(ceil(qwb[q])):Int32(floor(qwb[q+1]))])
            m_w_up_up = sum(asset_grid[Int32(ceil(qwb[q])):Int32(ceil(qwb[q+1]))] .* ad[Int32(ceil(qwb[q])):Int32(ceil(qwb[q+1]))]) / sum(ad[Int32(ceil(qwb[q])):Int32(ceil(qwb[q+1]))])
            qmw = (m_w_down_down+m_w_down_up+m_w_up_down+m_w_up_up)/4
            if occ==1
                quantile_mean_wealth[q,i] = qmw
            elseif occ==2
                quantile_mean_wealth_w[q,i] = qmw
            elseif occ==3
                quantile_mean_wealth_sp[q,i] = qmw
            else#if occ==4
                quantile_mean_wealth_emp[q,i] = qmw
            end
        end
    end

    #Q3-Q1
    interquartile_wealth[i] = evaluate(cdf_1,0.75) - evaluate(cdf_1,0.25)
    #Q8-Q2
    interquintile_wealth[i] = evaluate(cdf_1,0.8) - evaluate(cdf_1,0.2)
    #Q9-Q1
    interdecile_wealth[i] = evaluate(cdf_1,0.9) - evaluate(cdf_1,0.1)
    interquartile_wealth_w[i] = evaluate(cdf_w_1,0.75) - evaluate(cdf_w_1,0.25)
    interquintile_wealth_w[i] = evaluate(cdf_w_1,0.8) - evaluate(cdf_w_1,0.2)
    interdecile_wealth_w[i] = evaluate(cdf_w_1,0.9) - evaluate(cdf_w_1,0.1)
    interquartile_wealth_se[i] = evaluate(cdf_se_1,0.75) - evaluate(cdf_se_1,0.25)
    interquintile_wealth_se[i] = evaluate(cdf_se_1,0.8) - evaluate(cdf_se_1,0.2)
    interdecile_wealth_se[i] = evaluate(cdf_se_1,0.9) - evaluate(cdf_se_1,0.1)
    interquartile_wealth_emp[i] = evaluate(cdf_emp_1,0.75) - evaluate(cdf_emp_1,0.25)
    interquintile_wealth_emp[i] = evaluate(cdf_emp_1,0.8) - evaluate(cdf_emp_1,0.2)
    interdecile_wealth_emp[i] = evaluate(cdf_emp_1,0.9) - evaluate(cdf_emp_1,0.1)

    findnearest(A::AbstractArray,t) = findmin(abs.(A.-t))[2]

    #calc_share_of_wealth_between_shares_of_population
    function csowbsop(top_percent, bottom_percent, CDF_1, A_DISTR)
        a_top = evaluate(CDF_1, top_percent)
        a_bottom = evaluate(CDF_1, bottom_percent)
        #display([a_top, a_bottom])

        a_i_top = findnearest(asset_grid[1:max_a_i], a_top)
        a_i_bottom = findnearest(asset_grid[1:max_a_i], a_bottom)
        #display([a_i_top, a_i_bottom])

        #display([sum(asset_grid[a_i_bottom:a_i_top].*A_DISTR[a_i_bottom:a_i_top]),sum(asset_grid[1:max_a_i].*A_DISTR[1:max_a_i]), sum(asset_grid[a_i_bottom:a_i_top].*A_DISTR[a_i_bottom:a_i_top]) / sum(asset_grid[1:max_a_i].*A_DISTR[1:max_a_i])])
        return sum(asset_grid[a_i_bottom:a_i_top].*A_DISTR[a_i_bottom:a_i_top]) / sum(asset_grid[1:max_a_i].*A_DISTR[1:max_a_i])
    end

    palma_wealth[i] = csowbsop(1.0,0.9,cdf_1,asset_distr)/csowbsop(0.4,0.0,cdf_1,asset_distr)
    palma_wealth_w[i] = csowbsop(1.0,0.9,cdf_w_1,asset_distr_w)/csowbsop(0.4,0.0,cdf_w_1,asset_distr_w)
    palma_wealth_se[i] = csowbsop(1.0,0.9,cdf_se_1,asset_distr_se)/csowbsop(0.4,0.0,cdf_se_1,asset_distr_se)
    palma_wealth_emp[i] = csowbsop(1.0,0.9,cdf_emp_1,asset_distr_emp)/csowbsop(0.4,0.0,cdf_emp_1,asset_distr_emp)

    ratio2020[i] = csowbsop(1.0,0.8,cdf_1,asset_distr)/csowbsop(0.2,0.0,cdf_1,asset_distr)
    ratio2020_w[i] = csowbsop(1.0,0.8,cdf_w_1,asset_distr)/csowbsop(0.2,0.0,cdf_w_1,asset_distr)
    ratio2020_se[i] = csowbsop(1.0,0.8,cdf_se_1,asset_distr)/csowbsop(0.2,0.0,cdf_se_1,asset_distr)
    ratio2020_emp[i] = csowbsop(1.0,0.8,cdf_emp_1,asset_distr)/csowbsop(0.2,0.0,cdf_emp_1,asset_distr)

    gini_wealth[i] = sum([ asset_distr[y_i]*asset_distr[y_j]*abs(asset_grid[y_i]-asset_grid[y_j]) for y_i=1:max_a_i, y_j=1:max_a_i ]) / (2*sum(asset_grid[1:max_a_i].*asset_distr[1:max_a_i]))
    gini_wealth_w[i] = sum([ asset_distr_w[y_i]*asset_distr_w[y_j]*abs(asset_grid[y_i]-asset_grid[y_j]) for y_i=1:max_a_i, y_j=1:max_a_i ]) / (2*sum(asset_grid[1:max_a_i].*asset_distr_w[1:max_a_i]))
    gini_wealth_se[i] = sum([ asset_distr_se[y_i]*asset_distr_se[y_j]*abs(asset_grid[y_i]-asset_grid[y_j]) for y_i=1:max_a_i, y_j=1:max_a_i ]) / (2*sum(asset_grid[1:max_a_i].*asset_distr_se[1:max_a_i]))
    gini_wealth_emp[i] = sum([ asset_distr_emp[y_i]*asset_distr_emp[y_j]*abs(asset_grid[y_i]-asset_grid[y_j]) for y_i=1:max_a_i, y_j=1:max_a_i ]) / (2*sum(asset_grid[1:max_a_i].*asset_distr_emp[1:max_a_i]))

    unlimited_capital_choice = Array{Any}(undef,3)
    for occ in 1:3
        unlimited_capital_choice[occ] = capital_d[occ].<(lambdas[i].*(asset_grid.*ones(size(capital_d[occ]))))
    end

    for c = 1:4
        temp_occ = 4
        if c==1
            temp_occ = 1
        elseif c==2
            denumerator = sum([ sum(density_distr[occ]) for occ in 1:3])

            Y = sum([ sum(density_distr[occ].*output[occ]) for occ in 1:3])/denumerator
            K = sum([ sum(density_distr[occ].*capital_d[occ]) for occ in 1:3])/denumerator
            L = sum([ sum(density_distr[occ].*labour_d[occ]) for occ in 1:3])/denumerator
            TFPis[c-1,i] = Y/(K^eta*L^theta)
            TFPds[c-1,i] = Y/K^(eta/(theta+eta))

            mean_MPK[c-1,i] = sum([ sum(density_distr[occ].*(eta.*output[occ].*replace(capital_d[occ].^(-1.0),Inf=>0.0))) for occ in 1:3])/denumerator
            var_MPK[c-1,i] = sum([ sum(density_distr[occ].*(eta.*output[occ].*replace(capital_d[occ].^(-1.0),Inf=>0.0) .- mean_MPK[c-1,i]*denumerator).^2) for occ in 1:3])/denumerator

            mean_MPL[c-1,i] = sum([ sum(density_distr[occ].*(theta.*output[occ].*replace(labour_d[occ].^(-1.0),Inf=>0.0))) for occ in 1:3])/denumerator
            var_MPL[c-1,i] = sum([ sum(density_distr[occ].*(theta.*output[occ].*replace(labour_d[occ].^(-1.0),Inf=>0.0) .- mean_MPL[c-1,i]*denumerator).^2) for occ in 1:3])/denumerator

            #if c==2
            #    choice = ent_choice
            #end
            share_unbound[c-1,i] = sum([ sum(density_distr[occ].*unlimited_capital_choice[occ]) for occ in 1:3])/denumerator

            #[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
            #m = u_i,alpha_m_i,zeta_i
            m_skill_distr = permutedims(sum( sum(density_distr[2:3]), dims=[1,5])[1,:,:,:,1], [1,3,2])
            #w = u_i,alpha_m_i,alpha_w_i
            w_skill_distr = sum( sum(density_distr[2:3]), dims=[1,3])[1,:,1,:,:]

            avg_m_skill[c,i] = sum(m_skill_distr.*z_m_nodes)/denumerator
            avg_w_skill[c,i] = sum(w_skill_distr.*z_w_nodes)/denumerator

            var_m_skill[c,i] = sum(m_skill_distr.*(z_m_nodes .- avg_m_skill[c,i]*denumerator).^2)/denumerator
            var_w_skill[c,i] = sum(w_skill_distr.*(z_w_nodes .- avg_w_skill[c,i]*denumerator).^2)/denumerator

            sum_m_skill[c,i] = sum(m_skill_distr.*z_m_nodes)
            sum_w_skill[c,i] = sum(w_skill_distr.*z_w_nodes)

        elseif c==3
            temp_occ = 2
        elseif c==4
            temp_occ = 3
        else
            throw("c is out of bounds")
        end
        occ = temp_occ

        if c!=2
            denumerator = sum(density_distr[occ])
            if c>1
                #if c==2
                #    choice = ones(size(choice))
                #end
                #if c==2
                #    denumerator = sum(choice.*density_distr)
                #end
                Y = sum(density_distr[occ].*output[occ])/denumerator
                K = sum(density_distr[occ].*capital_d[occ])/denumerator
                L = sum(density_distr[occ].*labour_d[occ])/denumerator
                TFPis[c-1,i] = Y/(K^eta*L^theta)
                TFPds[c-1,i] = Y/K^(eta/(theta+eta))

                c_1 = capital_d[occ].^(-1.0)
                replace!(c_1,Inf=>0.0)
                mean_MPK[c-1,i] = sum(density_distr[occ].*(eta.*output[occ].*c_1))/denumerator
                var_MPK[c-1,i] = sum(density_distr[occ].*(eta.*output[occ].*c_1 .- mean_MPK[c-1,i]*denumerator).^2)/denumerator

                l_1 = labour_d[occ].^(-1.0)
                replace!(l_1,Inf=>0.0)
                mean_MPL[c-1,i] = sum(density_distr[occ].*(theta.*output[occ].*l_1))/denumerator
                var_MPL[c-1,i] = sum(density_distr[occ].*(theta.*output[occ].*l_1 .- mean_MPL[c-1,i]*denumerator).^2)/denumerator

                #if c==2
                #    choice = ent_choice
                #end
                share_unbound[c-1,i] = sum(density_distr[occ].*unlimited_capital_choice[occ])/denumerator
            end

            #[a_i,u_i,zeta_i,alpha_m_i,alpha_w_i]
            #m = u_i,alpha_m_i,zeta_i
            m_skill_distr = permutedims(sum(density_distr[occ], dims=[1,5])[1,:,:,:,1], [1,3,2])
            #w = u_i,alpha_m_i,alpha_w_i
            w_skill_distr = sum(density_distr[occ], dims=[1,3])[1,:,1,:,:]

            avg_m_skill[c,i] = sum(m_skill_distr.*z_m_nodes)/denumerator
            avg_w_skill[c,i] = sum(w_skill_distr.*z_w_nodes)/denumerator

            var_m_skill[c,i] = sum(m_skill_distr.*(z_m_nodes .- avg_m_skill[c,i]*denumerator).^2)/denumerator
            var_w_skill[c,i] = sum(w_skill_distr.*(z_w_nodes .- avg_w_skill[c,i]*denumerator).^2)/denumerator

            #display([avg_w_skill[c,i],avg_m_skill[c,i], var_w_skill[c,i],var_m_skill[c,i]])
        end
        next!(p)
    end

    old_TFPis[i] = ss_star[1][43][end-4]
    old_TFPds[i] = ss_star[1][43][end-3]
    Labour_productivity_s[i] = ss_star[1][43][end-2]
    Capital_productivity_s[i] = ss_star[1][43][end-1]
    Investment_productivity_s[i] = ss_star[1][43][end]
    #throw(error)
end

asset_grid = SSS[1][1][3]
number_non_zero_asset_grid = 0
number_u_nodes = SSS[1][4][7]
number_zeta_nodes = global_approx_params[4]
number_alpha_m_nodes = global_approx_params[5]
number_alpha_w_nodes = global_approx_params[6]
@showprogress for (occ,(u_i,(zeta_i,(alpha_m_i,alpha_w_i)))) in collect(Iterators.product(1:3,Iterators.product(1:number_u_nodes,Iterators.product(1:number_zeta_nodes,Iterators.product(1:number_alpha_m_nodes,1:number_alpha_w_nodes)))))
    global number_non_zero_asset_grid
    temp = findlast(round.(SSS[1][1][5][occ][:,u_i,zeta_i,alpha_m_i,alpha_w_i]; digits=7) .!= 0.0)
    if isnothing(temp)
        temp = 0
    end
    if temp > number_non_zero_asset_grid
        number_non_zero_asset_grid = temp
    end
end
#policy = SSS[l][1][4]
num_lambdas_sp = findfirst(x -> x>0.75, C_Ys)-1

#throw(error)

function create_plot(X::Vector{Float64},XLABEL::String,Y::Vector{Float64},YLABEL::String, IS_Y_PERCENTAGE::Bool=true, OCCUPATION::String="H", LOW_LIMIT=-Inf)
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
    YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=Int(round(num_lambdas/2;digits=0))))
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
                    ylabel=YLABEL,
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
function create_plots(X::Vector{Float64},XLABEL::String,Ys,YLABELs, IS_Y_PERCENTAGE::Bool=true, OCCUPATION=["H", "W", "SP", "EMP"], LOW_LIMIT=-Inf)
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
        YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=Int(round(num_lambdas/2;digits=0))))
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
                        xlabel = XLABEL,
                        ylabel = YLABEL,
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
function create_combined_plot(X::Vector{Float64},XLABEL::String,Ys,YLABELs,YLABEL, IS_Y_PERCENTAGE::Bool=true, OCCUPATION=["H", "W", "SP", "EMP"], LOW_LIMIT=-Inf)


    YLIMS1 = max(minimum(minimum.(Ys)), LOW_LIMIT)
    YLIMS2 = maximum(maximum.(Ys))
    #YLIMS=(YLIMS1-0.01, YLIMS2+0.01)
    YLIMMARGIN = abs(YLIMS2-YLIMS1)*0.015
    YLIMS=(YLIMS1-YLIMMARGIN, YLIMS2+YLIMMARGIN)
    YTICKS = collect(range(YLIMS1; stop=YLIMS2, length=Int(round(num_lambdas/2;digits=0))))
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
        elseif OCCUPATION[y_i]=="EMP"
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
                        legend=:outertopright,
                        xlabel=XLABEL,
                        ylabel=YLABEL,
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
plt = create_plot(C_Ys,"Credit/Output", Rs,"Interest Rate")
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_credit_to_gdp_interest_rate.png")
plt = create_plot(C_Ys,"Credit/Output", Ws,"Wage", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_credit_to_gdp_interest_wage.png")

#lambda to: Output, Consumption, Income/Output, Income
plt = create_plot(lambdas,"λ", Outputs,"Output", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Outputs.png")
plt = create_plot(lambdas,"λ", Consumptions,"Consumption", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Consumptions.png")
plt = create_plot(lambdas,"λ", Income_to_outputs,"Income/Output", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Income_to_outputs.png")
plt = create_plot(lambdas,"λ", Incomes,"Income", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_Incomes.png")
#Credit/Output to Income/Output
plt = create_plot(C_Ys,"Credit/Output", Income_to_outputs,"Income/Output", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_credit_to_gdp_Income_to_outputs.png")

#lambda to Credit/Output, Credit/Output to Capital/Output
plt = create_plot(lambdas,"λ", C_Ys,"Credit/Output", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_lambda_credit_to_gdp.png")
plt = create_plot(C_Ys,"Credit/Output", K_Ys,"Capital/Output", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_GENERAL)$(country)_credit_to_gdp_capital_to_gdp.png")

#=
LOCAL_DIR_OCCUPATION = "$(LOCAL_DIR)/Occupation/"
if Sys.iswindows()
    LOCAL_DIR_OCCUPATION = "$(LOCAL_DIR)\\Occupation\\"
end
mkpath(LOCAL_DIR_OCCUPATION)
#Crdit/Output to share of W, SP, EMP
plts = create_plots(C_Ys,"Credit/Output", [occENTs, occSEs, occEMPs, occWs],["Share of Entrepreneurs","Share of Sole Proprietors","Share of Employers","Share of Workers"],true,["ENT","SP","EMP","W"])
display(plts[1])
savefig(plts[1],"$(LOCAL_DIR_OCCUPATION)$(country)_credit_to_gdp_entrepreneurs.png")
display(plts[2])
savefig(plts[2],"$(LOCAL_DIR_OCCUPATION)$(country)_credit_to_gdp_sole_proprietors.png")
display(plts[3])
savefig(plts[3],"$(LOCAL_DIR_OCCUPATION)$(country)_credit_to_gdp_employers.png")
display(plts[4])
savefig(plts[4],"$(LOCAL_DIR_OCCUPATION)$(country)_credit_to_gdp_workers.png")
=#
LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)/Inequality/"
if Sys.iswindows()
    LOCAL_DIR_INEQUALITY = "$(LOCAL_DIR)\\Inequality\\"
end
mkpath(LOCAL_DIR_INEQUALITY)
#Credit/Output to: log-consumption, log-earnings
plt = create_plot(C_Ys,"Credit/Output", logcs,"Variance of log-consumption", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_log_consumption.png")

plt = create_plot(C_Ys,"Credit/Output", loges,"Variance of log-earnings", false)
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_log_earnings.png")

#mean quantiles of wealth
plt = create_combined_plot(C_Ys,"Credit/Output", [quantile_mean_wealth[1,:], quantile_mean_wealth[2,:], quantile_mean_wealth[3,:], quantile_mean_wealth[4,:], quantile_mean_wealth[5,:]],["1st","2nd","3rd","4th","5th"],"Mean of Wealth (quantiles)",false,["H","W","SP","EMP","ENT"])
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_combined_credit_to_gdp_mean_wealth_quantiles.png")
plt = create_combined_plot(C_Ys,"Credit/Output", [quantile_mean_wealth_w[1,:], quantile_mean_wealth_w[2,:], quantile_mean_wealth_w[3,:], quantile_mean_wealth_w[4,:], quantile_mean_wealth_w[5,:]],["1st","2nd","3rd","4th","5th"],"Mean of Workers' Wealth (quantiles)",false,["H","W","SP","EMP","ENT"])
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_combined_credit_to_gdp_mean_wealth_w_quantiles.png")
plt = create_combined_plot(C_Ys,"Credit/Output", [quantile_mean_wealth_sp[1,:], quantile_mean_wealth_sp[2,:], quantile_mean_wealth_sp[3,:], quantile_mean_wealth_sp[4,:], quantile_mean_wealth_sp[5,:]],["1st","2nd","3rd","4th","5th"],"Mean of Sole Proprietors' Wealth (quantiles)",false,["H","W","SP","EMP","ENT"])
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_combined_credit_to_gdp_mean_wealth_sp_quantiles.png")
plt = create_combined_plot(C_Ys,"Credit/Output", [quantile_mean_wealth_emp[1,:], quantile_mean_wealth_emp[2,:], quantile_mean_wealth_emp[3,:], quantile_mean_wealth_emp[4,:], quantile_mean_wealth_emp[5,:]],["1st","2nd","3rd","4th","5th"],"Mean of Employers' Wealth (quantiles)",false,["H","W","SP","EMP","ENT"])
display(plt)
savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_combined_credit_to_gdp_mean_wealth_emp_quantiles.png")

#Crdit/Output to Gini Workers and Entrepreneurs
plts = create_plots(C_Ys,"Credit/Output", [giniWs, giniEnts],["Gini for workers' income","Gini for entrepreneurs' income"],false,["W","ENT"])
for i in 1:length(plts)
    display(plts[i])
    if i==1
        savefig(plts[i],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_gini_workers.png")
    elseif i==2
        savefig(plts[i],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_gini_entrepreneurs.png")
    end
end
#=
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
=#

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

    # calculate mean
    #means[s,h,i]
    LABELS = ["Mean of $cn' $stat_name" for cn in CHOICE_NAMES]
    OCCS = ["H","W","SP","EMP","ENT"]
    #=
    ms = [means[s,3,1:num_lambdas][1:num_lambdas_sp];  fill(NaN, num_lambdas-num_lambdas_sp)]
    plts = create_plots(C_Ys,"Credit/Output", [means[s,1,:], means[s,2,:], ms, means[s,4,:]],LABELS,false,["H","W","SP","EMP"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_mean_$(choice_name)_$(stat_name).png")
    end
    =#
    plts = create_plots(C_Ys,"Credit/Output", [means[s,1,:], means[s,2,:], means[s,3,:], means[s,4,:]],LABELS,false,["H","W","SP","EMP"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_mean_$(choice_name)_$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", means[s,h,:],LABELS[h],false,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_mean_$(choice_name)_$(stat_name).png")
    end

    # calculate gini coefficent
    #ginis[s,h,i]
    LABELS = ["Gini of $cn' $stat_name" for cn in CHOICE_NAMES]
    #=
    gs = [ginis[s,3,1:num_lambdas][1:num_lambdas_sp];  fill(NaN, num_lambdas-num_lambdas_sp)]
    plts = create_plots(C_Ys,"Credit/Output", [ginis[s,1,:], ginis[s,2,:], gs, ginis[s,4,:]],LABELS,false,["H","W","SP","EMP"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_gini_$(choice_name)_$(stat_name).png")
    end
    =#
    plts = create_plots(C_Ys,"Credit/Output", [ginis[s,1,:], ginis[s,2,:], ginis[s,3,:], ginis[s,4,:]],LABELS,false,["H","W","SP","EMP"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_gini_$(choice_name)_$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", ginis[s,h,:],LABELS[h],false,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_gini_$(choice_name)_$(stat_name).png")
    end
    plts = create_plots(C_Ys,"Credit/Output", [ginis[s,1,:], ginis[s,2,:], ginis[s,3,:], ginis[s,4,:], ginis[s,5,:]],LABELS,false,["H","W","SP","EMP","ENT"])
    for h = 5
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_gini_$(choice_name)_$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", ginis[s,h,:],LABELS[h],false,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_gini_$(choice_name)_$(stat_name).png")

    end


    # calculate variance of log-s
    #avgs[s,h,i]
    LABELS = ["Average of $cn' $stat_name" for cn in CHOICE_NAMES]
    #=
    als = [avglogs[s,3,1:num_lambdas][1:num_lambdas_sp];  fill(NaN, num_lambdas-num_lambdas_sp)]
    plts = create_plots(C_Ys,"Credit/Output", [avglogs[s,1,:], avglogs[s,2,:], als, avglogs[s,4,:]],LABELS,false,["H","W","SP","EMP"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_log$(stat_name).png")
    end
    =#
    plts = create_plots(C_Ys,"Credit/Output", [avgs[s,1,:], avgs[s,2,:], avgs[s,3,:], avgs[s,4,:]],LABELS,false,["H","W","SP","EMP"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", avgs[s,h,:],LABELS[h],false,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_$(stat_name).png")
    end
    plts = create_plots(C_Ys,"Credit/Output", [avgs[s,1,:], avgs[s,2,:], avgs[s,3,:], avgs[s,4,:], avgs[s,5,:]],LABELS,false,["H","W","SP","EMP","ENT"])
    for h in 5
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", avgs[s,h,:],LABELS[h],false,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_$(stat_name).png")
    end

    #vars[s,h,i]
    LABELS = ["Variance of $cn' $stat_name" for cn in CHOICE_NAMES]
    #=
    vls = [varlogs[s,3,1:num_lambdas][1:num_lambdas_sp];  fill(NaN, num_lambdas-num_lambdas_sp)]
    plts = create_plots(C_Ys,"Credit/Output", [varlogs[s,1,:], varlogs[s,2,:], vls, varlogs[s,4,:]],LABELS,false,["H","W","SP","EMP"], 0.0)
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_log$(stat_name).png")
    end
    =#
    plts = create_plots(C_Ys,"Credit/Output", [vars[s,1,:], vars[s,2,:], vars[s,3,:], vars[s,4,:]],LABELS,false,["H","W","SP","EMP"], 0.0)
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", vars[s,h,:],LABELS[h],false,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_$(stat_name).png")
    end
    plts = create_plots(C_Ys,"Credit/Output", [vars[s,1,:], vars[s,2,:], vars[s,3,:], vars[s,4,:], vars[s,5,:]],LABELS,false,["H","W","SP","EMP","ENT"], 0.0)
    for h in 5
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", vars[s,h,:],LABELS[h],false,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_$(stat_name).png")
    end

    # calculate variance of log-s
    #avglogs[s,h,i]
    LABELS = ["Average of $cn' Log-$stat_name" for cn in CHOICE_NAMES]
    #=
    als = [avglogs[s,3,1:num_lambdas][1:num_lambdas_sp];  fill(NaN, num_lambdas-num_lambdas_sp)]
    plts = create_plots(C_Ys,"Credit/Output", [avglogs[s,1,:], avglogs[s,2,:], als, avglogs[s,4,:]],LABELS,false,["H","W","SP","EMP"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_log$(stat_name).png")
    end
    =#
    plts = create_plots(C_Ys,"Credit/Output", [avglogs[s,1,:], avglogs[s,2,:], avglogs[s,3,:], avglogs[s,4,:]],LABELS,false,["H","W","SP","EMP"])
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_log$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", avglogs[s,h,:],LABELS[h],false,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_log$(stat_name).png")
    end
    plts = create_plots(C_Ys,"Credit/Output", [avglogs[s,1,:], avglogs[s,2,:], avglogs[s,3,:], avglogs[s,4,:], avglogs[s,5,:]],LABELS,false,["H","W","SP","EMP","ENT"])
    for h in 5
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_log$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", avglogs[s,h,:],LABELS[h],false,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_avg_$(choice_name)_log$(stat_name).png")
    end

    #varlogs[s,h,i]
    if s!=3
        LABELS = ["Variance of $cn' Log-$stat_name" for cn in CHOICE_NAMES]
    else
        LABELS = ["Variance of $cn' $stat_name" for cn in CHOICE_NAMES]
    end
    #=
    vls = [varlogs[s,3,1:num_lambdas][1:num_lambdas_sp];  fill(NaN, num_lambdas-num_lambdas_sp)]
    plts = create_plots(C_Ys,"Credit/Output", [varlogs[s,1,:], varlogs[s,2,:], vls, varlogs[s,4,:]],LABELS,false,["H","W","SP","EMP"], 0.0)
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_log$(stat_name).png")
    end
    =#
    plts = create_plots(C_Ys,"Credit/Output", [varlogs[s,1,:], varlogs[s,2,:], varlogs[s,3,:], varlogs[s,4,:]],LABELS,false,["H","W","SP","EMP"], 0.0)
    for h in 1:length(plts)
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        if h==3
            choice_name = "SoleProprietors"
        end
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_log$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", varlogs[s,h,:],LABELS[h],false,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_log$(stat_name).png")
    end
    plts = create_plots(C_Ys,"Credit/Output", [varlogs[s,1,:], varlogs[s,2,:], varlogs[s,3,:], varlogs[s,4,:], varlogs[s,5,:]],LABELS,false,["H","W","SP","EMP","ENT"], 0.0)
    for h in 5
        display(plts[h])
        choice_name = CHOICE_NAMES[h]
        savefig(plts[h],"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_log$(stat_name)_full.png")

        plt = create_plot(C_Ys,"Credit/Output", varlogs[s,h,:],LABELS[h],false,OCCS[h])
        savefig(plt,"$(LOCAL_DIR_INEQUALITY)$(country)_credit_to_gdp_var_$(choice_name)_log$(stat_name).png")
    end
end

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
# TFP_ideal for SP,EMP
plts = create_plots(C_Ys,"Credit/Output", [TFPis[2,:], TFPis[3,:]],["TFP for Sole Proprietors","TFP for Employers"],false,["SP","EMP"])
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_tfp_ideal_sp.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_tfp_ideal_emp.png")
# TFP_data for SP,EMP
plts = create_plots(C_Ys,"Credit/Output", [TFPds[2,:], TFPds[3,:]],["TFP_data for Sole Proprietors","TFP_data for Employers"],false,["SP","EMP"])
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_tfp_data_sp.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_tfp_data_emp.png")


# mean of MPL for ENT, SP,EMP
plt = create_plot(C_Ys,"Credit/Output",mean_MPL[1,:],"Mean of MPL for Entrepreneurs", false, "ENT")
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_mean_mpl_ent.png")
plts = create_plots(C_Ys,"Credit/Output", [mean_MPL[2,:], mean_MPL[3,:]],["Mean of MPL for Sole Proprietors","Mean of MPL for Employers"],false,["SP","EMP"])
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_mean_mpl_sp.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_mean_mpl_emp.png")
LABELS = ["SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [mean_MPL[2,:], mean_MPL[3,:]],LABELS,"Mean of MPL",false,["SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_credit_to_gdp_mean_mpl.png")

# variance of MPL for ENT, SP,EMP
plt = create_plot(C_Ys,"Credit/Output",var_MPL[1,:],"Variance of MPL for Entrepreneurs", false, "ENT", 0.0)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_var_mpl_ent.png")
plts = create_plots(C_Ys,"Credit/Output", [var_MPL[2,:], var_MPL[3,:]],["Variance of MPL for Sole Proprietors","Variance of MPL for Employers"],false,["SP","EMP"],0.0)
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_var_mpl_sp.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_var_mpl_emp.png")
LABELS = ["SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [var_MPL[2,:], var_MPL[3,:]],LABELS,"Variance of MPL",false,["SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_credit_to_gdp_var_mpl.png")

# mean of MPK for ENT, SP,EMP
plt = create_plot(C_Ys,"Credit/Output",mean_MPK[1,:],"Mean of MPK for Entrepreneurs", false, "ENT")
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_mean_mpk_ent.png")
plts = create_plots(C_Ys,"Credit/Output", [mean_MPK[2,:], mean_MPK[3,:]],["Mean of MPK for Sole Proprietors","Mean of MPK for Employers"],false,["SP","EMP"])
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_mean_mpk_sp.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_mean_mpk_emp.png")
LABELS = ["SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [mean_MPK[2,:], mean_MPK[3,:]],LABELS,"Mean of MPK",false,["SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_credit_to_gdp_mean_mpk.png")

# variance of MPK for ENT, SP,EMP
plt = create_plot(C_Ys,"Credit/Output",var_MPK[1,:],"Variance of MPK for Entrepreneurs", false, "ENT", 0.0)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_var_mpk_ent.png")
plts = create_plots(C_Ys,"Credit/Output", [var_MPK[2,:], var_MPK[3,:]],["Variance of MPK for Sole Proprietors","Variance of MPK for Employers"],false,["SP","EMP"],0.0)
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_var_mpk_sp.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_var_mpk_emp.png")
LABELS = ["SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [var_MPK[2,:], var_MPK[3,:]],LABELS,"Variance of MPK",false,["SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_credit_to_gdp_var_mpk.png")

# share of unbound ENT, SP,EMP
plts = create_plots(C_Ys,"Credit/Output", [share_unbound[1,:], share_unbound[2,:], share_unbound[3,:]],["Share of Unconstrained Entrepreneurs","Share of Unconstrained Sole Proprietors","Share of Unconstrained Employers"],true,["ENT","SP","EMP"])
for plt in plts
    display(plt)
end
savefig(plts[1],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_ent_unbound.png")
savefig(plts[2],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_se_unbound.png")
savefig(plts[3],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_emp_unbound.png")

# average and variance of working and managerial skills across occupations
CHOICE_NAMES = ["Workers","Entrepreneurs","Sole Proprietors","Employers"]
LABELS = ["Average Managerial Skill across $(occup)" for occup in CHOICE_NAMES]
plts = create_plots(C_Ys,"Credit/Output", [avg_m_skill[1,:], avg_m_skill[2,:], avg_m_skill[3,:], avg_m_skill[4,:]],LABELS,false,["W","ENT","SP","EMP"])
for h in 1:length(plts)
    display(plts[h])
    occup = CHOICE_NAMES[h]
    if h==3
        occup = "SoleProprietors"
    end
    savefig(plts[h],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_avg_m_skill_$(occup).png")
end
LABELS = ["W","ENT","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [avg_m_skill[1,:], avg_m_skill[2,:], avg_m_skill[3,:], avg_m_skill[4,:]],LABELS,"Average Managerial Skill",false,["W","ENT","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_avg_m_skill_with_ENT.png")
LABELS = ["W","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [avg_m_skill[1,:], avg_m_skill[3,:], avg_m_skill[4,:]],LABELS,"Average Managerial Skill",false,["W","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_avg_m_skill.png")


LABELS = ["Average Working Skill across $(occup)" for occup in CHOICE_NAMES]
plts = create_plots(C_Ys,"Credit/Output", [avg_w_skill[1,:], avg_w_skill[2,:], avg_w_skill[3,:], avg_w_skill[4,:]],LABELS,false,["W","ENT","SP","EMP"])
for h in 1:length(plts)
    display(plts[h])
    occup = CHOICE_NAMES[h]
    if h==3
        occup = "SoleProprietors"
    end
    savefig(plts[h],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_avg_w_skill_$(occup).png")
end
LABELS = ["W","ENT","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [avg_w_skill[1,:], avg_w_skill[2,:], avg_w_skill[3,:], avg_w_skill[4,:]],LABELS,"Average Working Skill",false,["W","ENT","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_avg_w_skill_with_ENT.png")
LABELS = ["W","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [avg_w_skill[1,:], avg_w_skill[3,:], avg_w_skill[4,:]],LABELS,"Average Working Skill",false,["W","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_avg_w_skill.png")


LABELS = ["Variance of Managerial Skill across $(occup)" for occup in CHOICE_NAMES]
plts = create_plots(C_Ys,"Credit/Output", [var_m_skill[1,:], var_m_skill[2,:], var_m_skill[3,:], var_m_skill[4,:]],LABELS,false,["W","ENT","SP","EMP"], 0.0)
for h in 1:length(plts)
    display(plts[h])
    occup = CHOICE_NAMES[h]
    if h==3
        occup = "SoleProprietors"
    end
    savefig(plts[h],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_var_m_skill_$(occup).png")
end
LABELS = ["W","ENT","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [var_m_skill[1,:], var_m_skill[2,:], var_m_skill[3,:], var_m_skill[4,:]],LABELS,"Variance of Managerial Skill",false,["W","ENT","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_var_m_skill_with_ENT.png")
LABELS = ["W","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [var_m_skill[1,:], var_m_skill[3,:], var_m_skill[4,:]],LABELS,"Variance of Managerial Skill",false,["W","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_var_m_skill.png")


LABELS = ["Variance of Working Skill across $(occup)" for occup in CHOICE_NAMES]
plts = create_plots(C_Ys,"Credit/Output", [var_w_skill[1,:], var_w_skill[2,:], var_w_skill[3,:], var_w_skill[4,:]],LABELS,false,["W","ENT","SP","EMP"], 0.0)
for h in 1:length(plts)
    display(plts[h])
    occup = CHOICE_NAMES[h]
    if h==3
        occup = "SoleProprietors"
    end
    savefig(plts[h],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_var_w_skill_$(occup).png")
end
LABELS = ["W","ENT","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [var_w_skill[1,:], var_w_skill[2,:], var_w_skill[3,:], var_w_skill[4,:]],LABELS,"Variance of Working Skill",false,["W","ENT","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_var_w_skill_with_ENT.png")
LABELS = ["W","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [var_w_skill[1,:], var_w_skill[3,:], var_w_skill[4,:]],LABELS,"Variance of Working Skill",false,["W","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_var_w_skill.png")

LABELS = ["Sum of Managerial Skill across $(occup)" for occup in CHOICE_NAMES]
plts = create_plots(C_Ys,"Credit/Output", [sum_m_skill[1,:], sum_m_skill[2,:], sum_m_skill[3,:], sum_m_skill[4,:]],LABELS,false,["W","ENT","SP","EMP"])
for h in 1:length(plts)
    display(plts[h])
    occup = CHOICE_NAMES[h]
    if h==3
        occup = "SoleProprietors"
    end
    savefig(plts[h],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_sum_m_skill_$(occup).png")
end
LABELS = ["W","ENT","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [sum_m_skill[1,:], sum_m_skill[2,:], sum_m_skill[3,:], sum_m_skill[4,:]],LABELS,"Sum of Managerial Skill",false,["W","ENT","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_sum_m_skill_with_ENT.png")
LABELS = ["W","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [sum_m_skill[1,:], sum_m_skill[3,:], sum_m_skill[4,:]],LABELS,"Sum of Managerial Skill",false,["W","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_sum_m_skill.png")


LABELS = ["Sum of Working Skill across $(occup)" for occup in CHOICE_NAMES]
plts = create_plots(C_Ys,"Credit/Output", [sum_w_skill[1,:], sum_w_skill[2,:], sum_w_skill[3,:], sum_w_skill[4,:]],LABELS,false,["W","ENT","SP","EMP"])
for h in 1:length(plts)
    display(plts[h])
    occup = CHOICE_NAMES[h]
    if h==3
        occup = "SoleProprietors"
    end
    savefig(plts[h],"$(LOCAL_DIR_PRODUCTIVITY)$(country)_sum_w_skill_$(occup).png")
end
LABELS = ["W","ENT","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [sum_w_skill[1,:], sum_w_skill[2,:], sum_w_skill[3,:], sum_w_skill[4,:]],LABELS,"Sum of Working Skill",false,["W","ENT","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_sum_w_skill_with_ENT.png")
LABELS = ["W","SP","EMP"]
plt = create_combined_plot(C_Ys,"Credit/Output", [sum_w_skill[1,:], sum_w_skill[3,:], sum_w_skill[4,:]],LABELS,"Sum of Working Skill",false,["W","SP","EMP"])
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_combined_sum_w_skill.png")


# output by SP
plt = create_plot(C_Ys,"Credit/Output", output_SPs,"Aggregate output by Sole Proprietors", false, "SP")
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_output_sp.png")
# output by EMP
plt = create_plot(C_Ys,"Credit/Output", output_EMPs,"Aggregate output by Employers", false, "EMP")
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_output_emp.png")

# share of output by SP
plt = create_plot(C_Ys,"Credit/Output", share_of_output_SPs,"Share of output by Sole Proprietors", true, "SP")
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_of_output_sp.png")
# share of output by EMP
plt = create_plot(C_Ys,"Credit/Output", share_of_output_EMPs,"Share of output by Employers", true, "EMP")
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_share_of_output_emp.png")

# ratio of capital restricted to unrestricted
plt = create_plot(C_Ys,"Credit/Output", ratio_of_capitals,"Ratio of Restricted Capital to Unrestricted", true)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_ratio_of_capital.png")
# ratio of output restricted to unrestricted
plt = create_plot(C_Ys,"Credit/Output", ratio_of_outputs,"Ratio of Restricted Output to Unrestricted", true)
display(plt)
savefig(plt,"$(LOCAL_DIR_PRODUCTIVITY)$(country)_credit_to_gdp_ratio_of_output.png")


#TFP_ideal, TFP_data
plt = create_plot(C_Ys,"Credit/Output", old_TFPis,"TFP", false)
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
policy = ss_star[1][4]
earnings = ss_star[1][24]

wealth = Array{Any}(undef,3)
income = Array{Any}(undef,3)
consumpation = Array{Any}(undef,3)
for occ in 1:3
    wealth[occ] = ones(size(density_distr[occ])).*asset_grid
    income[occ] = ss_star[1][23][occ] .- wealth[occ]
    consumption[occ] = income[occ] .+ wealth[occ] .- policy[occ]
end

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

for i = 1:num_of_occs # H. W, SP, EMP
    if i == 1
        mean_consumption_to_income_wealth[i] = sum( [sum((consumption[occ]./ss_star[1][23][occ]).*density_distr[occ]) for occ=1:3])
        var_consumption_to_income_wealth[i] = sum( [sum((consumption[occ]./ss_star[1][23][occ] .- mean_consumption_to_income_wealth[i]).^2 .*density_distr[occ]) for occ=1:3])
        std_consumption_to_income_wealth[i] = sqrt(var_consumption_to_income_wealth[i])


        mean_income_to_income_wealth[i] = sum( [sum((income[occ]./ss_star[1][23][occ]).*density_distr[occ]) for occ=1:3])
        var_income_to_income_wealth[i] = sum( [sum((income[occ]./ss_star[1][23][occ] .- mean_income_to_income_wealth[i]).^2 .*density_distr[occ]) for occ=1:3])
        std_income_to_income_wealth[i] = sqrt(var_income_to_income_wealth[i])


        mean_consumption_to_income[i] = sum( [sum((consumption[occ]./income[occ]).*density_distr[occ]) for occ=1:3])
        var_consumption_to_income[i] = sum( [sum((consumption[occ]./income[occ] .- mean_consumption_to_income[i]).^2 .*density_distr[occ]) for occ=1:3])
        std_consumption_to_income[i] = sqrt(var_consumption_to_income[i])

        mean_savings_to_income[i] = sum( [sum(((policy[occ].-wealth[occ])./income[occ]).*density_distr[occ]) for occ=1:3])
        var_savings_to_income[i] = sum( [sum(((policy[occ].-wealth[occ])./income[occ] .- mean_savings_to_income[i]).^2 .*density_distr[occ]) for occ=1:3])
        std_savings_to_income[i] = sqrt(var_savings_to_income[i])

        mean_consumption_to_earnings[i] = sum( [sum((consumption[occ]./earnings[occ]).*density_distr[occ]) for occ=1:3])
        var_consumption_to_earnings[i] = sum( [sum((consumption[occ]./earnings[occ] .- mean_consumption_to_earnings[i]).^2 .*density_distr[occ]) for occ=1:3])
        std_consumption_to_earnings[i] = sqrt(var_consumption_to_earnings[i])

        mean_savings_to_earnings[i] = sum( [sum(((policy[occ].-wealth[occ])./earnings[occ]).*density_distr[occ]) for occ=1:3])
        var_savings_to_earnings[i] = sum( [sum(((policy[occ].-wealth[occ])./earnings[occ] .- mean_savings_to_earnings[i]).^2 .*density_distr[occ]) for occ=1:3])
        std_savings_to_earnings[i] = sqrt(var_savings_to_earnings[i])
    else

        if i==2
            choice = 1
        elseif i==3
            choice = 2
        elseif i==4
            choice = 3
        end

        mean_consumption_to_income_wealth[i] = sum((consumption[choice]./ss_star[1][23][choice]).*density_distr[choice])
        var_consumption_to_income_wealth[i] = sum((consumption[choice]./ss_star[1][23][choice] .- mean_consumption_to_income_wealth[i]).^2 .*density_distr[choice]./sum(density_distr[choice]))
        mean_consumption_to_income_wealth[i] /= sum(density_distr[choice])
        std_consumption_to_income_wealth[i] = sqrt(var_consumption_to_income_wealth[i])


        mean_income_to_income_wealth[i] = sum((income[choice]./ss_star[1][23][choice]).*density_distr[choice])
        var_income_to_income_wealth[i] = sum((income[choice]./ss_star[1][23][choice] .- mean_income_to_income_wealth[i]).^2 .*density_distr[choice]./sum(density_distr[choice]))
        mean_income_to_income_wealth[i] /= sum(density_distr[choice])
        std_income_to_income_wealth[i] = sqrt(var_income_to_income_wealth[i])


        mean_consumption_to_income[i] = sum((consumption[choice]./income[choice]).*density_distr[choice])
        var_consumption_to_income[i] = sum((consumption[choice]./income[choice] .- mean_consumption_to_income[i]).^2 .*density_distr[choice]./sum(density_distr[choice]))
        mean_consumption_to_income[i] /= sum(density_distr[choice])
        std_consumption_to_income[i] = sqrt(var_consumption_to_income[i])

        mean_savings_to_income[i] = sum(((policy[choice].-wealth[choice])./income[choice]).*density_distr[choice])
        var_savings_to_income[i] = sum(((policy[choice].-wealth[choice])./income[choice] .- mean_savings_to_income[i]).^2 .*density_distr[choice]./sum(density_distr[choice]))
        mean_savings_to_income[i] /= sum(density_distr[choice])
        std_savings_to_income[i] = sqrt(var_savings_to_income[i])

        mean_consumption_to_earnings[i] = sum((consumption[choice]./earnings[choice]).*density_distr[choice])
        var_consumption_to_earnings[i] = sum((consumption[choice]./earnings[choice] .- mean_consumption_to_earnings[i]).^2 .*density_distr[choice]./sum(density_distr[choice]))
        mean_consumption_to_earnings[i] /= sum(density_distr[choice])
        std_consumption_to_earnings[i] = sqrt(var_consumption_to_earnings[i])

        mean_savings_to_earnings[i] = sum(((policy[choice].-wealth[choice])./earnings[choice]).*density_distr[choice])
        var_savings_to_earnings[i] = sum(((policy[choice].-wealth[choice])./earnings[choice] .- mean_savings_to_earnings[i]).^2 .*density_distr[choice]./sum(density_distr[choice]))
        mean_savings_to_earnings[i] /= sum(density_distr[choice])
        std_savings_to_earnings[i] = sqrt(var_savings_to_earnings[i])
    end
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
